-- BFAST-irregular: version handling obscured observations (e.g., clouds)
--
-- We do not validate the "breaks", because they appear to be numerically quite unstable.
-- ==
-- compiled input @ data/sahara-cloudy.in.gz
-- output @ data/sahara-cloudy.out.gz
-- compiled input @ data/peru.in.gz

let iota32 n = (0..1..<i32.i64 n) :> [n]i32

let logplus (x: f32) : f32 =
  if x > (f32.exp 1)
  then f32.log x else 1

let adjustValInds [N] (n : i32) (ns : i32) (Ns : i32) (val_inds : [N]i32) (ind: i32) : i32 =
    if ind < Ns - ns then (#[unsafe] val_inds[ind+ns]) - n else -1

let filterPadWithKeys [n] 't
           (p : (t -> bool))
           (dummy : t)           
           (arr : [n]t) : ([n](t,i32), i32) =
  let tfs = map (\a -> if p a then 1 else 0) arr
  let isT = scan (+) 0 tfs
  let i   = last isT
  let inds= map2 (\a iT -> if p a then i64.i32 iT-1 else -1) arr isT
  let rs  = scatter (replicate n dummy) inds arr
  let ks  = scatter (replicate n 0) inds (iota32 n)
  in  (zip rs ks, i) 

-- | builds the X matrices; first result dimensions of size 2*k+2
let mkX_with_trend [N] (k2p2: i64) (f: f32) (mappingindices: [N]i32): [k2p2][N]f32 =
  map (\i ->
        map (\ind ->
                if i == 0 then 1f32
                else if i == 1 then f32.i32 ind
                else let (i', j') = (f32.i32 (i / 2), f32.i32 ind)
                     let angle = 2f32 * f32.pi * i' * j' / f 
                     in  if i % 2 == 0 then f32.sin angle 
                                       else f32.cos angle
            ) mappingindices
      ) (iota32 k2p2)

let mkX_no_trend [N] (k2p2m1: i64) (f: f32) (mappingindices: [N]i32): [k2p2m1][N]f32 =
  map (\ i ->
        map (\ind ->
                if i == 0 then 1f32
                else let i = i + 1
		     let (i', j') = (f32.i32 (i / 2), f32.i32 ind)
                     let angle = 2f32 * f32.pi * i' * j' / f 
                     in  if i % 2 == 0 then f32.sin angle 
                                       else f32.cos angle
            ) mappingindices
      ) (iota32 k2p2m1)

---------------------------------------------------
-- Adapted matrix inversion so that it goes well --
-- with intra-blockparallelism                   --
---------------------------------------------------

  let gauss_jordan [nm] (n:i32) (m:i32) (A: *[nm]f32): [nm]f32 =
    loop A for i < n do
      let v1 = A[i]
      let A' = map (\ind -> let (k, j) = (ind / m, ind % m)
                            in if v1 == 0.0 then #[unsafe] A[k*m+j] else
                            let x = #[unsafe] (A[j] / v1) in
                                if k < n-1  -- Ap case
                                then #[unsafe] ( A[(k+1)*m+j] - A[(k+1)*m+i] * x )
                                else x      -- irow case
                   ) (iota32 nm)
      in  scatter A (iota nm) A'

  let mat_inv [n] (A: [n][n]f32): [n][n]f32 =
    let m = 2*n
    let nm = n*m
    -- Pad the matrix with the identity matrix.
    let Ap = map (\ind -> let (i, j) = (ind / m, ind % m)
                          in  if j < n then #[unsafe] ( A[i,j] )
                                       else if j == n+i
                                            then 1.0
                                            else 0.0
                 ) (iota nm)
    let Ap' = gauss_jordan (i32.i64 n) (i32.i64 m) Ap
    -- Drop the identity matrix at the front!
    in (unflatten n m Ap')[0:n,n:2*n] :> [n][n]f32
--------------------------------------------------
--------------------------------------------------

let dotprod [n] (xs: [n]f32) (ys: [n]f32): f32 =
  reduce (+) 0.0 <| map2 (*) xs ys

let matvecmul_row [n][m] (xss: [n][m]f32) (ys: [m]f32) =
  map (dotprod ys) xss

let dotprod_filt [n] (vct: [n]f32) (xs: [n]f32) (ys: [n]f32) : f32 =
  f32.sum (map3 (\v x y -> x * y * if (f32.isnan v) then 0.0 else 1.0) vct xs ys)

let matvecmul_row_filt [n][m] (xss: [n][m]f32) (ys: [m]f32) =
    map (\xs -> map2 (\x y -> if (f32.isnan y) then 0 else x*y) xs ys |> f32.sum) xss

let matmul_filt [n][p][m] (xss: [n][p]f32) (yss: [p][m]f32) (vct: [p]f32) : [n][m]f32 =
  map (\xs -> map (dotprod_filt vct xs) (transpose yss)) xss

----------------------------------------------------
----------------------------------------------------

-- | implementation is in this entry point
--   the outer map is distributed directly
entry main [m][N] (trend: i32) (k: i32) (n: i32) (freq: f32)
                  (hfrac: f32) (lam: f32)
                  (mappingindices : [N]i32)
                  (images : [m][N]f32) =
  ----------------------------------
  -- 1. make interpolation matrix --
  ----------------------------------
  let k2p2 = i64.i32 (2*k + 2)
  let k2p2' = if trend > 0 then k2p2 else k2p2-1
  let X = opaque <|
	  if trend > 0
          then mkX_with_trend k2p2' freq mappingindices
	  else mkX_no_trend   k2p2' freq mappingindices
  

  -- PERFORMANCE BUG: instead of `let Xt = copy (transpose X)`
  --   we need to write the following ugly thing to force manifestation:
  let zero = f32.i64 <| (N*N + 2*N + 1) / (N + 1) - N - 1
  let Xt  = opaque <|
            map (map (+zero)) (copy (transpose X))

  let Xh  =  (X[:,:i64.i32 n])
  let Xth =  (Xt[:i64.i32 n,:])
  let Yh  =  (images[:,:i64.i32 n])
  
  ----------------------------------
  -- 2. mat-mat multiplication    --
  ----------------------------------
  let Xsqr = opaque <|
             map (matmul_filt Xh Xth) Yh

  ----------------------------------
  -- 3. matrix inversion          --
  ----------------------------------
  let Xinv = opaque <|
             map mat_inv Xsqr
  ---------------------------------------------
  -- 4. several matrix-vector multiplication --
  ---------------------------------------------
  let beta0  = map (matvecmul_row_filt Xh) Yh   -- [2k+2]
               |> opaque

  let beta   = map2 matvecmul_row Xinv beta0    -- [2k+2]
               |> opaque -- ^ requires transposition of Xinv
                                    --   unless all parallelism is exploited

  let y_preds= map (matvecmul_row Xt) beta      -- [N]
               |> opaque -- ^ requires transposition of Xt (small)
                                    --   can be eliminated by passing
                                    --   (transpose X) instead of Xt

  ---------------------------------------------
  -- 5. filter etc.                          --
  ---------------------------------------------
  let (Nss, y_errors, val_indss) = ( opaque <| unzip3 <|
    map2 (\y y_pred ->
            let y_error_all = zip y y_pred |>
                map (\(ye,yep) -> if !(f32.isnan ye) 
                                  then ye-yep else f32.nan )
            let (tups, Ns) = filterPadWithKeys (\y -> !(f32.isnan y)) (f32.nan) y_error_all
            let (y_error, val_inds) = unzip tups
            in  (Ns, y_error, val_inds)
         ) images y_preds )

  ---------------------------------------------
  -- 6. ns and sigma                         --
  ---------------------------------------------
  let (hs, nss, sigmas) = opaque <| unzip3 <|
    map2 (\yh y_error ->
            let ns    = map (\ye -> if !(f32.isnan ye) then 1 else 0) yh
                        |> reduce (+) 0
            let sigma = map (\i -> if i < ns then #[unsafe] y_error[i] else 0.0)
                            (iota32 (i64.i32 n))
                        |> map (\ a -> a*a ) |> reduce (+) 0.0
            let sigma = f32.sqrt ( sigma / (f32.i32 (ns-i32.i64 k2p2)) )
            let h     = i32.f32 ( (f32.i32 ns) * hfrac )
            in  (h, ns, sigma)
         ) Yh y_errors

  ---------------------------------------------
  -- 7. moving sums first and bounds:        --
  ---------------------------------------------
  let hmax = reduce_comm (i32.max) 0 hs
  let MO_fsts = zip3 y_errors nss hs |>
    map (\(y_error, ns, h) -> #[unsafe]
            map (\i -> if i < h then #[unsafe] y_error[i + ns-h+1] else 0.0)
                (iota32 (i64.i32 hmax))
            |> reduce (+) 0.0 
        ) |> opaque

  let Nmn = N-i64.i32 n
  let BOUND = map (\q -> let t   = n+1+q
                         let time = #[unsafe] mappingindices[t-1]
                         let tmp = logplus ((f32.i32 time) / (f32.i32 mappingindices[N-1]))
                         in  lam * (f32.sqrt tmp)
                  ) (iota32 Nmn)

  ---------------------------------------------
  -- 8. moving sums computation:             --
  ---------------------------------------------
  let (_MOs, _MOs_NN, _breaks, means) = zip (zip4 Nss nss sigmas hs) (zip3 MO_fsts y_errors val_indss) |>
    map (\ ( (Ns,ns : i32,sigma, h), (MO_fst,y_error,val_inds) ) ->
            let MO = map (\j -> if j >= Ns-ns then 0.0
                                else if j == 0 then MO_fst
                                else #[unsafe] (-y_error[ns-h+j] + y_error[ns+j])
                         ) (iota32 Nmn) |> scan (+) 0.0
	    
            let MO' = map (\mo -> mo / (sigma * (f32.sqrt (f32.i32 ns))) ) MO
	        let (is_break, fst_break) = 
		    map3 (\mo b j ->  if j < Ns - ns && !(f32.isnan mo)
				      then ( (f32.abs mo) > b, j )
				      else ( false, j )
		         ) MO' BOUND (iota32 Nmn)
		        |> reduce_comm (\ (b1,i1) (b2,i2) -> 
               				    if b1 then (b1,i1) 
              			          else if b2 then (b2, i2)
    				              else (b1,i1) 
              	      	     ) (false, -1)
	        let mean = map2 (\x j -> if j < Ns - ns then x else 0.0 )
                                MO' (iota32 Nmn)
			    |> reduce (+) 0.0

	        let fst_break' = if !is_break then -1
                             else let adj_break = adjustValInds n ns Ns val_inds fst_break
                                  in  ((adj_break-1) / 2) * 2 + 1  -- Cosmin's validation hack
            let fst_break' = if ns <=5 || Ns-ns <= 5 then -2 else fst_break'

            let val_inds' = map (adjustValInds n ns Ns val_inds)
                                (iota32 Nmn)
            let MO'' = scatter (replicate Nmn f32.nan)
                               (map i64.i32 val_inds') MO'
            in (MO'', MO', fst_break', mean)
        ) |> unzip4

  in means
-- For Fabian: with debugging info, replace the result with the next line
--in (MO_fsts, Nss, nss, sigmas, _MOs, _MOs_NN, BOUND, breaks, means, y_errors, y_preds)


-- gcc -O2 --std=c99 bfast-cloudy-wip.c -lOpenCL -lm
-- FUTHARK_INCREMENTAL_FLATTENING=1 ~/WORK/gits/futhark/tools/futhark-autotune --compiler=futhark-opencl --pass-option=--default-tile-size=8 --stop-after 1500 --calc-timeout bfast-cloudy-wip.fut --compiler=futhark-opencl

