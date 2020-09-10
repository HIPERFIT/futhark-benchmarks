-- BFAST-irregular: version handling obscured observations (e.g., clouds)
-- ==
-- compiled input @ data/sahara.in.gz
-- output @ data/sahara.out.gz

import "lib/github.com/diku-dk/linalg/linalg"

module f32linalg = mk_linalg f32

-- f32.round
let rint (r: f32) : f32 =
  let r' = f32.i64 (i64.f32 r)
  in  if r - r' < 0.5
      then r' else r' + 1

let logplus (x: f32) : f32 =
  if x > (f32.exp 1)
  then f32.log x else 1

-- | builds the X matrices; first result dimensions of size 2*k+2
let mkX (k2p2: i64) (N: i64) (f: f32) : [k2p2][N]f32 =
  [ replicate N 1     -- first   row
  , map (f32.i64 >-> (+1)) (iota N) -- second  row
  ] ++
  ( map (\ i ->       -- sin/cos rows
          map (\j -> let i' = f32.i64 (i / 2)
                     let j' = f32.i64 j
                     let angle = 2 * f32.pi * i' * j' / f
                     in  if i % 2 == 0 then f32.sin angle else f32.cos angle
              ) (map (+1) (iota N))
        ) (2 ... k2p2-1)
  ) :> [k2p2][N]f32

-- | compute the actual number of values that precisely
--   encapsulate the first n valid values of y.
let findSplit [N] (n: i32) (y: [N]f32) : i32 =
  let flgs = map (\v -> if f32.isnan v then 0 else 1) y
  let scnf = scan (+) 0 flgs
  let pairs= map (\(i,v)->(i,v==n))
                 (zip (map i32.i64 (iota N)) scnf)
  let (m,b)= reduce (\(i1,b1) (i2,b2) ->
                        if b1 then (i1, b1)
                              else (i2, b2)
                    ) (0,false) pairs
  in  if b then m+1 else i32.i64 N

---------------------------------------------------
-- Adapted matrix inversion so that it goes well --
-- with intra-blockparallelism                   --
---------------------------------------------------

  let gauss_jordan [nm] (n:i64) (A: *[nm]f32): [nm]f32 =
    let m = nm / n in
    loop A for i < n do
      let v1 = #[unsafe] A[i]

      let A' = map (\ind -> let (k, j) = (ind / i32.i64 m, ind % i32.i64 m)
                            let x = #[unsafe] (A[j] / v1) in
                                if k < i32.i64 n-1  -- Ap case
                                then #[unsafe] ( A[(k+1)*i32.i64 m+j] -
                                                 A[(k+1)*i32.i64 m+i32.i64 i] * x )
                                else x      -- irow case
                   ) (map i32.i64 (iota (n*m)))
      in  scatter A (iota (n*m)) A'

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
    let Ap' = unflatten n m (gauss_jordan n Ap)

    -- Drop the identity matrix at the front.
    in Ap'[0:n,n:n * 2] :> [n][n]f32
--------------------------------------------------
--------------------------------------------------

let dotprod_filt [n] (flgs: [n]bool) (xs: [n]f32) (ys: [n]f32) : f32 =
    f32.sum (map3 (\flg x y -> if flg then x*y else 0) flgs xs ys)

let matsqr [m] [n] (flgs: [n]bool) (xss: [m][n]f32) : *[m][m]f32 =
    map (\xs -> map (\ys -> dotprod_filt flgs xs ys) xss) xss

let matvecmul_row_filt [n][m] (flgs: [m]bool) (xss: [n][m]f32) (ys: [m]f32) =
    map (dotprod_filt flgs ys) xss

-- | The core of the alg: the computation for a time series
--   for one pixel.
let bfast [N] (Nmn: i64) (f: f32) (k: i32) (n: i32)
              (hfrac: f32) (lam: f32)
              (y: [N]f32) :
              [Nmn]f32 =
  -- it's ok, don't panick: whatever is invariant to the
  -- outer map is gonna be hoisted out!
  -- (Just to advertize the compiler a bit)
  -- let period = rint ( (f32.i64 N) / (f32.i64 n) )
  let X = mkX (i64.i32 (2*k+2)) N f

  let m    = n -- findSplit n y   -- n
  let flgs = map (\v -> !(f32.isnan v)) y

  let flgsh = #[unsafe] (flgs[:i64.i32 m])
  let yh    = #[unsafe] (y[:i64.i32 m])
  let Xh    = #[unsafe] (X[:,:i64.i32 m])

  -- line 2, beta-hat computation
  -- fit linear regression model:
  let Xsqr = matsqr flgsh Xh              -- [2k+2][2k+2]
  let Xinv = mat_inv    Xsqr              -- [2k+2][2k+2]
  let beta0= matvecmul_row_filt flgsh Xh yh     -- [2k+2]
  let beta = f32linalg.matvecmul_row Xinv beta0 -- [2k+2]

  -- line 3, y-hat computation: get predictions and residuals.
  -- the entries in y_pred corresponding to indices `i` such
  -- that `flgs[i]==false` are invalid.
  let y_pred = f32linalg.matvecmul_row (transpose X) beta -- [N]


  -- line 4: error
  let y_error= map3 (\flg ye yep -> if flg then (ye-yep) else 0) flgs y y_pred -- [N]

  -- moving sums:
  let h = i32.f32 ( (f32.i32 m) * hfrac )
  let MO_fst = reduce (+) 0 ( #[unsafe] (y_error[i64.i32 (n-h+1) : i64.i32 (n+1)]) )
  let MO_ini = map (\j ->
                      if j == 0 then MO_fst
                      else  #[unsafe] (-y_error[n-h+j] + y_error[n+j])
                   ) (map i32.i64 (iota Nmn))
  let MO = scan (+) 0 MO_ini

  -- line 5: sigma
  let tmp = map (\ a -> a*a ) (#[unsafe] (y_error[:i64.i32 m]))

  let sigma = reduce (+) 0 tmp
  let sigma = f32.sqrt ( sigma / (f32.i32 (n-2*k-2)) )

  let MO = map (\mo -> mo / (sigma * (f32.sqrt (f32.i32 n))) ) MO

  -- line 10: BOUND computation (will be hoisted)
  let BOUND = map (\q -> let t   = n+1+q
                         let tmp = logplus ((f32.i32 t) / (f32.i32 n))
                         in  lam * (f32.sqrt tmp)
                  ) (map i32.i64 (iota Nmn))

  let breaks = map2 (\m b -> if (f32.isnan m) || (f32.isnan b) then 0 else (f32.abs m) - b) MO BOUND
  in  breaks

-- | entry point
entry main [m][N] (k: i32) (n: i32) (freq: f32)
                  (hfrac: f32) (lam: f32)
                  (images : [m][N]f32) :
                  ([m][]f32) =
  let res = map (bfast (N-i64.i32 n) freq k n hfrac lam) images
  in  res
