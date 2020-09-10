-- Parallel blocked LU-decomposition.
--
-- ==
-- input @ data/16by16.in
-- output @ data/16by16.out
-- compiled input @ data/64.in
-- output @ data/64.out
-- compiled input @ data/256.in
-- output @ data/256.out
-- compiled input @ data/512.in.gz
-- output @ data/512.out.gz
-- compiled input @ data/2048.in.gz
-- output @ data/2048.out.gz


let dotprod [n] (a: [n]f32) (b: [n]f32): f32 =
  map2 (*) a b
       |> reduce (+) 0


-------------------------------------------
---- lud_diagonal_omp -----
--
-- This version should be tuned and run with FUTHARK_INCREMENTAL_FLATTENING=1
-------------------------------------------
---------------------------------------------------------------------
---------------------------------------------------------------------
let lud_diagonal [b] (a: [b][b]f32): *[b][b]f32 =
  map1 (\mat ->
          let mat = copy mat
          in loop (mat: *[b][b]f32) for i < b-1 do
             let col = map (\j -> if j > i then
                                    #[unsafe] (mat[j,i] - (dotprod mat[j,:i] mat[:i,i])) / mat[i,i]
                                  else
                                    mat[j,i])
                           (iota b)
            let mat[:,i] = col

            let row = map (\j -> if j > i then
                                   mat[i+1, j] - (dotprod mat[:i+1, j] mat[i+1, :i+1])
                                 else
                                   mat[i+1, j])
                          (iota b)
            let mat[i+1] = row

            in mat
       ) (unflatten (opaque 1) b a)
       |> head

------------------------------
------------------------------
---- LUD Perimeter Upper -----
------------------------------
------------------------------

let lud_perimeter_upper [m][b] (diag: [b][b]f32, a0s: [m][b][b]f32): *[m][b][b]f32 =
    let a1s = map (\ (x: [b][b]f32): [b][b]f32  -> transpose(x)) a0s in
    let a2s =
        map  (\a1: [b][b]f32  ->
              map  (\row0: [b]f32  ->   -- Upper
                    loop row=replicate b 0.0f32 for i < b do
                    let sum = (loop sum=0.0f32 for k < i do sum + diag[i,k] * row[k])
                    let row[i] = row0[i] - sum
                    in  row
                   ) a1
             ) a1s
    in map (\x: [b][b]f32 -> transpose(x)) a2s


------------------------------
------------------------------
---- LUD Perimeter Lower -----
------------------------------
------------------------------

let lud_perimeter_lower [b][m] (diag: [b][b]f32, mat: [m][b][b]f32): *[m][b][b]f32 =
  map (\blk: [b][b]f32  ->
        map  (\ (row0: [b]f32): *[b]f32  ->   -- Lower
                loop row=replicate b 0.0f32 for j < b do
                        let sum = loop sum=0.0f32 for k < j do
                            sum + diag[k,j] * row[k]
                        let row[j] = (row0[j] - sum) / diag[j,j]
                        in  row
            ) blk
      ) mat


------------------------------
------------------------------
----     LUD Internal    -----
------------------------------
------------------------------

let lud_internal [m][b] (top_per: [m][b][b]f32, lft_per: [m][b][b]f32, mat_slice: [m][m][b][b]f32 ): *[m][m][b][b]f32 =
  let top_slice = map transpose top_per in
  map (\(mat_arr: [m][b][b]f32, lft: [b][b]f32): [m][b][b]f32  ->
        map (\ (mat_blk: [b][b]f32, top: [b][b]f32): [b][b]f32  ->
                map  (\ (mat_row: [b]f32, lft_row: [b]f32): [b]f32  ->
                        map  (\(mat_el, top_row)  ->
                                let prods = map2 (*) lft_row top_row
                                let sum   = f32.sum prods
                                in mat_el - sum
                             ) (zip (mat_row) top)
                    ) (zip (mat_blk) lft )
           ) (zip (mat_arr) (top_slice) )
     ) (zip (mat_slice) (lft_per) )

let block_size: i64 = 32

let pad_to [n] 'a (m: i64) (x: a) (arr: [n]a) : [m]a =
  arr ++ replicate (m - n) x :> [m]a

--------------------------------------------
---- Main Driver:
--------------------------------------------
let main [m] (mat: [m][m]f32): [m][m]f32 =
    let b = block_size
    let num_blocks = (m+b-1) / b -- rounding up
    let n = b * num_blocks
    -- Maybe pad the input to be a multiple of the block size.
    let padding = n - m
    let mat = if padding != 0
              then map (pad_to n 0) mat ++
                   replicate padding (replicate n 0f32)
              else mat :> [n][n]f32
    -------------------------------------------------
    ---- transform matrix in [n/b,n/b,b,b] block ----
    ---- versions for upper and lower parts      ----
    ---- the blocks of the lower part            ----
    -------------------------------------------------
    let matb =
        map  (\i_b: [num_blocks][b][b]f32  ->
                map  (\j_b: [b][b]f32  ->
                        map (\i: [b]f32  ->
                                map  (\j: f32  ->
                                        #[unsafe] mat[i_b*b+i, j_b*b + j]
                                    ) (iota(b) )
                           ) (iota(b) )
                    ) (iota(num_blocks) )
            ) (iota(num_blocks) )

    --------------------------------------
    ---- sequential tiled loop driver ----
    --------------------------------------
    let matb = loop(matb) for step < ((n / b) - 1) do
        -----------------------------------------------
        ---- 1. compute the current diagonal block ----
        -----------------------------------------------
        let diag = lud_diagonal(matb[step,step]) in

        ----------------------------------------
        ---- 2. compute the top  perimeter  ----
        ----------------------------------------
        let row_slice = matb[step,step+1:num_blocks]
        let top_per_irreg = lud_perimeter_upper(diag, row_slice)
        
        ----------------------------------------
        ---- 3. compute the left perimeter  ----
        ----    and update matrix           ----
        ----------------------------------------
        let col_slice = matb[step+1:num_blocks,step]
        let lft_per_irreg = lud_perimeter_lower(diag, col_slice)
        
        ----------------------------------------
        ---- 4. compute the internal blocks ----
        ----------------------------------------
        let inner_slice = matb[step+1:num_blocks,step+1:num_blocks]
        let internal = lud_internal(top_per_irreg, lft_per_irreg, inner_slice)

        ----------------------------------------
        ---- 5. update matrix in place      ----
        ----------------------------------------
        let matb[step,step] = diag
        let matb[step, step+1:num_blocks] = top_per_irreg
        let matb[step+1:num_blocks, step] = lft_per_irreg
        let matb[step+1:num_blocks, step+1:num_blocks] = internal
        
        in matb
    ---------------------
    -- LOOP ENDS HERE! --
    ---------------------

    let last_step = (n / b) - 1 in
    let matb[last_step,last_step] = 
            lud_diagonal( matb[last_step, last_step] )

    let ret_padded = map (\i_ind  ->
                          map  (\j_ind  ->
                                let (ii, jj) = (i_ind/b, j_ind/b)
                                let ( i,  j) = (i_ind - ii*b, j_ind - jj*b)
                                in  #[unsafe] matb[ii,jj,i,j]
                               ) (iota n)
                         ) (iota n)
    in take m (map (take m) ret_padded)
