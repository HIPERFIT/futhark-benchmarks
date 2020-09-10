let replicate32 n = replicate (i64.i32 n)

let bins k is = map (%k) is |> map i64.i32

-- Simple cases with addition operator, which can be translated
-- directly to atomic addition.
-- ==
-- entry: sum_i32
-- random input { 10 [1000000]i32 [1000000]i32 } auto output
-- random input { 100 [1000000]i32 [1000000]i32 } auto output
-- random input { 1000 [1000000]i32 [1000000]i32 } auto output
-- random input { 10000 [1000000]i32 [1000000]i32 } auto output
-- random input { 100000 [1000000]i32 [1000000]i32 } auto output

entry sum_i32 [n] (k: i32) (is : [n]i32) (vs : [n]i32) : []i32 =
  reduce_by_index (replicate32 k 0) (+) 0 (bins k is) vs

-- An f32 requires a little more work from the compiler.
-- ==
-- entry: sum_f32
-- random input { 10 [1000000]i32 [1000000]f32 } auto output
-- random input { 100 [1000000]i32 [1000000]f32 } auto output
-- random input { 1000 [1000000]i32 [1000000]f32 } auto output
-- random input { 10000 [1000000]i32 [1000000]f32 } auto output
-- random input { 100000 [1000000]i32 [1000000]f32 } auto output

entry sum_f32 [n] (k: i32) (is : [n]i32) (vs : [n]f32) : []f32 =
  reduce_by_index (replicate32 k 0) (+) 0 (bins k is) vs

-- Do both!
-- ==
-- entry: sum_i32_f32
-- random input { 10 [1000000]i32 [1000000]i32 [1000000]f32 } auto output
-- random input { 100 [1000000]i32 [1000000]i32 [1000000]f32 } auto output
-- random input { 1000 [1000000]i32 [1000000]i32 [1000000]f32 } auto output
-- random input { 10000 [1000000]i32 [1000000]i32 [1000000]f32 } auto output
-- random input { 100000 [1000000]i32 [1000000]i32 [1000000]f32 } auto output

entry sum_i32_f32 [n] (k: i32) (is : [n]i32) (vs1 : [n]i32) (vs2 : [n]f32) : ([]i32, []f32) =
  (reduce_by_index (replicate32 k 0) (+) 0 (bins k is) vs1,
   reduce_by_index (replicate32 k 0) (+) 0 (bins k is) vs2)

-- Now a fancier operator, but because the payload is an i32, an
-- efficient implementation is possible.
-- ==
-- entry: absmax_i32
-- random input { 10 [1000000]i32 [1000000]i32 } auto output
-- random input { 100 [1000000]i32 [1000000]i32 } auto output
-- random input { 1000 [1000000]i32 [1000000]i32 } auto output
-- random input { 10000 [1000000]i32 [1000000]i32 } auto output
-- random input { 100000 [1000000]i32 [1000000]i32 } auto output

let absmax (x: i32) (y: i32): i32 =
  if i32.abs x < i32.abs y then y else x

entry absmax_i32 [n] (k: i32) (is : [n]i32) (vs : [n]i32) : []i32 =
  reduce_by_index (replicate32 k 0) absmax 0 (bins k is) vs

-- Now a vectorised operator.  If the compiler is clever, it can
-- compile this quite efficiently.
-- ==
-- entry: sum_vec_i32
-- random input { 10 [10000]i32 [1000000]i32 } auto output
-- random input { 10 [1000]i32 [1000000]i32 } auto output
-- random input { 10000 [10000]i32 [1000000]i32 } auto output
-- random input { 10000 [1000]i32 [1000000]i32 } auto output

entry sum_vec_i32 [n][m] (k: i32) (is : [m]i32) (vs : [n]i32) : [][]i32 =
  let l = n/m
  let vs' = unflatten m l vs
  in reduce_by_index (replicate32 k (replicate l 0)) (map2 (+))
                     (replicate l 0) (bins k is) vs'

-- An operator that the compiler really cannot do anything clever
-- about - a locking-based approach is needed.
-- -- ==
-- entry: argmax_i32
-- random input { 10 [1000000]i32 [1000000]i32 } auto output
-- random input { 100 [1000000]i32 [1000000]i32 } auto output
-- random input { 1000 [1000000]i32 [1000000]i32 } auto output
-- random input { 10000 [1000000]i32 [1000000]i32 } auto output
-- random input { 100000 [1000000]i32 [1000000]i32 } auto output

let argmax_op ((x: i32), (i: i32)) ((y: i32), (j: i32)) =
  if y > x then (y, j)
  else if y < x then (x, i)
  else if i > j then (y, j)
  else (x, i)

entry argmax_i32 [n] (k: i32) (is : [n]i32) (vs : [n]i32) : ([]i32, []i32) =
  reduce_by_index (replicate32 k (i32.lowest, -1))
                  argmax_op (i32.lowest, -1)
                  (bins k is) (zip vs (map i32.i64 (iota n)))
  |> unzip
