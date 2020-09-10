let array f n m =
  let n64 = i64.i32 n
  let m64 = i64.i32 m
  in map (\i -> replicate m64 (f i)) (iota n64)

-- ==
-- entry: sum_iota_i32
-- only_c input { 1 10000000 } auto output
-- only_c input { 10 1000000 } auto output
-- only_c input { 100 100000 } auto output
-- only_c input { 1000 10000 } auto output
-- only_c input { 10000 1000 } auto output
-- only_c input { 100000 100 } auto output
-- only_c input { 1000000 10 } auto output
-- only_c input { 10000000 1 } auto output

entry sum_iota_i32 n m = array i32.i64 n m |> map (scan (+) 0)

-- ==
-- entry: sum_i32
-- only_c random input { [1][10000000]i32 } auto output
-- only_c random input { [10][1000000]i32 } auto output
-- only_c random input { [100][100000]i32 } auto output
-- only_c random input { [1000][10000]i32 } auto output
-- only_c random input { [10000][1000]i32 } auto output
-- only_c random input { [100000][100]i32 } auto output
-- only_c random input { [1000000][10]i32 } auto output
-- only_c random input { [10000000][1]i32 } auto output

entry sum_i32 = map (scan (+) 0i32)
