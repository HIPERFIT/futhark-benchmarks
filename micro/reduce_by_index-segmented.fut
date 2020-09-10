let bins k is = map (%k) is |> map i64.i32

-- Histogram-like computation.
-- ==
-- entry: sum_i32
-- random input { 10 [256][4000]i32 [256][4000]i32 } auto output
-- random input { 1000 [256][4000]i32 [256][4000]i32 } auto output
-- random input { 100000 [256][4000]i32 [256][4000]i32 } auto output

entry sum_i32 [n][m] (k: i32) (iss: [n][m]i32) (vss: [n][m]i32) : [][]i32 =
  let k64 = i64.i32 k
  in map2 (\is vs -> reduce_by_index (replicate k64 0) (+) 0 (bins k is) vs) iss vss
