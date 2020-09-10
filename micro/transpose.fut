-- ==
-- entry: transpose_i8
-- only_c random input { 1 100000000 [100000000]i8 }
-- only_c random input { 2  50000000 [100000000]i8 }
-- only_c random input { 4  25000000 [100000000]i8 }
-- only_c random input { 10 10000000 [100000000]i8 }
-- only_c random input { 100 1000000 [100000000]i8 }
-- only_c random input { 1000 100000 [100000000]i8 }
-- only_c random input { 10000 10000 [100000000]i8 }
-- only_c random input { 100000 1000 [100000000]i8 }
-- only_c random input { 1000000 100 [100000000]i8 }
-- only_c random input { 10000000 10 [100000000]i8 }
-- only_c random input { 25000000  4 [100000000]i8 }
-- only_c random input { 50000000  2 [100000000]i8 }
-- only_c random input { 100000000 1 [100000000]i8 }

-- ==
-- entry: transpose_i32
-- only_c random input { 1 100000000 [100000000]i32 }
-- only_c random input { 2  50000000 [100000000]i32 }
-- only_c random input { 4  25000000 [100000000]i32 }
-- only_c random input { 10 10000000 [100000000]i32 }
-- only_c random input { 100 1000000 [100000000]i32 }
-- only_c random input { 1000 100000 [100000000]i32 }
-- only_c random input { 10000 10000 [100000000]i32 }
-- only_c random input { 100000 1000 [100000000]i32 }
-- only_c random input { 1000000 100 [100000000]i32 }
-- only_c random input { 10000000 10 [100000000]i32 }
-- only_c random input { 25000000  4 [100000000]i32 }
-- only_c random input { 50000000  2 [100000000]i32 }
-- only_c random input { 100000000 1 [100000000]i32 }

-- ==
-- entry: transpose_i64
-- only_c random input { 1 100000000 [100000000]i64 }
-- only_c random input { 2  50000000 [100000000]i64 }
-- only_c random input { 4  25000000 [100000000]i64 }
-- only_c random input { 10 10000000 [100000000]i64 }
-- only_c random input { 100 1000000 [100000000]i64 }
-- only_c random input { 1000 100000 [100000000]i64 }
-- only_c random input { 10000 10000 [100000000]i64 }
-- only_c random input { 100000 1000 [100000000]i64 }
-- only_c random input { 1000000 100 [100000000]i64 }
-- only_c random input { 10000000 10 [100000000]i64 }
-- only_c random input { 25000000  4 [100000000]i64 }
-- only_c random input { 50000000  2 [100000000]i64 }
-- only_c random input { 100000000 1 [100000000]i64 }

let unpack_2d n m = unflatten (i64.i32 n) (i64.i32 m)

entry transpose_i8  n m xs: [][]i8  = unpack_2d n m xs |> transpose
entry transpose_i32 n m xs: [][]i32 = unpack_2d n m xs |> transpose
entry transpose_i64 n m xs: [][]i64 = unpack_2d n m xs |> transpose

-- ==
-- entry: map_transpose_i8
-- only_c random input { 1 1 100000000 [100000000]i8 }
-- only_c random input { 10 1 10000000 [100000000]i8 }
-- only_c random input { 1000 1 100000 [100000000]i8 }
-- only_c random input { 1 1000 100000 [100000000]i8 }
-- only_c random input { 10 1000 10000 [100000000]i8 }
-- only_c random input { 1000 1000 100 [100000000]i8 }
-- only_c random input { 1 10000000 10 [100000000]i8 }
-- only_c random input { 10 10000000 1 [100000000]i8 }
-- only_c random input { 1000 100000 1 [100000000]i8 }

-- ==
-- entry: map_transpose_i32
-- only_c random input { 1 1 100000000 [100000000]i32 }
-- only_c random input { 10 1 10000000 [100000000]i32 }
-- only_c random input { 1000 1 100000 [100000000]i32 }
-- only_c random input { 1 1000 100000 [100000000]i32 }
-- only_c random input { 10 1000 10000 [100000000]i32 }
-- only_c random input { 1000 1000 100 [100000000]i32 }
-- only_c random input { 1 10000000 10 [100000000]i32 }
-- only_c random input { 10 10000000 1 [100000000]i32 }
-- only_c random input { 1000 100000 1 [100000000]i32 }

-- ==
-- entry: map_transpose_i64
-- only_c random input { 1 1 100000000 [100000000]i64 }
-- only_c random input { 10 1 10000000 [100000000]i64 }
-- only_c random input { 1000 1 100000 [100000000]i64 }
-- only_c random input { 1 1000 100000 [100000000]i64 }
-- only_c random input { 10 1000 10000 [100000000]i64 }
-- only_c random input { 1000 1000 100 [100000000]i64 }
-- only_c random input { 1 10000000 10 [100000000]i64 }
-- only_c random input { 10 10000000 1 [100000000]i64 }
-- only_c random input { 1000 100000 1 [100000000]i64 }

let unpack_3d k n m = unflatten_3d (i64.i32 k) (i64.i32 n) (i64.i32 m)

entry map_transpose_i8  k n m xs:  [][][]i8 = unpack_3d k n m xs |> map transpose
entry map_transpose_i32 k n m xs: [][][]i32 = unpack_3d k n m xs |> map transpose
entry map_transpose_i64 k n m xs: [][][]i64 = unpack_3d k n m xs |> map transpose
