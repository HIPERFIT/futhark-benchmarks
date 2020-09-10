let array f n m =
  let n64 = i64.i32 n
  let m64 = i64.i32 m
  in map (\i -> replicate m64 (f i)) (iota n64)

-- ==
-- entry: sum_iota_i32
-- only_c random input { 1 100000000 }
-- only_c random input { 10 10000000 }
-- only_c random input { 100 1000000 }
-- only_c random input { 1000 100000 }
-- only_c random input { 10000 10000 }
-- only_c random input { 100000 1000 }
-- only_c random input { 1000000 100 }
-- only_c random input { 10000000 10 }
-- only_c random input { 100000000 1 }

entry sum_iota_i32 n m: []i32 = array i32.i64 n m |> map i32.sum

-- ==
-- entry: sum_i8
-- only_c random input { 1 100000000 [100000000]i8 }
-- only_c random input { 10 10000000 [100000000]i8 }
-- only_c random input { 100 1000000 [100000000]i8 }
-- only_c random input { 1000 100000 [100000000]i8 }
-- only_c random input { 10000 10000 [100000000]i8 }
-- only_c random input { 100000 1000 [100000000]i8 }
-- only_c random input { 1000000 100 [100000000]i8 }
-- only_c random input { 10000000 10 [100000000]i8 }
-- only_c random input { 100000000 1 [100000000]i8 }

entry sum_i8 n m xs: []i8 = unflatten (i64.i32 n) (i64.i32 m) xs |> map i8.sum

-- ==
-- entry: sum_i16
-- only_c random input { 1 100000000 [100000000]i16 }
-- only_c random input { 10 10000000 [100000000]i16 }
-- only_c random input { 100 1000000 [100000000]i16 }
-- only_c random input { 1000 100000 [100000000]i16 }
-- only_c random input { 10000 10000 [100000000]i16 }
-- only_c random input { 100000 1000 [100000000]i16 }
-- only_c random input { 1000000 100 [100000000]i16 }
-- only_c random input { 10000000 10 [100000000]i16 }
-- only_c random input { 100000000 1 [100000000]i16 }

entry sum_i16 n m xs: []i16 = unflatten (i64.i32 n) (i64.i32 m) xs |> map i16.sum

-- ==
-- entry: sum_i32
-- only_c random input { 1 100000000 [100000000]i32 }
-- only_c random input { 10 10000000 [100000000]i32 }
-- only_c random input { 100 1000000 [100000000]i32 }
-- only_c random input { 1000 100000 [100000000]i32 }
-- only_c random input { 10000 10000 [100000000]i32 }
-- only_c random input { 100000 1000 [100000000]i32 }
-- only_c random input { 1000000 100 [100000000]i32 }
-- only_c random input { 10000000 10 [100000000]i32 }
-- only_c random input { 100000000 1 [100000000]i32 }

entry sum_i32 n m xs: []i32 = unflatten (i64.i32 n) (i64.i32 m) xs |> map i32.sum

-- ==
-- entry: sum_i64
-- only_c random input { 1 100000000 [100000000]i64 }
-- only_c random input { 10 10000000 [100000000]i64 }
-- only_c random input { 100 1000000 [100000000]i64 }
-- only_c random input { 1000 100000 [100000000]i64 }
-- only_c random input { 10000 10000 [100000000]i64 }
-- only_c random input { 100000 1000 [100000000]i64 }
-- only_c random input { 1000000 100 [100000000]i64 }
-- only_c random input { 10000000 10 [100000000]i64 }
-- only_c random input { 100000000 1 [100000000]i64 }

entry sum_i64 n m xs: []i64 = unflatten (i64.i32 n) (i64.i32 m) xs |> map i64.sum

-- Now for some non-commutative reductions.

-- ==
-- entry: prod_mat4_i32
-- no_python random input { 1 10000000 [10000000]i32 [10000000]i32 [10000000]i32 [10000000]i32 }
-- no_python random input { 10 1000000 [10000000]i32 [10000000]i32 [10000000]i32 [10000000]i32 }
-- no_python random input { 100 100000 [10000000]i32 [10000000]i32 [10000000]i32 [10000000]i32 }
-- no_python random input { 1000 10000 [10000000]i32 [10000000]i32 [10000000]i32 [10000000]i32 }
-- no_python random input { 10000 1000 [10000000]i32 [10000000]i32 [10000000]i32 [10000000]i32 }
-- no_python random input { 100000 100 [10000000]i32 [10000000]i32 [10000000]i32 [10000000]i32 }
-- no_python random input { 1000000 10 [10000000]i32 [10000000]i32 [10000000]i32 [10000000]i32 }
-- no_python random input { 10000000 1 [10000000]i32 [10000000]i32 [10000000]i32 [10000000]i32 }

let mat4_mul (+) (*) (a0,a1,a2,a3) (b0,b1,b2,b3) =
  ((b0*a0 + b1*a2),
   (b0*a1 + b1*a3),
   (b2*a0 + b3*a2),
   (b2*a1 + b3*a3))

let mat4' f a b c d = (f a, f b, f c, f d)

entry prod_mat4_i32 n m as bs cs ds =
  let ass = unflatten (i64.i32 n) (i64.i32 m) as
  let bss = unflatten (i64.i32 n) (i64.i32 m) bs
  let css = unflatten (i64.i32 n) (i64.i32 m) cs
  let dss = unflatten (i64.i32 n) (i64.i32 m) ds
  in map4 (map4 (mat4' id)) ass bss css dss
     |> map (reduce (mat4_mul (i32.+) (i32.*)) (1, 0, 0, 1) >-> (.0))
