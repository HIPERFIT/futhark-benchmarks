-- Benchmarks for the scan SOAC.
--
-- Some benchmarks do not perform validation for floating-point
-- datasets because the sequential reference is too inaccurate, or
-- because the random data is so badly conditioned that the results
-- are fragile (FP is not associative nor commutative, so we are
-- cheating anyway).

-- ==
-- entry: sum_iota_i8 sum_iota_i32
-- input { 10000 } auto output
-- input { 100000 } auto output
-- input { 1000000 } auto output
-- only_c input { 10000000 } auto output
-- only_c input { 100000000 } auto output

-- ==
-- entry: sum_iota_f32 sum_iota_f64
-- input { 10000 }
-- input { 100000 }
-- input { 1000000 }
-- only_c input { 10000000 }
-- only_c input { 100000000 }

entry sum_iota_i8  n = i64.i32 n |> iota |> map i8.i64 |> scan (+) 0
entry sum_iota_i32 n = i64.i32 n |> iota |> map i32.i64 |> scan (+) 0
entry sum_iota_f32 n = i64.i32 n |> iota |> map f32.i64 |> scan (+) 0
entry sum_iota_f64 n = i64.i32 n |> iota |> map f64.i64 |> scan (+) 0

-- ==
-- entry: sum_i8 sum_i32
-- random input { [10000]i32 } auto output
-- random input { [100000]i32 } auto output
-- random input { [1000000]i32 } auto output
-- only_c no_python random input { [10000000]i32 } auto output
-- only_c no_python random input { [100000000]i32 } auto output

-- ==
-- entry: sum_f32 sum_f64
-- random input { [10000]i32 }
-- random input { [100000]i32 }
-- random input { [1000000]i32 }
-- only_c no_python random input { [10000000]i32 }
-- only_c no_python random input { [100000000]i32 }

entry sum_i8  = map  i8.i32 >-> scan (+) 0
entry sum_i32 = map i32.i32 >-> scan (+) 0
entry sum_f32 = map f32.i32 >-> scan (+) 0
entry sum_f64 = map f64.i32 >-> scan (+) 0

-- ==
-- entry: sum_scaled_i8 sum_scaled_i32
-- random input { [10000]i32 } auto output
-- random input { [100000]i32 } auto output
-- random input { [1000000]i32 } auto output
-- only_c only_c random input { [10000000]i32 } auto output
-- only_c only_c random input { [100000000]i32 } auto output

-- ==
-- entry: sum_scaled_f32 sum_scaled_f64
-- random input { [10000]i32 }
-- random input { [100000]i32 }
-- random input { [1000000]i32 }
-- only_c only_c random input { [10000000]i32 }
-- only_c only_c random input { [100000000]i32 }

entry sum_scaled_i8 xs = let ys = map (i8.i32 >-> (*2)) xs
                         in (scan (+) 0 ys, ys)
entry sum_scaled_i32 xs = let ys = map (i32.i32 >-> (*2)) xs
                          in (scan (+) 0 ys, ys)
entry sum_scaled_f32 xs = let ys = map (f32.i32 >-> (*2)) xs
                          in (scan (+) 0 ys, ys)
entry sum_scaled_f64 xs = let ys = map (f64.i32 >-> (*2)) xs
                          in (scan (+) 0 ys, ys)

-- Now for some non-commutative reductions.

-- ==
-- entry: prod_iota_mat4_i8 prod_iota_mat4_i32
-- input { 10000 } auto output
-- input { 100000 } auto output
-- input { 1000000 } auto output
-- only_c input { 10000000 } auto output

-- ==
-- entry: prod_iota_mat4_f32 prod_iota_mat4_f64
-- input { 10000 }
-- input { 100000 }
-- input { 1000000 }
-- only_c input { 10000000 }

-- ==
-- entry: prod_mat4_i8 prod_mat4_i32
-- random input { [10000]i32 [10000]i32 [10000]i32 [10000]i32 } auto output
-- random input { [100000]i32 [100000]i32 [100000]i32 [100000]i32 } auto output
-- only_c random input { [1000000]i32 [1000000]i32 [1000000]i32 [1000000]i32 } auto output

-- ==
-- entry: prod_mat4_f32 prod_mat4_f64
-- random input { [10000]i32 [10000]i32 [10000]i32 [10000]i32 }
-- random input { [100000]i32 [100000]i32 [100000]i32 [100000]i32 }
-- only_c random input { [1000000]i32 [1000000]i32 [1000000]i32 [1000000]i32 }

let mat4_mul (+) (*) (a0,a1,a2,a3) (b0,b1,b2,b3) =
  ((b0*a0 + b1*a2),
   (b0*a1 + b1*a3),
   (b2*a0 + b3*a2),
   (b2*a1 + b3*a3))

let mat4 x = (x, x, x, x)
let mat4' f a b c d = (f a, f b, f c, f d)

entry prod_iota_mat4_i8  n =
 i64.i32 n |> iota |> map (i8.i64 >-> mat4)  |> scan (mat4_mul  (i8.+)  (i8.*)) (1, 0, 0, 1) |> map (.0)
entry prod_iota_mat4_i32 n =
 i64.i32 n |> iota |> map (i32.i64 >-> mat4) |> scan (mat4_mul (i32.+) (i32.*)) (1, 0, 0, 1) |> map (.0)
entry prod_iota_mat4_f32 n =
 i64.i32 n |> iota |> map (f32.i64 >-> mat4) |> scan (mat4_mul (f32.+) (f32.*)) (1, 0, 0, 1) |> map (.0)
entry prod_iota_mat4_f64 n =
 i64.i32 n |> iota |> map (f64.i64 >-> mat4) |> scan (mat4_mul (f64.+) (f64.*)) (1, 0, 0, 1) |> map (.0)

entry prod_mat4_i8 as bs cs ds  =
 map4 (mat4' i8.i32) as bs cs ds |> scan (mat4_mul  (i8.+)  (i8.*)) (1, 0, 0, 1) >-> map (.0)
entry prod_mat4_i32 as bs cs ds =
 map4 (mat4' i32.i32) as bs cs ds |> scan (mat4_mul (i32.+) (i32.*)) (1, 0, 0, 1) >-> map (.0)
entry prod_mat4_f32 as bs cs ds =
 map4 (mat4' f32.i32) as bs cs ds |> scan (mat4_mul (f32.+) (f32.*)) (1, 0, 0, 1) >-> map (.0)
entry prod_mat4_f64 as bs cs ds =
 map4 (mat4' f64.i32) as bs cs ds |> scan (mat4_mul (f64.+) (f64.*)) (1, 0, 0, 1) >-> map (.0)

-- Try a non-commutative reduction with a beefy operator containing lots of control flow.

let lss 't (t: t) (pred1: t -> bool) (pred2: t -> t -> bool) (xs: []t) =
  let max = i32.max

  let rop (lssx, lisx, lcsx, tlx, firstx, lastx)
          (lssy, lisy, lcsy, tly, firsty, lasty) =
    let connect = pred2 lastx firsty || tlx == 0 || tly == 0
    let newlss = if connect then max (lcsx + lisy) (max lssx lssy)
                            else max lssx lssy
    let newlis = if lisx == tlx && connect then lisx + lisy else lisx
    let newlcs = if lcsy == tly && connect then lcsy + lcsx else lcsy
    let first = if tlx == 0 then firsty else firstx
    let last  = if tly == 0 then lastx else lasty
    in (newlss, newlis, newlcs, tlx+tly, first, last)

    let mop x =
      let xmatch = if pred1 x then 1 else 0
      in (xmatch, xmatch, xmatch, 1, x, x)

  in (scan rop (0,0,0,0,t,t) (map mop xs)) |> map (.0)

-- ==
-- entry: lss_iota_i8 lss_iota_i32 lss_iota_f32 lss_iota_f64
-- input { 10000 } auto output
-- input { 100000 } auto output
-- input { 1000000 } auto output
-- only_c input { 10000000 } auto output

entry lss_iota_i8  n = i64.i32 n |> iota |> map  i8.i64 |> lss 0 (const true) (<=)
entry lss_iota_i32 n = i64.i32 n |> iota |> map i32.i64 |> lss 0 (const true) (<=)
entry lss_iota_f32 n = i64.i32 n |> iota |> map f32.i64 |> lss 0 (const true) (<=)
entry lss_iota_f64 n = i64.i32 n |> iota |> map f64.i64 |> lss 0 (const true) (<=)

-- ==
-- entry: lss_i8 lss_i32 lss_f32 lss_f64
-- random input { [10000]i32 } auto output
-- random input { [100000]i32 } auto output
-- random input { [1000000]i32 } auto output

entry lss_i8  = map  i8.i32 >-> lss 0 (const true) (<=)
entry lss_i32 = map i32.i32 >-> lss 0 (const true) (<=)
entry lss_f32 = map f32.i32 >-> lss 0 (const true) (<=)
entry lss_f64 = map f64.i32 >-> lss 0 (const true) (<=)
