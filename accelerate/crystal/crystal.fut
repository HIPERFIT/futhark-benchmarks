-- The first data set is equal to the Accelerate default.
--
-- ==
-- tags { futhark-c futhark-opencl }
-- bench input {
--   200
--   30.0f32
--   5
--   1
--   1.0f32
-- }
--
-- nobench compiled input {
--   20
--   30.0f32
--   5
--   50
--   0.5f32
-- }
--
-- nobench compiled input {
--   40
--   30.0f32
--   5
--   50
--   0.5f32
-- }
--
-- nobench compiled input {
--   40
--   30.0f32
--   50
--   50
--   0.5f32
-- }
--
-- input {
--   2000
--   30.0f32
--   50
--   1
--   1.0f32
-- }
--
-- input {
--   4000
--   30.0f32
--   50
--   1
--   1.0f32
-- }

import "lib/github.com/athas/matte/colour"

let odd(n: i32): bool = (n & 1) == 1

let point(scale: f32, x: f32, y: f32): (f32, f32) =
  (x * scale, y * scale)

let rampColour(v: f32): argb.colour =
  argb.from_rgba 1.0 (0.4 + (v * 0.6)) v 0.0

let wrap(n: f32): f32 =
  let n' = n - f32.i64(i64.f32(n))
  let odd_in_int = i64.f32(n) & 1
  let even_in_int = 1 - odd_in_int
  in f32.i64(odd_in_int) * (1.0 - n') + f32.i64(even_in_int) * n'

let wave(th: f32, x: f32, y: f32): f32 =
  let cth = f32.cos(th)
  let sth = f32.sin(th)
  in (f32.cos(cth * x + sth * y) + 1.0) / 2.0

let waver(th: f32, x: f32, y: f32, n: i64): f32 =
  f32.sum (tabulate n (\i -> wave(f32.i64(i) * th, x, y)))

let waves(degree: i64, phi: f32, x: f32, y: f32): f32 =
  let th = f32.pi / phi
  in wrap(waver(th, x, y, degree))

let quasicrystal(scale: f32, degree: i64, time: f32, x: f32, y: f32): argb.colour =
  let phi = 1.0 + (time ** 1.5) * 0.005
  let (x', y') = point(scale, x, y)
  in rampColour(waves(degree, phi, x', y'))

let normalize_index(i: i64, field_size: i64): f32 =
  f32.i64 (i) / f32.i64 (field_size)

entry render_frame (field_size: i64) (scale: f32) (degree: i64) (time: f32)
                  : [field_size][field_size]argb.colour =
  tabulate_2d field_size field_size
  (\y x -> quasicrystal(scale, degree, time,
                        normalize_index(x, field_size),
                        normalize_index(y, field_size)))

entry main(field_size: i32) (scale: f32) (degree: i32)
          (n_steps: i32) (time_delta: f32) =
  let field_size = i64.i32 field_size
  let degree = i64.i32 degree
  -- Hack to avoid returning something gigantic.
  let frames = tabulate (i64.i32 n_steps) (\step_i ->
                      let time = f32.i64 step_i * time_delta
                      in render_frame field_size scale degree time)
  in map (\frame -> [i32.u32 (frame[0,0] % u32.i64 (length (flatten frame)))])
         frames
