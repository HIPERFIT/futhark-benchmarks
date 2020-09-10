-- Fluid simulation library.
--
-- A port of Accelerate's version:
-- https://github.com/AccelerateHS/accelerate-examples/tree/master/examples/fluid
-- (mostly based on the C version, since it was simpler).
--
-- Variable naming conventions:
--
--   + `g`: The grid resolution.  The fluid simulation works on a square grid of
--          size `g * g`.
--   + `i`: the horizontal index.
--   + `j`: the vertical index.
--   + `S*`: Some array, purpose unknown.
--   + `U*`: Array of horizontal forces.
--   + `V*`: Array of vertical forces.
--   + `D*`: Array of densities.
--
-- The simulation handles the 1-size border around the `(g - 2) * (g - 2)` grid
-- differently from the inner content.  The border depends on the outermost
-- inner values.  The original C version writes this border after writing the
-- inner values, by reading from the array.  We would like to handle this
-- pattern with a single map, so the Futhark port instead *calculates* the
-- outermost inner values when it needs them for the outer bound values, which
-- means that a few calculations are done twice.  The alternative would be to
-- first calculate all the inner values, and then write the outer values
-- afterwards.
--
-- ==
-- compiled input @ benchmarking/medium.in

-- medium.in attributes: n_steps=1, n_solver_steps=40, grid_res=100

------------------------------------------------------------
-- General helper functions.
------------------------------------------------------------

let clamp (x: f32): i8 =
  if x < 0.0
  then 0i8
  else if x > 255.0
  then 255i8
  else i8.f32 x

let inside
  (i: i32)
  (j: i32)
  (g: i32):
  bool =
  i >= 1 && i <= g - 2
  && j >= 1 && j <= g - 2

let in_outside_corner
  (i: i32)
  (j: i32)
  (g: i32):
  bool =
  (i == 0 || i == g - 1) && (j == 0 || j == g - 1)

module type edge_handling_mapper = {
  type info [g]

  val inner: i32 -> i32 -> (g: i64) -> info [g] -> f32
}

module edge_handling (mapper: edge_handling_mapper) = {
  let corner_index_neighbors
    (i: i32) (j: i32) (g: i32): (i32, i32, i32, i32) =
    if i == 0 && j == 0
    then (1, 0, 0, 1)
    else if i == 0 && j == g - 1
    then (1, g - 1, 0, g - 2)
    else if i == g - 1 && j == 0
    then (g - 2, 0, g - 1, 1)
    else (g - 2, g - 1, g - 1, g - 2)

  let outermost_inner_index
    (i: i32) (j: i32) (g: i32) (b: i32): (i32, i32, f32) =
    if i == 0
    then (1, j, if b == 1 then -1.0 else 1.0)
    else if i == g - 1
    then (g - 2, j, if b == 1 then -1.0 else 1.0)
    else if j == 0
    then (i, 1, if b == 2 then -1.0 else 1.0)
    else if j == g - 1
    then (i, g - 2, if b == 2 then -1.0 else 1.0)
    else (0, 0, 0.0) -- This is not supposed to happen.

  let handle
    (i: i32) (j: i32) (g: i64) (b: i32)
    (info: mapper.info []): f32 =
    if inside i j (i32.i64 g)
    then mapper.inner i j g info
    else let base (i': i32) (j': i32): f32 =
           let (i1, j1, f) = outermost_inner_index i' j' (i32.i64 g) b
           in f * mapper.inner i1 j1 g info

         in if in_outside_corner i j (i32.i64 g)
            then let (i1, j1, i2, j2) = corner_index_neighbors i j (i32.i64 g)
                 in 0.5 * (base i1 j1 + base i2 j2)
            else base i j
}


------------------------------------------------------------
-- Implementation.
------------------------------------------------------------

module edge_handling_lin_solve = edge_handling({
  type info [g] = ([g][g]f32, [g][g]f32, f32, f32)

  let inner
    (i: i32)
    (j: i32)
    (g: i64)
    ((s0, s1, a, c): info [g]):
    f32 =
    -- A stencil.
    #[unsafe] ((s0[i, j] + a *
             (  s1[(i - 1), j]
              + s1[(i + 1), j]
              + s1[i, (j - 1)]
              + s1[i, (j + 1)])) / c)

})

let lin_solve [g]
  (n_solver_steps: i32)
  (s0: [g][g]f32)
  (b: i32)
  (a: f32)
  (c: f32):
  [g][g]f32 =
  loop s1 = replicate g (replicate g 0.0) for _k < n_solver_steps do
    map (\i ->
         map (\j ->
              edge_handling_lin_solve.handle (i32.i64 i) (i32.i64 j) g b (s0, s1, a, c))
             (0..<g))
        (0..<g)


let diffuse [g]
  (s: [g][g]f32)
  (b: i32)
  (n_solver_steps: i32)
  (diffusion_rate_or_viscosity: f32)
  (time_step: f32):
  [g][g]f32 =
  let a = (time_step * diffusion_rate_or_viscosity
           * f32.i64 (g - 2) * f32.i64 (g - 2))
  in lin_solve n_solver_steps s b a (1.0 + 4.0 * a)


module edge_handling_advect = edge_handling({
  type info [g] = ([g][g]f32, [g][g]f32, [g][g]f32, f32)

  let inner
    (i: i32)
    (j: i32)
    (g: i64)
    ((s0, u, v, time_step0): info [g]):
    f32 =
    let x = f32.i32 i - time_step0 * #[unsafe] u[i, j]
    let y = f32.i32 j - time_step0 * #[unsafe] v[i, j]

    let x = if x < 0.5 then 0.5 else x
    let x = if x > f32.i64 (g - 2) + 0.5 then f32.i64 (g - 2) + 0.5 else x
    let i0 = i32.f32 x
    let i1 = i0 + 1

    let y = if y < 0.5 then 0.5 else y
    let y = if y > f32.i64 (g - 2) + 0.5 then f32.i64 (g - 2) + 0.5 else y
    let j0 = i32.f32 y
    let j1 = j0 + 1

    let s1 = x - f32.i32 i0
    let s0' = 1.0 - s1
    let t1 = y - f32.i32 j0
    let t0 = 1.0 - t1

    in #[unsafe] (s0' * (t0 * s0[i0, j0] + t1 * s0[i0, j1])
               + s1 * (t0 * s0[i1, j0] + t1 * s0[i1, j1]))
})

let advect [g]
  (s0: [g][g]f32)
  (u: [g][g]f32)
  (v: [g][g]f32)
  (b: i32)
  (time_step: f32):
  *[g][g]f32 =

  let time_step0 = time_step * f32.i64 (g - 2)
  in map (\i -> map (\j ->
                     edge_handling_advect.handle (i32.i64 i) (i32.i64 j) g b (s0, u, v, time_step0))
                    (0..<g))
         (0..<g)


module edge_handling_project_top = edge_handling({
  type info [g] = ([g][g]f32, [g][g]f32)

  let inner
    (i: i32)
    (j: i32)
    (g: i64)
    ((u0, v0): info [g]):
    f32 =
    #[unsafe] (-0.5 * (  u0[i + 1, j]
                    - u0[i - 1, j]
                    + v0[i, j + 1]
                    - v0[i, j - 1]) / f32.i64 g)
})

module edge_handling_project_bottom = edge_handling({
  type info [g] = ([g][g]f32, [g][g]f32, i32, i32, i32, i32)

  let inner
    (i: i32)
    (j: i32)
    (g: i64)
    ((p0, s0, i0d, j0d, i1d, j1d): info [g]):
    f32 =
    #[unsafe] (s0[i, j] - 0.5 * f32.i64 (g - 2)
            * (p0[i + i0d, j + j0d] - p0[i + i1d, j + j1d]))
})

let project [g]
  (n_solver_steps: i32)
  (u0: [g][g]f32)
  (v0: [g][g]f32):
  (*[g][g]f32, *[g][g]f32) =

  let project_top: [g][g]f32 =
    map (\i -> map (\j ->
                     edge_handling_project_top.handle (i32.i64 i) (i32.i64 j) g 0 (u0, v0))
                   (0..<g))
        (0..<g)

  let project_bottom
    (p0: [g][g]f32)
    (s0: [g][g]f32)
    (b: i32)
    (i0d: i32)
    (j0d: i32)
    (i1d: i32)
    (j1d: i32):
    *[g][g]f32 =
    map (\i -> map (\j ->
                    edge_handling_project_bottom.handle (i32.i64 i) (i32.i64 j) g b
                      (p0, s0, i0d, j0d, i1d, j1d))
                   (0..<g))
        (0..<g)

  let div0 = project_top
  let p0 = lin_solve n_solver_steps div0 0 1.0 4.0
  let u1 = project_bottom p0 u0 1 1 0 (-1) 0
  let v1 = project_bottom p0 v0 2 0 1 0 (-1)
  in (u1, v1)


------------------------------------------------------------
-- Step function.
------------------------------------------------------------

let step [g]
  (u0: [g][g]f32)
  (v0: [g][g]f32)
  (d0: [g][g]f32)
  (n_solver_steps: i32)
  (time_step: f32)
  (diffusion_rate: f32)
  (viscosity: f32):
  (*[g][g]f32, *[g][g]f32, *[g][g]f32) =

  let (u1, v1) =
    let u1 = diffuse u0 1 n_solver_steps viscosity time_step
    let v1 = diffuse v0 2 n_solver_steps viscosity time_step
    let (u2, v2) = project n_solver_steps u1 v1
    let u3 = advect u2 u2 v2 1 time_step
    let v3 = advect v2 u2 v2 2 time_step
    let (u4, v4) = project n_solver_steps u3 v3
    in (u4, v4)

  let dens_step =
    let d1 = diffuse d0 0 n_solver_steps diffusion_rate time_step
    let d2 = advect d1 u0 v0 0 time_step
    in d2

  let d1 = dens_step
  in (u1, v1, d1)


------------------------------------------------------------
-- Visualisation.
------------------------------------------------------------

let draw_densities [g]
  (ds: [g][g]f32)
  (g_minus_two: i64):
  [g_minus_two][g_minus_two][3]i8 =
  let ks = 1..2...g_minus_two
  in map (\i: [g_minus_two][3]i8  ->
            map (\j: [3]i8  ->
                   let value = clamp (255.0 * #[unsafe] ds[i, j])
                   in [value, value, value]) ks) ks

let draw_one_frame [g]
  (u0: [g][g]f32)
  (v0: [g][g]f32)
  (d0: [g][g]f32)
  (n_solver_steps: i32)
  (time_step: f32)
  (diffusion_rate: f32)
  (viscosity: f32)
  (g_minus_two: i64):
  ([g_minus_two][g_minus_two][3]i8,
   [g][g]f32, [g][g]f32, [g][g]f32) =
  let (u1, v1, d1) = step u0 v0 d0 n_solver_steps
                          time_step diffusion_rate viscosity
  in (draw_densities d1 g_minus_two, u1, v1, d1)

entry draw_one_frame_raw [g]
  (u0: [g][g]f32)
  (v0: [g][g]f32)
  (d0: [g][g]f32)
  (n_solver_steps: i32)
  (time_step: f32)
  (diffusion_rate: f32)
  (viscosity: f32): ([][][3]i8,
                     [g][g]f32,
                     [g][g]f32,
                     [g][g]f32) =
  draw_one_frame u0 v0 d0 n_solver_steps time_step diffusion_rate
                 viscosity (g - 2)

------------------------------------------------------------
-- Benchmarking.
------------------------------------------------------------

let get_end_frame [g]
  (u0: [g][g]f32)
  (v0: [g][g]f32)
  (d0: [g][g]f32)
  (n_steps: i32)
  (n_solver_steps: i32)
  (time_step: f32)
  (diffusion_rate: f32)
  (viscosity: f32):
  ([g][g]f32, [g][g]f32, [g][g]f32) =
  loop (u0, v0, d0) for _i < n_steps do
    step u0 v0 d0 n_solver_steps time_step
         diffusion_rate viscosity

let main [g]
  (u0: [g][g]f32)
  (v0: [g][g]f32)
  (d0: [g][g]f32)
  (n_steps: i32)
  (n_solver_steps: i32)
  (time_step: f32)
  (diffusion_rate: f32)
  (viscosity: f32): ([g][g]f32,
                    [g][g]f32,
                    [g][g]f32) =
  get_end_frame u0 v0 d0 n_steps n_solver_steps
                time_step diffusion_rate viscosity

entry poke [g] (a: *[g][g]f32) (i: i32) (j: i32) (v: f32) = a with [i,j] = v
