-- Simple seven-point stencil - in three dimensions, which is at least
-- a little interesting.
--
-- Out input data is just random garbage of the same size as the
-- Parboil data sets.  The running time is not data-sensitive here
-- anyway.
--
-- Actually, we don't even have garbage of the same size (too big for
-- Github), so we replicate a smaller slice.  And slice it for the
-- output.  Should not affect running time much.
--
-- ==
-- compiled input @ data/small.in
-- output @ data/small.out
-- compiled input @ data/default.in
-- output @ data/default.out

let next [nx][ny][nz] (c0: f32) (c1: f32) (A0: [nx][ny][nz]f32) (i: i32) (j: i32) (k: i32): f32 =
  -- We map the edge to itself.
  #[unsafe] if i == 0 || i == i32.i64 nx-1 || j == 0 || j == i32.i64 ny-1 || k == 0 || k == i32.i64 nz-1
  then A0[i,j,k]
  else (A0[i,j,k+1] + A0[i,j,k-1] +
        A0[i,j-1,k] + A0[i,j+1,k] +
        A0[i-1,j,k] + A0[i+1,j,k]) * c1
       + A0[i,j,k] * c0


let main [ny][nz] (iterations: i32) (nx: i32) (A0_slice: [ny][nz]f32): [ny][nz]f32 =
  -- c0 and c1 are also hardcoded in Rodinia.
  let c0 = 1.0f32/6.0f32
  let c1=1.0f32/6.0f32/6.0f32
  let nx = i64.i32 nx

  let A0 = loop A0 = replicate nx A0_slice for _i < iterations do
             tabulate_3d nx ny nz
                         (\i j k -> next c0 c1 A0 (i32.i64 i) (i32.i64 j) (i32.i64 k))
  in A0[0]
