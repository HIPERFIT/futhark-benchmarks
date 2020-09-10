-- Simple game of life implementation with a donut world.
--
-- ==
-- nobench compiled input {
--   [[0, 0, 0, 0, 0],
--    [0, 0, 1, 0, 0],
--    [0, 0, 0, 1, 0],
--    [0, 1, 1, 1, 0],
--    [0, 0, 0, 0, 0]]
--   10000
--   4i64
--   }

import "genlife"
import "conway"
import "quad"
import "quad2"
import "rule101"

entry conway_init [n][m] (bs: [n][m]bool) =
 conway.init bs
entry conway_steps [n][m] (k: i32) (world: [n][m]conway.cell) =
 loop (world) for _i < k do conway.step world
entry conway_render [n][m] (world: [n][m]conway.cell) =
 conway.render world
entry conway_uninit [n][m] (world: [n][m]conway.cell) =
 conway.uninit world

type fading_world [n][m] = [n][m]conway_fading.cell
entry conway_fading_init [n][m] (bs: [n][m]bool): fading_world [n][m] =
 conway_fading.init bs
entry conway_fading_steps (k: i32) (world: fading_world [][]): fading_world [][] =
 loop (world) for _i < k do conway_fading.step world
entry conway_fading_render (world: fading_world [][]) =
 conway_fading.render world
entry conway_fading_uninit (world: fading_world [][]) =
 conway_fading.uninit world

type quad_world [n][m] = [n][m]quad.cell
entry quad_init [n][m] (bs: [n][m]bool): quad_world [n][m] =
 quad.init bs
entry quad_steps (k: i32) (world: quad_world [][]): quad_world [][] =
 loop (world) for _i < k do quad.step world
entry quad_render (world: quad_world [][]) =
 quad.render world
entry quad_uninit (world: quad_world [][]) =
 quad.uninit world

type quad_fading_world [n][m] = [n][m]quad_fading.cell
entry quad_fading_init [n][m] (bs: [n][m]bool): quad_fading_world [n][m] =
 quad_fading.init bs
entry quad_fading_steps (k: i32) (world: quad_fading_world [][]): quad_fading_world [][] =
 loop (world) for _i < k do quad_fading.step world
entry quad_fading_render (world: quad_fading_world [][]) =
 quad_fading.render world
entry quad_fading_uninit (world: quad_fading_world [][]) =
 quad_fading.uninit world

type quad2_world [n][m] = [n][m]quad2.cell
entry quad2_init [n][m] (bs: [n][m]bool): quad2_world [n][m] =
 quad2.init bs
entry quad2_steps (k: i32) (world: quad2_world [][]): quad2_world [][] =
 loop (world) for _i < k do quad2.step world
entry quad2_render (world: quad2_world [][]) =
 quad2.render world
entry quad2_uninit (world: quad2_world [][]) =
 quad2.uninit world

type rule101_world [n][m] = [n][m]rule101.cell
entry rule101_init [n][m] (bs: [n][m]bool): rule101_world [n][m] =
 rule101.init bs
entry rule101_steps (k: i32) (world: rule101_world [][]): rule101_world [][] =
 loop (world) for _i < k do rule101.step world
entry rule101_render (world: rule101_world [][]) =
 rule101.render world
entry rule101_uninit (world: rule101_world [][]) =
 rule101.uninit world

-- Just a simple test to force the program to be compiled.
let main [n][m] (base_pattern: [n][m]i32) (_repeats: i32) (k: i64) =
  let l = k * m
  let pattern = unflatten (n*k) (m*k)
                (flatten_3d (replicate k
                             (map (\row -> flatten_to l (replicate k (map bool.i32 row)))
                              base_pattern)))
  let world =
    loop world = conway.init pattern for _i < k do conway.step world
  in i32.sum (map i32.bool (flatten world))
