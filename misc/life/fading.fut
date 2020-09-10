-- Add history-based fading to any Game of Life ruleset.  Only affects
-- the visualisation; the underlying rules are the same.
--
-- When a cell goes "dead", it does not immediately switch to the
-- colour of dead cells.  Instead, it will switch to the "dying"
-- colour, and fade to the dead colour over some period of steps.

import "genlife"

import "lib/github.com/athas/matte/colour"

module type fading = {
  val dying_speed: f32
}

module fading_life (F: fading) (R: rules) (V: visuals with cell = R.cell)
       : vis_game_of_life =
{
-- The cell, and the amount of time it has been dead.
type cell = (R.cell, argb.colour, i32)

let dead (c: R.cell): bool = V.uninit c == false

  open (gen_life_vis {
type cell = cell
let value ((c,_,_): cell) = R.value c
let weights = R.weights
let step ((c,col,h): cell) (neighbours: i32) =
  let c' = R.step c neighbours
  let died = ! (dead c) && dead c'
  in (c',
      if died then V.colour c else col,
      if died then 0          else h + 1)
} {
  type cell = cell

let init (b: bool): cell = (V.init b, argb.black, 10000)

let uninit ((c,_,_): cell) = V.uninit c

let colour ((c,col,h): cell) =
  let normal = V.colour c in
  if dead c
  then argb.mix 1f32 col (f32.i32 h * F.dying_speed) normal
  else normal
})

}

module slow_fader = fading_life { let dying_speed = 0.1f32 }
