import "lib/github.com/athas/vector/vspace"
import "physics"
module bruteforce = import "nbody"
module bh = import "nbody-bh"

type rotation = ((f32, f32, f32),
                 (f32, f32, f32),
                 (f32, f32, f32))

let rotate_point (rotation: rotation) ({x,y,z}: position): position =
  {x= x*rotation.0.0 + y*rotation.1.0 + z*rotation.2.0,
   y= x*rotation.0.1 + y*rotation.1.1 + z*rotation.2.1,
   z= x*rotation.0.2 + y*rotation.1.2 + z*rotation.2.2}

let rotate_points [n] (rotation: rotation) (ps: [n]position): [n]position =
  map (rotate_point rotation) ps

let rotate_x_matrix (angle: f32): rotation =
  ((1f32,        0f32,         0f32),
   (0f32, f32.cos angle, -f32.sin angle),
   (0f32, f32.sin angle,  f32.cos angle))

let rotate_y_matrix (angle: f32): rotation =
  ((f32.cos angle,  0f32, f32.sin angle),
   (0f32,           1f32, 0f32),
   (-f32.sin angle, 0f32, f32.cos angle))

let rotmult (x: rotation) (y: rotation): rotation =
  let dot a b = a.0 * b.0 + a.1 * b.1 + a.2 * b.2
  in ((dot (x.0.0, x.1.0, x.2.0) y.0,
       dot (x.0.1, x.1.1, x.2.1) y.0,
       dot (x.0.2, x.1.2, x.2.2) y.0),
      (dot (x.0.0, x.1.0, x.2.0) y.1,
       dot (x.0.1, x.1.1, x.2.1) y.1,
       dot (x.0.2, x.1.2, x.2.2) y.1),
      (dot (x.0.0, x.1.0, x.2.0) y.2,
       dot (x.0.1, x.1.1, x.2.1) y.2,
       dot (x.0.2, x.1.2, x.2.2) y.2))

let rotation (x_rotation: f32) (y_rotation: f32): rotation =
  rotmult (rotate_x_matrix x_rotation) (rotate_y_matrix y_rotation)

let inverse_rotation (x_rotation: f32) (y_rotation: f32): rotation =
  rotmult (rotate_y_matrix y_rotation) (rotate_x_matrix x_rotation)

import "lib/github.com/athas/matte/colour"

-- FIXME: This rendering is terrible because it is not FoV-aware.
let render_point (h: i64) (w: i64)
                 {x=x_ul: f32, y=y_ul: f32} {x=x_br: f32, y=y_br: f32}
                 (max_mass: f32) ({x,y,z=_}:position) (m: f32):
                 (i64, argb.colour) =
  -- Draw nothing if the point is outside the viewport.
  if x < x_ul || x > x_br || y < y_ul || y > y_br then (-1, 0)
  else
  -- Normalise x,y to positions in interval (0,1) within the viewport.
  let x' = (x-x_ul) / (x_br-x_ul)
  let y' = (y-y_ul) / (y_br-y_ul)
  -- Convert x',y' to screen coordinate space.
  let x'' = i64.f32(x' * f32.i64(w))
  let y'' = i64.f32(y' * f32.i64(h))
  let intensity = if m >= max_mass
                  then 255
                  else 128 + u32.f32((m / max_mass) * 128)
  let colour = intensity * 0x10000 +
               intensity * 0x100 +
               0xFF
  in (y''*w + x'', colour)

import "lib/github.com/diku-dk/lys/lys"
import "lib/github.com/diku-dk/cpprandom/random"

module rnge = pcg32
module dist = uniform_real_distribution f32 rnge
module norm = normal_distribution f32 rnge

let mk_body_cloud (radius: f32) (rng: rnge.rng) : (rnge.rng, body) =
  let (rng, position) =
    let (rng, angle) = dist.rand (0, 2*f32.pi) rng
    let (rng, z_angle) = dist.rand (0, 2*f32.pi) rng
    let (rng, length) = dist.rand (0, radius) rng
    let x = f32.cos angle * length
    let y = f32.sin angle * length
    let z = f32.sin z_angle * length
    in (rng, {x,y,z})
  let (rng, mass) = dist.rand (0, 1) rng
  let velocity = {x=0, y=0, z=0}
  in (rng, {position, mass, velocity})

let mk_body_spiral (radius: f32) (rng: rnge.rng) : (rnge.rng, body) =
  let (rng, position) =
    let (rng, length) = dist.rand (0, radius) rng
    let angle = length / 3*f32.pi
    let x = f32.sin angle * length
    let y = f32.cos angle * length
    let (rng, z) = dist.rand (0, radius/10) rng
    in (rng, {x,y,z})
  let (rng, mass) = dist.rand (0, 1) rng
  let velocity = {x=0, y=0, z=0}
  in (rng, {position, mass, velocity})

let mk_body_orbit (radius: f32) (rng: rnge.rng) : (rnge.rng, body) =
  let (rng, position) =
    let (rng, x) = norm.rand {mean=0, stddev=radius} rng
    let (rng, y) = norm.rand {mean=0, stddev=radius/100} rng
    let (rng, z) = norm.rand {mean=0, stddev=radius/100} rng
    let y = y * (x/100)
    let z = z * (x/10)
    in (rng, {x,y,z})
  let (rng, mass) = dist.rand (0, 1) rng
  let velocity = {x=0, y=position.x/radius*10, z=0}
  in (rng, {position, mass, velocity})

let mk_bodies (mk: (rnge.rng -> (rnge.rng, body))) (rng: rnge.rng) (n: i64)
            : (rnge.rng, [n]body) =
  let rngs = rnge.split_rng n rng
  let (rngs, bodies) = unzip (map mk rngs)
  in (rnge.join_rng rngs, bodies)

type opt 'a = #none | #some a
type solver = #bruteforce | #bh
type info = (f32, i32, i32, f32)
module lys : lys with text_content = info = {
  type~ state = { bodies: []body
                , height: i64
                , width: i64
                , ul: {x: f32, y: f32}
                , br: {x: f32, y: f32}
                , invert: bool
                , max_mass: f32
                , rotation: {x: f32, y: f32}
                , rotating: {x: f32, y: f32}
                , zooming: f32
                , rng : rnge.rng
                , paused : bool
                , solver: solver
                , theta: f32
                , attractor : opt vec3.vector
                }

  let epsilon : f32 = 50

  let screen_space_to_world_space x y (s: state) =
    let x_dist = s.br.x - s.ul.x
    let y_dist = s.br.y - s.ul.y
    let x = s.ul.x + (f32.i64 x / f32.i64 s.width) * x_dist
    let y = s.ul.y + (f32.i64 y / f32.i64 s.height) * y_dist
    in rotate_point (inverse_rotation (-s.rotation.x) (-s.rotation.y)) {x,y,z=0}

  let blob x y (s: state) : state =
    let {x,y,z} = screen_space_to_world_space x y s
    let blob_N = 100
    let move (body: body) = body with position.x = body.position.x + x
                                 with position.y = body.position.y + y
                                 with position.z = body.position.z + z
    let (rng, bodies) = mk_bodies (mk_body_cloud 3) s.rng blob_N
    in s with rng = rng
         with bodies = s.bodies ++ map move bodies

  let maybe_attract [n] (opt_attractor: opt vec3.vector) (td: f32) (bodies: [n]body)
                      : [n]body =
    match opt_attractor
    case #none -> bodies
    case #some position ->
      let pm = {position, mass = 100000}
      let accels = map (\b -> accel epsilon (pointmass b) pm) bodies
      in map2 (advance_body td) bodies accels

  let init (seed: u32) (height: i64) (width: i64) : state =
    let n = 10000
    let rng = rnge.rng_from_seed [i32.u32 seed]
    let (rng, bodies) = mk_bodies (mk_body_cloud (f32.min (f32.i64 height) (f32.i64 width))) rng n
    let max_mass = f32.maximum (map (.mass) bodies)
    in { bodies
       , height, width
       , ul = {x = f32.i64 width / -2,
               y = f32.i64 height / -2}
       , br = {x = f32.i64 width / 2,
               y = f32.i64 height / 2}
       , invert = false
       , max_mass
       , rotation = {x=0, y=0}
       , rotating = {x=0, y=0}
       , zooming = 0
       , rng
       , paused = false
       , solver = #bruteforce
       , theta = 0.5
       , attractor = #none
       }

  let move_bodies td (s: state) =
    s with bodies = if s.paused
                    then s.bodies
                    else maybe_attract s.attractor td
                         (match s.solver
                          case #bruteforce ->
                            bruteforce.advance_bodies epsilon 1 s.bodies
                          case #bh ->
                            bh.advance_bodies s.theta epsilon 1 s.bodies)

  let handle_rotation td (s: state) =
    s with rotation.x = s.rotation.x + td * s.rotating.x
      with rotation.y = s.rotation.y + td * s.rotating.y

  let handle_zoom td (s: state) =
    let x_dist = s.br.x - s.ul.x
    let y_dist = s.br.y - s.ul.y
    in s with br = {x = s.br.x - td * x_dist * s.zooming,
                    y = s.br.y - td * y_dist * s.zooming}
         with ul = {x = s.ul.x + td * x_dist * s.zooming,
                    y = s.ul.y + td * y_dist * s.zooming}

  let event (e: event) (s: state) =
    match e

    case #step td -> s
                     |> move_bodies td
                     |> handle_rotation td
                     |> handle_zoom td

    case #keydown {key} ->
      if key == SDLK_RIGHT then s with rotating.y = 1
      else if key == SDLK_LEFT then s with rotating.y = -1
      else if key == SDLK_UP then s with rotating.x = 1
      else if key == SDLK_DOWN then s with rotating.x = -1
      else if key == SDLK_SPACE then s with paused = !s.paused
      else if key == SDLK_PAGEUP then s with zooming = 1
      else if key == SDLK_PAGEDOWN then s with zooming = -1
      else if key == SDLK_m then s with solver = match s.solver
                                                 case #bruteforce -> #bh
                                                 case #bh -> #bruteforce
      else if key == SDLK_1 then s with theta = f32.max 0 (s.theta * 0.9)
      else if key == SDLK_2 then s with theta = f32.min 100 (s.theta * 1.1)

      else if key == SDLK_s
      then let radius =
             f32.min (f32.i64 s.height) (f32.i64 s.width) / 2
           let (rng, bodies) =
             mk_bodies (mk_body_spiral radius) s.rng (length s.bodies)
           in s with rng = rng with bodies = bodies

      else if key == SDLK_o
      then let radius = f32.min (f32.i64 s.height) (f32.i64 s.width) / 2
           let (rng, bodies) =
             mk_bodies (mk_body_orbit radius) s.rng (length s.bodies)
           in s with rng = rng with bodies = bodies

      else if key == SDLK_r
      then let radius = f32.min (f32.i64 s.height) (f32.i64 s.width) / 2
           let (rng, bodies) =
             mk_bodies (mk_body_cloud radius) s.rng (length s.bodies)
           in s with rng = rng with bodies = bodies

      else s

    case #keyup {key} ->
      if key == SDLK_RIGHT || key == SDLK_LEFT then s with rotating.y = 0
      else if key == SDLK_UP || key == SDLK_DOWN then s with rotating.x = 0
      else if key == SDLK_PAGEUP || key == SDLK_PAGEDOWN then s with zooming = 0
      else s

    case #mouse {buttons, x, y} ->
      if (buttons & 1) != 0
      then blob (i64.i32 x) (i64.i32 y) s
      else if (buttons & 4) != 0
      then s with attractor = #some (screen_space_to_world_space (i64.i32 x) (i64.i32 y) s)
      else s with attractor = #none

    case _ -> s

  let resize h w (s: state) =
    s with height = h with width = w

  let render (s: state) =
    let background = if s.invert then argb.white else argb.black
    let (is, vs) =
      unzip (map2
              (render_point s.height s.width s.ul s.br s.max_mass)
              (rotate_points (rotation s.rotation.x s.rotation.y)
                             (map (.position) s.bodies))
              (map (.mass) s.bodies))
    let vs' = map (\x -> if s.invert then !x else x) vs
    in unflatten s.height s.width
       (scatter (replicate (s.height*s.width) background) is vs')

  type text_content = info
  let text_content (fps: f32) (s: state) : info =
    (fps,
     i32.i64 (length s.bodies),
     match s.solver case #bruteforce -> 0
                    case #bh -> 1,
     s.theta)
  let grab_mouse = false
  let text_format () = "FPS: %f\nN: %d\nSolver: %[brute force|Barnes-Hut]\nTheta: %f"
  let text_colour (s: state) = if s.invert then argb.black else argb.white
}
