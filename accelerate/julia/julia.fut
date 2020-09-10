import "lib/github.com/diku-dk/complex/complex"
import "lib/github.com/athas/matte/colour"

module c32 = mk_complex f32
type c32 = c32.complex

let dot ((x,y): c32) = x*x + y*y

let julia (screenX: i64) (screenY: i64) (next: f32 -> c32 -> c32)
          (t: f32) (x0: f32) (y0: f32) (width: f32)
          (limit: i32) (radius: f32)
        : [screenY][screenX](c32, i32) =
  let complex_of_pixel x y =
    let height = f32.i64 screenY / f32.i64 screenX * width
    let xmin   = x0 - width  / 2
    let ymin   = y0 - height / 2

    let re     = xmin + (f32.i64 x * width)  / f32.i64 screenX
    let im     = ymin + (f32.i64 y * height) / f32.i64 screenY
    in c32.mk re im

  let pixel y x = let z0 = complex_of_pixel x y
                  in loop (z, i) = (z0, 0) while i < limit
                                                 && dot z < radius
                     do (next t z, i + 1)
  in tabulate_2d screenY screenX pixel

let escape_to_rgba (limit: i32) (palette: []argb.colour) ((z, n): (c32, i32)): argb.colour =
  if n == limit then argb.black
  else let scale   = 256
       let shift   = 1664
       let mag     = c32.mag z
       let smooth  = f32.log2 (f32.log2 mag)
       let ix      = i64.f32 (f32.sqrt (f32.i32 n + 1 - smooth) * scale + shift) % length palette
    in #[unsafe] palette[ix]

let render (palette: []argb.colour)
           (screenX: i64) (screenY: i64)
           (next: f32 -> c32 -> c32) (t: f32)
           (x0: f32) (y0: f32) (width: f32)
           (limit: i32) (radius: f32)
         : [screenY][screenX]argb.colour =
  julia screenX screenY next t x0 y0 width limit radius
  |> map (map (escape_to_rgba limit palette))

import "lib/github.com/diku-dk/lys/lys"

-- | Cubic interpolation.
let cubic (x0:f32,x1:f32) (y0:f32,y1:f32) (m0:f32,m1:f32) (x:f32) =
  let h    = x1 - x0
  let t    = (x - x0) / h
  let h_00 = (1.0 + 2.0*t) * (1.0 - t) ** 2.0
  let h_10 = t * (1.0 - t) ** 2.0
  let h_01 = t ** 2.0 * (3.0 - 2.0 * t)
  let h_11 = t ** 2.0 * (t - 1.0)
  in y0 * h_00 + h * m0 * h_10 + y1 * h_01 + h * m1 * h_11

let interp (x0:f32, x1:f32)
           (y0:argb.colour, y1:argb.colour)
           ((mr0,mg0,mb0):(f32,f32,f32),
            (mr1,mg1,mb1):(f32,f32,f32))
           (x: f32) =
  let (r0,g0,b0,_) = argb.to_rgba y0
  let (r1,g1,b1,_) = argb.to_rgba y1
  in argb.from_rgba (cubic (x0,x1) (r0,r1) (mr0,mr1) x)
                    (cubic (x0,x1) (g0,g1) (mg0,mg1) x)
                    (cubic (x0,x1) (b0,b1) (mb0,mb1) x)
                    0.0

let ultra (p: f32) =
  let p0 = 0.0
  let p1 = 0.16
  let p2 = 0.42
  let p3 = 0.6425
  let p4 = 0.8575
  let p5 = 1.0

  let rgb8 (r: i32) (g: i32) (b: i32) =
    argb.from_rgba (f32.i32 r / 255.0) (f32.i32 g / 255.0) (f32.i32 b / 255.0) 0.0

  let c0 = rgb8 0   7   100
  let c1 = rgb8 32  107 203
  let c2 = rgb8 237 255 255
  let c3 = rgb8 255 170 0
  let c4 = rgb8 0   2   0
  let c5 = c0

  let m0 = (0.7843138, 2.4509804,  2.52451)
  let m1 = (1.93816,   2.341629,   1.6544118)
  let m2 = (1.7046283, 0.0,        0.0)
  let m3 = (0.0,       -2.2812111, 0.0)
  let m4 = (0.0,       0.0,        0.0)
  let m5 = m0

  in
  if p <= p1 then interp (p0,p1) (c0,c1) (m0,m1) p else
  if p <= p2 then interp (p1,p2) (c1,c2) (m1,m2) p else
  if p <= p3 then interp (p2,p3) (c2,c3) (m2,m3) p else
  if p <= p4 then interp (p3,p4) (c3,c4) (m3,m4) p else
                  interp (p4,p5) (c4,c5) (m4,m5) p


let ultra_palette (points: i64) : [points]argb.colour =
  tabulate points (\ix -> ultra (f32.i64 ix / f32.i64 points))

let golden = (1 + f32.sqrt 5) / 2

-- All the 'next' functions are of the form
--
--   f(t, z) = z ** 2 + a * cos t + b * sin t
--
-- for various 'a' and 'b'.
type next_fn = {a: f32, b: f32}

let next ({a,b}: next_fn) (t: f32) (z: c32) =
  z c32.* z c32.+ c32.mk (a * f32.cos t) (b * f32.sin t)

type preset = (next_fn, f32, f32, f32, i32, f32)
let presets : []preset =
  map (\(a,b) -> ({a,b}, 0, 0, 4, 244, 15))
      ([(0.7885, 0.7885),
        (1-golden, 0),
        (0.285, 0.1),
        (0.45, 0.1428),
        (-0.70176, -0.3842),
        (-0.835, -0.2321),
        (-0.8, 0.156),
        (-0.7269, 0.1889),
        (0, -0.8)])

module lys : lys with text_content = () = {
  type state =
    { screenX: i64
    , screenY: i64
    , palette: [2048]argb.colour
    , time: f32
    , speed: f32
    , posX: f32
    , posY: f32
    , width: f32
    , radius: f32
    , iters: i32
    , next_fn: next_fn
    , animating: bool
    , zooming: f32
    , panning: (f32, f32)
    }

  let load_preset (s: state) (p: preset) =
    let (f, posX, posY, width, iters, radius) = p
    in s with posX = posX
         with posY = posY
         with next_fn = f
         with width = width
         with iters = iters
         with radius = radius

  let init _ h w : state =
    let s = { screenX = w
            , screenY = h
            , palette = ultra_palette 2048
            , time = 0
            , speed = 1
            , posX = 0
            , posY = 0
            , width = 4
            , radius = 16
            , iters = 255
            , next_fn = {a=0, b=0}
            , animating = true
            , zooming = 0
            , panning = (0, 0)
            }
    in load_preset s presets[0]

  let event (e: event) (s: state) =
    match e

    case #keydown {key} ->
      if key >= SDLK_1 && key <= SDLK_9 then
        load_preset s presets[key-SDLK_1]
      else if key == SDLK_z then s with radius = s.radius * 0.5
      else if key == SDLK_c then s with radius = s.radius * 2
      else if key == SDLK_a then s with iters = i32.f32 (f32.i32 s.iters * 0.8)
      else if key == SDLK_d then s with iters = i32.f32 (f32.i32 s.iters * 1.2)
      else if key == SDLK_q then s with speed = s.speed * 0.75
      else if key == SDLK_e then s with speed = s.speed * 1.25

      else if key == SDLK_w then s with zooming = -0.25
      else if key == SDLK_s then s with zooming = 0.25

      else if key == SDLK_RIGHT then s with time = s.time + s.speed / 60
      else if key == SDLK_LEFT then s with time = s.time - s.speed / 60

      else if key == SDLK_SPACE then s with animating = !s.animating
      else s

    case #keyup {key} ->
      if key == SDLK_w || key == SDLK_s then s with zooming = 0
      else s

    case #mouse {buttons, x, y} ->
      let dx = (s.panning.0-f32.i32 x) * s.width / f32.i64 s.screenX
      let dy = (s.panning.1-f32.i32 y) * s.width / f32.i64 s.screenY
      let (x', y') = if buttons != 0
                     then (dx + s.posX, dy + s.posY)
                     else (s.posX, s.posY)
      in s with panning = (f32.i32 x, f32.i32 y)
           with posX = x'
           with posY = y'

    case #step td ->
      s with time = s.time + f32.bool s.animating * td * s.speed
        with width = (1 + td * s.zooming) * s.width

    case _ -> s

  let resize h w (s: state) =
    s with screenX = w with screenY = h

  let render (s: state) =
    render s.palette s.screenX s.screenY (next s.next_fn) s.time s.posX s.posY
           s.width s.iters s.radius

  type text_content = ()
  let grab_mouse = false
  let text_format () = ""
  let text_content _ _ = ()
  let text_colour _ = argb.white
}
