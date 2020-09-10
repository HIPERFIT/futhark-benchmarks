module mandelbrot32 = import "mandelbrot32"
module mandelbrot64 = import "mandelbrot64"

import "lib/github.com/diku-dk/lys/lys"
import "lib/github.com/athas/matte/colour"

type config = { xcentre: f64
              , ycentre: f64
              , width: f64
              , limit: i32
              , radius: f64 }

let cfg (xcentre, ycentre, width, limit, radius) : config =
  {xcentre, ycentre, width, limit, radius}

let presets: []config =
  [cfg (-0.7,                   0.0,                             3.067,                  100,   16.0),
   cfg (0.20508818500545423,    0.9014915666351141   * 900/1440, 6.375321937544527e-6,   629,   256.0),
   cfg (0.4510757067879078,     0.6144133202705898   * 900/1440, 7.632248223018773e-5,   399,   4.0),
   cfg (0.3469337523117071,     0.6866350870407725   * 900/1440, 3.508380713647269e-5,   505,   1048576.0),
   cfg (-0.7902001921590814,    0.24910667566731381  * 900/1440, 5.071115028132377e-4,   1176,  3.4359738368e10),
   cfg (2.3127178455019423e-2, -1.301205470975472    * 900/1440, 3.6349313304610088e-9,  566,   4.0),
   cfg (2.3127176148480418e-2, -1.3012054707668765   * 900/1440, 2.71444790387451e-10,   604,   4.0),
   cfg (2.3127176156746785e-2, -1.301205470242045    * 900/1440, 4.49615119202067e-12,   2000,  4.0),
   cfg (0.2550376327692795,     8.962363618058007e-4 * 900/1440, 7.351698819132829e-5,   1412,  256.0),
   cfg (0.25498593633806477,    8.726424280526077e-4 * 900/1440, 1.6858526052251987e-10, 10492, 4.0)]

type text_content = (i32, i32, i32, f32)
module lys : lys with text_content = text_content = {
  type state = { config: config
               , precision: i32
               , height: i64
               , width: i64
               , mouse: (i32, i32) }

  let init _ height width : state =
    {config = presets[0], precision = 32,
     height, width, mouse = (0,0)}

  let resize height width (s: state) = s with height = height with width = width

  let diff (x1: i32, y1: i32) (x2, y2) = (x2 - x1, y2 - y1)
  let move (s: state) (dx: i32, dy: i32) =
    let aspect_ratio = f64.i64 s.width / f64.i64 s.height
    let x_per_pixel = s.config.width / f64.i64 s.width
    let height = s.config.width*(1/aspect_ratio)
    let y_per_pixel = height / f64.i64 s.height
    in s.config with xcentre = s.config.xcentre - x_per_pixel * f64.i32 dx
                with ycentre = s.config.ycentre - y_per_pixel * f64.i32 dy
  let grab_mouse = false

  let event (e: event) (s: state) =
    match e
    case #keydown {key} ->
      if key >= '0' && key <= '9'
      then s with config = presets[key-'0']
      else if key == 'p'
      then s with precision = if s.precision == 32 then 64 else 32
      else if key == 'q'
      then s with config.limit = s.config.limit - 1
      else if key == 'e'
      then s with config.limit = s.config.limit + 1
      else if key == 'z'
      then s with config.radius = s.config.radius * 0.99
      else if key == 'c'
      then s with config.radius = s.config.radius * 1.01
      else s
    case #mouse {buttons, x, y} ->
    s with mouse = (x,y) with config = if buttons != 0 then move s (diff s.mouse (x,y))
                                       else s.config
    case #wheel {dx=_, dy} ->
      s with config.width = s.config.width * (1 - 0.01 * f64.i32 dy)
    case _ -> s

  type text_content = text_content
  let text_format () = "FPS: %d; bits: %d; iterations: %d; radius: %.1f"
  let text_content (fps: f32) (s: state) = (i32.f32 (f32.round fps), s.precision,
                                            s.config.limit, f32.f64 s.config.radius)
  let text_colour _ = argb.white

  let render (s: state) =
    if s.precision == 32
    then mandelbrot32.render_mandelbrot s.width s.height
                                        s.config.xcentre s.config.ycentre
                                        s.config.width
                                        s.config.limit
                                        s.config.radius
    else mandelbrot64.render_mandelbrot s.width s.height
                                        s.config.xcentre s.config.ycentre
                                        s.config.width
                                        s.config.limit
                                        s.config.radius
}
