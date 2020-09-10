-- A wrapper around the Pan library that provides actual screen
-- frames.

import "lib/github.com/athas/matte/colour"
import "pan"

let fcolorToColour (r,g,b,a) = argb.from_rgba r g b a

let boolToColour: bool -> argb.colour = boolToFColor >-> fcolorToColour

type^ argb_image = img argb.colour

let cimage_to_argb: cimage -> argb_image = (>-> fcolorToColour)

let region_to_argb: region -> argb_image = (>-> boolToColour)

-- Fractal definition from "Composing Fractals" by Mark P. Jones.  The
-- main changes are due to the fact that Futhark is not lazy.
let fairlyClose ((u,v): point) = (u*u + v*v) < 100f32

import "lib/github.com/diku-dk/complex/complex"

module complex = mk_complex f32
type complex = complex.complex

let mandelbrotNext p1 p1' = complex.(p1 + p1' * p1')

let mandelbrotEscapes (n: i32): img (complex, i32) =
  \p1 -> let continue (p, i) = i < n && fairlyClose p
         let next (p1', i) = (mandelbrotNext p1 p1', i+1)
         in iterate_while continue next (p1, 0)

let juliaEscapes (p0: point) (n: i32): img (complex, i32) =
  \p1 -> let continue (p, i) = i < n && fairlyClose p
         let next (p1, i) = (mandelbrotNext p0 p1, i+1)
         in iterate_while continue next (p1, 0)

let mandelbrotSet (n: i32): region =
  mandelbrotEscapes n >-> (\(_, i) -> i==n)

let mandelbrotImage 'color (n: i32) (colourise: (complex, i32) -> color): img color =
  mandelbrotEscapes n >-> colourise

let mandelbrotGreyscale (n: i32): cimage =
  mandelbrotImage n <| \(_, i) -> let i' = f32.i32 i / f32.i32 n
                                  in (i', i', i', 1f32)

let visualise_argb_image (img: argb_image)
                         (screen_width: i64) (screen_height: i64)
                         (width: f32) (height: f32)
                         (xcentre: f32) (ycentre: f32) =
  let aspect_ratio = width / height
  -- Project physical pixel coordinate to position in plane.
  let (xmin,ymin) = ((xcentre - width/2f32),
                     (ycentre - (1f32/aspect_ratio)*width/2f32))
  let (xmax,ymax) = ((xcentre + width/2f32),
                     (ycentre + (1f32/aspect_ratio)*width/2f32))
  let sizex = xmax - xmin
  let sizey = ymax - ymin
  let p (y,x) = (xmin + (f32.i64 x * sizex) / f32.i64 screen_width,
                 ymin + (f32.i64 y * sizey) / f32.i64 screen_height)
  in tabulate_2d screen_height screen_width (curry (p >-> img))

let visualise_cimage (img: cimage) = visualise_argb_image (cimage_to_argb img)

let juliaGreyscale (p0: point) (n: i32): cimage =
  juliaEscapes p0 n >-> \(_, i) -> let i' = f32.i32 i / f32.i32 n
                                   in (i', i', i', 1f32)

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

-- the ultraPalette from Accelerate.
let mk_palette (points: i32) (ix: i32): argb.colour =
  let p = f32.i32 ix / f32.i32 points

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

let log2 (x: f32) = f32.log x / f32.log 2.0

let escape_to_colour (limit: i32) (points: i32)
                     (z: complex, n: i32): argb.colour =
  if limit == n then argb.black
  else let smooth = log2 (log2 (complex.mag z))
       let scale = 256.0
       let shift = 1664.0
       let ix = i32.f32 (f32.sqrt (f32.i32 n + 1.0 - smooth) * scale + shift)
       in mk_palette points (ix %% points)

let mandelbrot_greyscale t =
  mandelbrotGreyscale 100 |> rotate' t |> visualise_cimage

let julia_greyscale (xuser: f32) (yuser: f32) t =
  juliaGreyscale (xuser, yuser) 100 |> rotate' t |> visualise_cimage

let mandelbrot_colour t =
  mandelbrotImage 100 (escape_to_colour 100 2048) |>
  rotate' t |> visualise_argb_image

let julia_colour (user_x: f32) (user_y: f32) t =
  juliaEscapes (user_x, user_y) 100 |> rotate' t
  |> (>-> escape_to_colour 100 2048) |> visualise_argb_image

let figure_7_15 t =
  altRings |> shiftXor t |> uscale 0.05f32
  |> (>-> boolToColour) |> visualise_argb_image

let fancy t =
  cond (altRings |> shiftXor t |> uscale 0.05f32)
       (mandelbrotImage 100 (escape_to_colour 100 2048) |> rotate' t)
       (mandelbrotGreyscale 100 |> rotate' t |> cimage_to_argb)
  |> visualise_argb_image

import "lib/github.com/diku-dk/lys/lys"

type text_content = i32
module lys: lys with text_content = text_content = {
  let grab_mouse = false

  type mode = #mandelbrot_greyscale
            | #julia_greyscale
            | #mandelbrot_colour
            | #julia_colour
            | #figure_7_15
            | #fancy

  type state = {screen_size: {height:i64, width:i64},
                image_size: {height:f32, width:f32},
                centre: (f32, f32),
                userpos: (f32, f32),
                mouse: (i32, i32),
                mode: mode,
                time: f32,
                zooming: f32,
                paused: bool}

  let init _ h w: state = {
      screen_size = {height=h, width=w},
      image_size = {height=10, width=10},
      centre = (0, 0),
      userpos = (0.5, 0.5),
      mouse = (0, 0),
      mode = #mandelbrot_greyscale,
      time = 0,
      zooming = 0,
      paused = false
    }

  let resize h w (s: state) =
    let h_shrink = f32.i64 s.screen_size.height / f32.i64 h
    let w_shrink = f32.i64 s.screen_size.width / f32.i64 w
    in s with screen_size = {height=h, width=w}
         with image_size = {width=s.image_size.width * (1/w_shrink),
                            height=s.image_size.height * (1/h_shrink)}

  let keydown k (s: state) =
    if      k == '1' then s with mode = #mandelbrot_greyscale
    else if k == '2' then s with mode = #julia_greyscale
    else if k == '3' then s with mode = #mandelbrot_colour
    else if k == '4' then s with mode = #julia_colour
    else if k == '5' then s with mode = #figure_7_15
    else if k == '6' then s with mode = #fancy
    else if k == SDLK_KP_PLUS then s with zooming = -0.01
    else if k == SDLK_KP_MINUS then s with zooming = 0.01
    else if k == ' ' then s with paused = !s.paused
    else s

  let keyup k (s: state) =
    if      k == SDLK_KP_PLUS then s with zooming = 0
    else if k == SDLK_KP_MINUS then s with zooming = 0
    else s

  let diff (x1: i32, y1: i32) (x2, y2) = (x2 - x1, y2 - y1)

  let move_pixels (s: state) (dx, dy): state =
    let x_per_pixel = s.image_size.width / f32.i64 s.screen_size.width
    let y_per_pixel = s.image_size.height / f32.i64 s.screen_size.height
    in s with centre = (s.centre.0 + x_per_pixel * f32.i32 dx,
                        s.centre.1 + y_per_pixel * f32.i32 dy)

  let mouse (buttons: i32) x y (s: state) =
    let dpos = diff (x,y) s.mouse
    let s = s with mouse = (x,y)
    let s = if (buttons & 1) == 1
            then move_pixels s dpos
            else s
    let s = if (buttons & 4) == 4
            then s with userpos = (f32.i32 x / f32.i64 s.screen_size.width,
                                   f32.i32 y / f32.i64 s.screen_size.height)
            else s
    in s

  let do_zoom factor (s: state) =
    s with image_size = {width = s.image_size.width * (1+factor),
                         height = s.image_size.height * (1+factor)}

  let step td (s: state) =
    let td' = if s.paused then 0 else td
    in do_zoom s.zooming s with time = s.time + td'

  let event (e: event) (s: state) =
    match e
    case #step td -> step td s
    case #keydown {key} -> keydown key s
    case #keyup {key} -> keyup key s
    case #mouse {buttons, x, y} -> mouse buttons x y s
    case #wheel {dx=_, dy} -> do_zoom (-(f32.i32 dy)/100) s
    case _ -> s

  let render (s: state): [][]argb.colour =
    let render f = f
                   s.time
                   s.screen_size.width s.screen_size.height
                   s.image_size.width s.image_size.height
                   s.centre.0 s.centre.1
    in match s.mode
       case #mandelbrot_greyscale -> render mandelbrot_greyscale
       case #julia_greyscale -> render (julia_greyscale s.userpos.0 s.userpos.1)
       case #mandelbrot_colour -> render mandelbrot_colour
       case #julia_colour -> render (julia_colour s.userpos.0 s.userpos.1)
       case #figure_7_15 -> render figure_7_15
       case #fancy -> render fancy

  type text_content = text_content

  let text_format () = "FPS: %d"

  let text_content (fps: f32) (_: state): text_content =
    i32.f32 fps

  let text_colour = const argb.yellow
}
