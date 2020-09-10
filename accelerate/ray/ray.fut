import "lib/github.com/diku-dk/lys/lys"
import "lib/github.com/athas/matte/colour"
module trace = import "trace"
import "types"


type text_content = (i32, i32)
module lys: lys with text_content = text_content = {
  type eye = {pos: position, a: f32, b: f32}

  type state = { height: i64
               , width: i64
               , fov: i32
               , eye: eye
               , forward_sgn: i32
               , sideway_sgn: i32
               , limit: i32
               , time: f32
   }

  let init _ h w: state =
    { height = h, width = w,
      fov = 100, limit = 4,
      eye = {pos={x=50,y= -100,z= -700}, a=f32.pi/2, b=0},
      forward_sgn = 0,
      sideway_sgn = 0,
      time = 0}

  let resize h w (s: state): state =
    s with height = h with width = w

  let keydown (k: i32) (s: state) =
    if      k == 'w' then s with forward_sgn = 1
    else if k == 's' then s with forward_sgn = -1
    else if k == 'a' then s with sideway_sgn = -1
    else if k == 'd' then s with sideway_sgn = 1
    else if k == 'z' then s with limit = i32.max (s.limit - 1) 0
    else if k == 'x' then s with limit = s.limit + 1
    else s

  let keyup (k: i32) (s: state) =
    if      k == 'w' then s with forward_sgn = 0
    else if k == 's' then s with forward_sgn = 0
    else if k == 'a' then s with sideway_sgn = 0
    else if k == 'd' then s with sideway_sgn = 0
    else s

  let grab_mouse = true

  let move_speed: f32 = 1000
  let forwards td ({pos=_, a,b}: eye) (s: i32) =
    let amount = move_speed * f32.i32 s * td
    in {x = amount * f32.cos(a) * f32.cos(b),
        y = amount * f32.sin(b),
        z = amount * f32.sin(a) * f32.cos(b)}

  let sideways td ({pos=_, a,b=_}: eye) (s: i32) =
    let amount = move_speed * f32.i32 s * td
    in {x = amount * f32.cos(a + f32.pi/2),
        y = 0 : f32,
        z = amount * f32.sin(a + f32.pi/2)}

  let event (e: event) (s: state) =
    match e
    case #keydown {key} -> keydown key s
    case #keyup {key} -> keyup key s
    case #step td ->
      s with time = s.time + td
        with eye.pos = s.eye.pos
                       |> vec3.((+forwards td s.eye s.forward_sgn))
                       |> vec3.((+sideways td s.eye s.sideway_sgn))
    case #mouse {buttons, x, y} ->
      if buttons == 0 then
        s with eye.a = s.eye.a + f32.i32 x/f32.i64 s.width
          with eye.b = f32.min (f32.max (s.eye.b + f32.i32 y/f32.i64 s.height)
                                        (-f32.pi/2+0.001))
                               (f32.pi/2-0.001)
      else s
    case _ -> s


  let render (s: state): [][]argb.colour =
    trace.main (i32.i64 s.width) (i32.i64 s.height) s.fov
               s.eye.pos.x s.eye.pos.y s.eye.pos.z
               s.eye.a s.eye.b
               s.limit s.time

  type text_content = text_content

  let text_format () = "FPS: %d; rendering limit: %d"

  let text_content (fps: f32) (s: state): text_content =
    (i32.f32 fps, s.limit)

  let text_colour = const argb.yellow

}
