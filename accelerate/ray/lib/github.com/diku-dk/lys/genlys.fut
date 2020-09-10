-- | ignore

-- This file exists as a wrapper that defines entry points in the
-- specific form that liblys.c requires.  It is copied into place and
-- modified by the rules in common.mk.

module m = import "lys"

type^ state = m.lys.state

entry init (seed: u32) (h: i32) (w: i32): state =
  m.lys.init seed (i64.i32 h) (i64.i32 w)

entry grab_mouse: bool =
  m.lys.grab_mouse

entry resize (h: i32) (w: i32) (s: state): state =
  m.lys.resize (i64.i32 h) (i64.i32 w) s

entry key (e: i32) (key: i32) (s: state): state =
  let e' = if e == 0 then #keydown {key} else #keyup {key}
  in m.lys.event e' s

entry mouse (buttons: i32) (x: i32) (y: i32) (s: state): state =
  m.lys.event (#mouse {buttons, x, y}) s

entry wheel (dx: i32) (dy: i32) (s: state): state =
  m.lys.event (#wheel {dx, dy}) s

entry step (td: f32) (s: state): state =
  m.lys.event (#step td) s

entry render (s: state) = m.lys.render s

entry text_colour (s: state): u32 =
  m.lys.text_colour s

entry text_format: []u8 = m.lys.text_format ()

entry text_content (render_duration: f32) (s: state) =
  m.lys.text_content render_duration s
