-- Mandelbrot visualisation with 32-bit floats.

import "/futlib/math"
import "generic_mandelbrot"

module mandelbrot = mandelbrot f32

entry render_mandelbrot (screenX: i32) (screenY: i32)
                        (xcentre: f32) (ycentre: f32) (width: f32)
                        (limit: i32) (radius: f32)
                        : [screenX][screenY][3]u8 =
  mandelbrot.render_mandelbrot screenX screenY xcentre ycentre width limit radius
