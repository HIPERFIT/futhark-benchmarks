import "lib/github.com/athas/matte/colour"
import "types"

module type object = {
  type object
  val distance_to: object -> position -> direction -> (bool, f32)
}

type sphere = {position: position,
               colour: argb.colour,
               shine: f32,
               radius: f32}

module sphere: object with object = sphere = {
  type object = sphere

  let distance_to (sphere: sphere)
                  (origin: position)
                  (direction: direction): (bool, f32) =
    let pos = sphere.position
    let radius = sphere.radius
    let p = vec3.(origin + (scale (dot (pos - origin) direction) direction))
    let d_cp = vec3.(norm (p - pos))
    let sep = vec3.(p - origin)
    let miss = d_cp >= radius || vec3.dot sep direction <= 0.0
    in if miss
       then (false, f32.inf)
       else (true, vec3.norm sep - f32.sqrt (radius*radius - d_cp * d_cp))
}

let sphere_normal (sphere: sphere) (point: position): direction =
  vec3.normalise (point vec3.- sphere.position)

type plane = {position: position,
              normal: direction,
              colour: argb.colour,
              shine: f32}

type plane_check = plane

module plane: object with object = plane = {
  type object = plane

  let distance_to (plane: plane)
                  (origin: position)
                  (direction: direction): (bool, f32) =
    let pos = plane.position
    let normal = plane.normal
    let theta = vec3.dot direction normal
    in if theta >= 0.0
       then (false, f32.inf)
       else (true, vec3.(dot (pos - origin) normal) / theta)
}

let checkers (c: argb.colour) ({x,y=_,z}: position): argb.colour =
  let v1 = i32.f32 (x/100.0) % 2
  let v2 = i32.f32 (z/100.0) % 2
  let v3 = i32.bool (x < 0.0)
  let v4 = i32.bool (z < 0.0)
  in if v1 ^ v2 ^ v3 ^ v4 == 1
     then c
     else c |> argb.dim |> argb.dim |> argb.dim

type objects [nsphere] [nplane] = {spheres: [nsphere]sphere, planes: [nplane]plane}
