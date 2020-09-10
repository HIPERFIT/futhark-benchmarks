-- Does not exactly reproduce the same results as the original PBBS
-- implementation.  Reasons unclear, but it's not due to numerical
-- issues.  The same geometry calculations work fine in other ray
-- tracers.
-- ==
-- compiled input @ data/angel.in.gz
-- compiled input @ data/dragon.in.gz
-- compiled input @ data/happy.in.gz

import "prim"
import "bvh"

type point = vec3.vector
type triangle = {a: point, b: point, c: point, index: i32}
type ray = {orig: point, dir: point}

type^ bvh [n] = bvh [n] triangle

let aabb_hit ({min, max}: aabb) ({orig, dir}: ray) (tmin: f32) (tmax: f32): bool =
  -- Unrolled loop.
  let iter min' max' origin' dir' tmin tmax =
    let invD = 1 / dir'
    let t0 = (min' - origin') * invD
    let t1 = (max' - origin') * invD
    let (t0, t1) = if invD < 0 then (t1, t0) else (t0, t1)
    let tmin = f32.max t0 tmin
    let tmax = f32.min t1 tmax
    in (tmin, tmax)
  let (tmin, tmax) =
    iter min.x max.x orig.x dir.x tmin tmax
  in if tmax <= tmin then false
     else let (tmin, tmax) =
            iter min.y max.y orig.y dir.y tmin tmax
          in if tmax <= tmin then false
             else let (tmin, tmax) =
                    iter min.z max.z orig.z dir.z tmin tmax
                  in !(tmax <= tmin)

let EPSILON: f32 = 0.0001

let triangle_hit ({orig=o,dir=d}: ray) ({a,b,c,index=_}: triangle) : f32 =
  let e1 = b vec3.- a
  let e2 = c vec3.- a
  let pvec = vec3.cross d e2
  let det = vec3.dot e1 pvec
  -- If determinant is zero then ray is parallel with the triangle
  -- plane.
  in if det > -EPSILON && det < EPSILON
     then 0
     else
     -- Calculate distance from 'a' to origin
     let invDet = 1/det
     -- Calculate distance from 'a' to origin
     let tvec = o vec3.- a
     -- u and v are the barycentric coordinates in triangle if u >= 0,
     -- v >= 0 and u + v <= 1
     let u = vec3.dot tvec pvec * invDet
     -- check against one edge and opposite point
     in if u < 0.0 || u > 1.0 then 0
        else
        let qvec = vec3.cross tvec e1
        let v = vec3.dot d qvec * invDet
        -- Check against other edges
        in if v < 0.0 || u + v > 1.0 then 0
           else -- Distance along the ray, i.e. intersect at o + t * d
             vec3.dot e2 qvec * invDet

let triangle_aabb ({a,b,c,index=_}: triangle) : aabb =
  {min = {x=a.x `f32.min` b.x `f32.min` c.x,
          y=a.y `f32.min` b.y `f32.min` c.y,
          z=a.z `f32.min` b.z `f32.min` c.z},
   max = {x=a.x `f32.max` b.x `f32.max` c.x,
          y=a.y `f32.max` b.y `f32.max` c.y,
          z=a.z `f32.max` b.z `f32.max` c.z}}

let bvh_hit [n] (bvh: bvh [n]) (r: ray) : i32 =
  let contains aabb = aabb_hit aabb r 0 f32.inf
  let closest_hit (j, t_max) i t =
    let h = triangle_hit r t
    in if h < t_max
       then (i, h)
       else (j, t_max)
  let (j, _) = bvh_fold contains closest_hit (-1, f32.inf) bvh
  in if j == -1 then -1 else bvh.L[j].index

let main [npoints][ntriangles][nrays]
         (points_enc: [npoints][3]f32)
         (triangles_enc: [ntriangles][3]i32)
         (rays_enc: [nrays][2][3]f32)
         : [nrays]i32 =
  -- The input data is in an encoded format, so the first thing we
  -- need to do is give it a nicer representation.

  let to_point p = {x=p[0], y=p[1], z=p[2]}
  let points = map to_point points_enc

  -- Despite what the documentation says, the point indices are
  -- actually 1-based.
  let to_triangle t i =
    {a=points[t[0]-1], b=points[t[1]-1], c=points[t[2]-1], index=i32.i64 i}
  let triangles = map2 to_triangle triangles_enc (iota ntriangles)

  let to_ray r = {orig = to_point r[0], dir = to_point r[1]}
  let rays = map to_ray rays_enc

  let bvh = bvh_mk triangle_aabb triangles

  in map (bvh_hit bvh) rays
