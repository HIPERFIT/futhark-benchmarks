-- multi-dimensional spatial Euclid distance square
--
-- Much inspiration has been gained from the kmeans-vector package on
-- Hackage.
--
-- ==
-- tags { futhark-opencl futhark-c }
-- nobench input @ data/trivial.in
-- output @ data/trivial.out
-- compiled input @ data/100.in
-- output @ data/100.out
-- input @ data/204800.in.gz
-- output @ data/204800.out
-- input @ data/kdd_cup.in.gz
-- output @ data/kdd_cup.out

let euclid_dist_2 [d] (pt1: [d]f32) (pt2: [d]f32): f32 =
  f32.sum (map (\x->x*x) (map2 (-) pt1 pt2))

let closest_point (p1: (i32,f32)) (p2: (i32,f32)): (i32,f32) =
  if p1.1 < p2.1 then p1 else p2


let find_nearest_point [k][d] (pts: [k][d]f32) (pt: [d]f32): i32 =
  let (i, _) = foldl (\acc (i, p) -> closest_point acc (i32.i64 i, euclid_dist_2 pt p))
                     (0, f32.inf)
                     (zip (indices pts) pts)
  in i

let add_centroids [d] (x: [d]f32) (y: [d]f32): *[d]f32 =
  map2 (+) x y

let centroids_of [n][d] (k: i64) (points: [n][d]f32) (membership: [n]i32): [k][d]f32 =
  let points_in_clusters =
    reduce_by_index (replicate k 0) (+) 0 (map i64.i32 membership) (replicate n 1)

  let cluster_sums =
    reduce_by_index (replicate k (replicate d 0)) (map2 (+)) (replicate d 0)
                    (map i64.i32 membership)
                    points

  in map2 (\point n -> map (/r32 (if n == 0 then 1 else n)) point)
          cluster_sums points_in_clusters

let main [n][d]
        (threshold: i32) (k: i32) (max_iterations: i32)
        (points: [n][d]f32): ([][]f32, i32) =
  let k = i64.i32 k

  -- Assign arbitrary initial cluster centres.
  let cluster_centres = take k points
  -- Also assign points arbitrarily to clusters.
  let membership = map i32.i64 (map (%k) (iota n))
  let delta = threshold + 1
  let i = 0
  let (_,cluster_centres,_,i) =
    loop (membership, cluster_centres, delta, i)
    while delta > threshold && i < max_iterations do
      -- For each point, find the cluster with the closest centroid.
      let new_membership = map (find_nearest_point cluster_centres) points
      -- Then, find the new centres of the clusters.
      let new_centres = centroids_of k points new_membership
      let delta = i32.sum (map (\b -> if b then 0 else 1)
                               (map2 (==) membership new_membership))
      in (new_membership, new_centres, delta, i+1)
  in (cluster_centres, i)
