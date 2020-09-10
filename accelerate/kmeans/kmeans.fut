-- This is an implementation of K-means specialised for only two
-- dimensions.  See [0] for an implementation that works for arbitrary
-- dimensionality - it's not all that more complicated.
--
-- [0]: https://github.com/diku-dk/futhark-benchmarks/tree/master/rodinia/kmeans
-- ==
-- compiled input @ data/trivial.in
-- output @ data/trivial.out
-- compiled input @ data/k5_n50000.in
-- output @ data/k5_n50000.out
-- compiled input @ data/k5_n200000.in
-- output @ data/k5_n200000.out

type point = (f32,f32)

let add_points ((x1,y1): point) ((x2,y2): point): point =
  (x1+x2, y1+y2)

let scale_point ((x,y): point) (s: f32) : point =
  (x*s, y*s)

let euclid_dist_2 ((x1,y1): point) ((x2,y2): point): f32 =
  (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)

let closest_point (p1: (i32, f32)) (p2: (i32, f32)): (i32, f32) =
  if p1.1 < p2.1 then p1 else p2

let find_nearest_point [k] (pts: [k]point) (pt: point): i32 =
  let (i, _) = foldl (\acc (i, p) -> closest_point acc (i32.i64 i, euclid_dist_2 pt p))
                     (0, f32.inf)
                     (zip (indices pts) pts)
  in i

let centroids_of [n] (k: i64) (points: [n]point) (membership: [n]i32): [k]point =
  let cluster_sizes =
    reduce_by_index (replicate k 0) (+) 0 (map i64.i32 membership) (replicate n 1)
  let cluster_sums =
    reduce_by_index (replicate k (0,0)) add_points (0,0) (map i64.i32 membership) points
  in map2 scale_point cluster_sums (map (1/) (map r32 cluster_sizes))

let continue [k] (old_centres: [k]point) (cur_centres: [k]point): bool =
  let changed ((x1,y1), (x2,y2)) =
    f32.abs (x1-x2) > 0.01 || f32.abs(y1-y2) > 0.01
  in (any changed (zip old_centres cur_centres))

let main [n] (k: i32) (points_in: [n][2]f32): ([][2]f32, i32) =
  let k = i64.i32 k

  -- Transform from 2D-array to array of pairs.
  let points = map (\point -> (point[0], point[1])) points_in

  -- Assign arbitrary initial cluster centres.
  let initial_cluster_centres = take k points

  let recentre (cluster_centres: [k]point) =
        -- For each point, find the cluster with the closest centroid.
        let new_membership = map (find_nearest_point cluster_centres) points
        -- Then find the new centres of the clusters.
        in centroids_of k points new_membership

  let (_, final_cluster_centres, i) =
    loop (old_centres, cur_centres, i) =
      (initial_cluster_centres, recentre initial_cluster_centres, 0)
    while continue old_centres cur_centres do
      (cur_centres, recentre cur_centres, i+1)
  in (map (\(x, y) -> [x,y]) final_cluster_centres, i)
