import "buildKDtree"
import "knn-iteration"

let k = 5i64

let propagate [m][q][d][n] (ref_pts: [m][d]f32) (indir: [m]i32)
                           (kd_tree: [q](i32,f32,i32))
                           (queries: [n][d]f32) : ([][k]i32, [][k]f32, []i32, i32, i32) = 
  let (median_dims, median_vals, _) = unzip3 kd_tree
  let num_nodes  = q
  let num_leaves = num_nodes + 1
  let h = (log2 (i32.i64 num_leaves)) - 1
  let ppl = m / num_leaves
  let leaves = unflatten num_leaves ppl ref_pts

  let query_leaves = map (findLeaf median_dims median_vals h) queries
  let (query_leaves, query_inds) = sortQueriesByLeavesRadix (h+1) query_leaves
  let queries = gather2D queries query_inds

  let knns     = replicate n (replicate k (-1i32, f32.inf))
  let ord_knns = replicate n (replicate k (0, f32.inf))
  let STEP = 64
  let visited = replicate (num_leaves/STEP) (-1i32) --replicate ((h+1)*32) (-1i32)

  let (_, _, _, _, _, ord_knns, _, loop_count, visited) =
      loop (ongoing_queries, knns, last_leaves,  stacks,           dists,              ord_knns, query_inds, i,    visited) =
           (queries,         knns, query_leaves, replicate n 0i32, replicate n 0.0f32, ord_knns, copy query_inds, 0i32, visited)
      while (length ongoing_queries > 0) do
        iterationSorted (i32.i64 STEP) h kd_tree leaves ongoing_queries knns last_leaves stacks dists ord_knns query_inds i visited

  let (knn_inds0, knn_dsts) = unzip <| map unzip <| ord_knns
  let knn_inds = map (\kinds -> map (\ind -> indir[ind]) kinds) knn_inds0
  in (knn_inds, knn_dsts, visited, i32.i64 num_leaves, loop_count)

-- ==
-- entry: main
--
-- no_gtx780 no_k40 compiled random input { 256i32 [2097152][7]f32 [10000000][7]f32 }
-- output @ valid-data/knn-ppl-256-m-2097152-n-10000000-d-7-k-5.out.gz

entry main [m][d][n] (defppl: i32) (refs: [m][d]f32) (queries: [n][d]f32) =
    let (height, num_inner_nodes, _ppl, m') =
      computeTreeShape (i32.i64 m) defppl
    let (refs_pts, indir, median_dims, median_vals, clanc_eqdim) =
      mkKDtree height (i64.i32 num_inner_nodes) (i64.i32 m') refs
    let (_knn_inds, knn_dsts, _visited, _num_leaves, _loop_count) =
            propagate refs_pts indir (zip3 median_dims median_vals clanc_eqdim) queries
    in  knn_dsts[:(n/64)]
