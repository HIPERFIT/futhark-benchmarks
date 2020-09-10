import "lib/github.com/diku-dk/sorts/radix_sort"

import "brute-force"
import "kd-traverse"

let sortQueriesByLeavesRadix [n] (num_bits: i32) (leaves: [n]i32) : ([n]i32, [n]i32) =
  unzip <| radix_sort_by_key (\(l,_) -> l) num_bits i32.get_bit (zip leaves (map i32.i64 (iota n)))

let iterationSorted [q][n][d][k][N][num_leaves][ppl]
            (STEP: i32) (h: i32)
            (kd_tree: [q](i32,f32,i32))
            (leaves:  [num_leaves][ppl][d]f32)
            -- ^ invariant
            (ongoing_queries: [n][d]f32)
            (knns:        [n][k](i32, f32))
            (last_leaves: [n]i32)
            (stacks:      [n]i32)
            (dists:       [n]f32)
            (ord_knns:   *[N][k](i32,f32))
            (query_inds:  [n]i32)
            (i: i32) (visited: *[]i32)
          : ([][d]f32, [][k](i32, f32), []i32, []i32, []f32, *[N][k](i32,f32), []i32, i32, *[]i32) =

  -- start at old leaf and find a new leaf, until done!
  let (knns', new_leaves, new_stacks, new_dists) = unzip4 <|
    map5 (\ query knn leaf_ind stack dist ->
            let knn' = bruteForce query knn (leaf_ind*i32.i64 ppl,
                                             leaves[leaf_ind])
            let wnn  = knn'[k-1].1
            let (leaf_ind', stack', dist') =
              traverseOnce h kd_tree (query, wnn)
                           (leaf_ind, stack, dist)
            in  (knn', leaf_ind', stack', dist')
         ) ongoing_queries knns last_leaves stacks dists

  let (new_leaves_all, sort_inds) = sortQueriesByLeavesRadix (h+2) new_leaves
  -- we need (h+2) bits because the finish leaf is represented by num_leaves

  let num_valid = map (\l -> if l < i32.i64 num_leaves then 1 else 0) new_leaves_all
                  |> reduce_comm (+) 0i32
                  |> i64.i32

  let ongoing_queries'= (gather2D ongoing_queries sort_inds)[:num_valid]
  let query_inds_all  = gather1D query_inds sort_inds
  let knns_all        = gather2D knns' sort_inds
  let (valid, terminated) = split num_valid (zip query_inds_all knns_all)
  let (query_inds', knns'')= unzip valid
  let (terminated_inds, terminated_knns) = unzip terminated
  let ord_knns' = scatter2D ord_knns terminated_inds terminated_knns
  let (new_stacks', new_dists') = unzip <|
      (gather1D (zip new_stacks new_dists) sort_inds)[:num_valid]

  let visited = if (i != 0) && ((i%STEP) == 0)
                then let visited[i/STEP - 1] = i32.i64 num_valid in visited
                else visited

  let new_leaves' = new_leaves_all[:num_valid]
                        
  in  (ongoing_queries', knns'', new_leaves', new_stacks', new_dists', ord_knns', query_inds', i+1, visited)

