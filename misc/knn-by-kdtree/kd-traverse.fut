import "util"

-- give the median dimensions and values for each internal k-d tree node,
--   finds the leaf to which the query naturally belongs to.
let findLeaf [q][d] (median_dims: [q]i32) (median_vals: [q]f32)
                    (height: i32) (query: [d]f32) =
  let leaf =
    loop (node_index) = (0)
      while !(isLeaf height node_index) do
        if query[median_dims[node_index]] <= median_vals[node_index]
        then (node_index+1)*2-1
        else (node_index+1)*2
  in leaf - i32.i64 q

-- This is implemented for 1-dim
-- height: the height of the tree (without leaves)
-- query: our query (1-dim)
-- wnn:    the worse nearest neighbor (distance)
-- median_dims: the dimension that was split in each node of the k-d tree
-- median_vals: the median value in each node of the k-d tree
-- last_leaf: last visited leaf
-- stack:  array of booleans (which should be represented as an int)
--         which denotes whether the `second child` of a node has
--         been visited already
-- Results: the index of the new leaf and the new stack
--
let traverseOnce [q] [d] (height: i32)
                 (kd_tree: [q](i32,f32,i32))
                 (query: [d]f32, wnn: f32)
                 (last_leaf: i32, stack: i32, dist: f32) : (i32, i32, f32) =

  let (median_dims, median_vals, clanc_eqdim) = unzip3 kd_tree
  let last_leaf = last_leaf + i32.i64 q
  let no_leaf   = 2*q + 1

  -- helper functions for reading and writing the stack,
  --   which is maintained as an int
  let getPackedInd (stk: i32) (ind: i32) : bool = 
    let b = stk & (1<<ind) in b != 0
  let setPackedInd (stk: i32) (ind: i32) (v: bool) =
    let fst = stk & ((1<<ind)-1)
    let snd = (stk >> (ind+1)) << (ind+1)
    let mid = if v then (1 << ind) else 0
    in  ( (fst | snd) | mid )

  let getLevel (node_idx: i32) : i32 = log2 (node_idx+1) 
  let getAncSameDimContrib (q_m_i: f32) (node_stack: i32) (node: i32) : f32 =
    (loop (idx, res) = (node, 0.0f32) 
      while (idx >= 0) do
        let anc = clanc_eqdim[idx] in
        if anc == (-1i32) then (-1i32, 0.0f32)
        else 
          let anc_lev = getLevel anc
          let is_anc_visited = getPackedInd node_stack anc_lev
          in  if !is_anc_visited then (anc, res)
              else (-1i32, median_vals[anc] - q_m_i)
    ).1

  -- go back on the stack and find the new node to be visited
  let (parent_rec, stack, count, dist, rec_node) =
      loop (node_index, stack, count, dist, rec_node) =
           (last_leaf, stack, height, dist, -1)
            while (node_index != 0) && (rec_node < 0) do
                let parent = getParent node_index
                let scnd_visited = getPackedInd stack count --stack[count]

                -- visiting condition
                let q_m_d   = query[median_dims[parent]]
                let cur_med_dst = median_vals[parent] - q_m_d
                let cur_med_sqr = cur_med_dst * cur_med_dst

                let prv_med_dst = getAncSameDimContrib q_m_d stack parent
                let prv_med_sqr = prv_med_dst * prv_med_dst

                let dist_minu = f32.abs(dist - cur_med_sqr + prv_med_sqr)
                let dist_plus = f32.abs(dist - prv_med_sqr + cur_med_sqr)

                in if scnd_visited
                   then -- continue backing-up towards the root
                        (parent, stack, count-1, dist_minu, -1)

                   else -- the node_index is actually the `first` child of parent,
                        -- since the `second` has not been visited yet
                       let to_visit = (f32.sqrt dist_plus) < wnn
                       --let to_visit = (f32.abs cur_med_dst) < wnn
                       in  if !to_visit
                           then (parent, stack, count-1, dist, -1)
                           else -- update the stack
                                let fst_node = node_index
                                let snd_node = if (fst_node % 2) == 0 then fst_node-1 else fst_node+1
                                let stack = setPackedInd stack count true 
                                -- let stack[count] = true
                                in  (parent, stack, count, dist_plus, snd_node)
  -- find a new leaf
  let (new_leaf, new_stack, _) =
      if parent_rec == 0 && rec_node == -1
      then -- we are done, we are at the root node
           -- and its second child has been visited
           (i32.i64 no_leaf, stack, 0)

      else -- now traverse downwards by computing `first`
           -- and stop when you discovered a new leaf
           loop (node_index, stack, count) = 
                (rec_node, stack, count)
           while !(isLeaf height node_index) do
              let count = count+1
              let stack = setPackedInd stack count false
              -- let stack[count] = false
              let node_index =
                  if query[median_dims[node_index]] <= median_vals[node_index]
                  then (node_index+1)*2-1
                  else (node_index+1)*2
              in (node_index, stack, count)

  in (new_leaf-i32.i64 q, new_stack, dist)
