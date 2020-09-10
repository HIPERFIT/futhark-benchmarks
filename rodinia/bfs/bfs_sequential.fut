-- A naive, sequential version of BFS.  Its purpose is to exist as a simple,
-- working solution.
-- ==
-- tags { nobench }
--
-- compiled input @ data/4096nodes.in
-- output @ data/4096nodes.out
-- compiled input @ data/512nodes_high_edge_variance.in.gz
-- output @ data/512nodes_high_edge_variance.out
-- input @ data/graph1MW_6.in.gz
-- output @ data/graph1MW_6.out.gz
-- input @ data/64kn_32e-var-1-256-skew.in.gz
-- output @ data/64kn_32e-var-1-256-skew.out

let node_work [n][e] (tid: i32,
                      cost: *[n]i32,
                      nodes_start_index: [n]i32,
                      nodes_n_edges: [n]i32,
                      edges_dest: [e]i32,
                      graph_visited: [n]bool,
                      graph_mask: *[n]bool,
                      updating_graph_mask: *[n]bool): (*[n]i32, *[n]bool, *[n]bool) =
  let start_index = nodes_start_index[tid]
  let n_edges = nodes_n_edges[tid]
  let graph_mask[tid] = false
  let (cost, updating_graph_mask) =
    loop ((cost, updating_graph_mask)) for i in start_index..<start_index+n_edges do
      let id = edges_dest[i]
      let visited = graph_visited[id]
      in if ! visited
         then
           let cost[id] = cost[tid] + 1
           let updating_graph_mask[id] = true
           in (cost, updating_graph_mask)
         else
           (cost, updating_graph_mask)
    in (cost, graph_mask, updating_graph_mask)

let step [n][e] (cost: *[n]i32)
                (nodes_start_index: [n]i32)
                (nodes_n_edges: [n]i32)
                (edges_dest: [e]i32)
                (graph_visited: [n]bool)
                (graph_mask: *[n]bool)
                (updating_graph_mask: *[n]bool): (*[n]i32, *[n]bool, *[n]bool) =
  let (active_indices, _) = unzip (filter (.1) (zip (map i32.i64 (iota n)) graph_mask))

  -- This loop is a kernel in Rodinia.  Futhark's regularity makes this a bit
  -- tricky to express as a map.
  in loop ((cost, graph_mask, updating_graph_mask))
     for indices_i < length active_indices do
       let i = active_indices[indices_i]
       in node_work(i, cost, nodes_start_index, nodes_n_edges,
                    edges_dest, graph_visited,
                    graph_mask, updating_graph_mask)

import "common"

let main = common_main step
