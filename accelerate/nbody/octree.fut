import "lib/github.com/athas/vector/vspace"
import "lib/github.com/diku-dk/sorts/radix_sort"
import "radixtree"

module vec3 = mk_vspace_3d f32
type vec3 = vec3.vector

let vec (x, y, z) : vec3 = {x,y,z}

type mass = f32
type position = vec3.vector
type velocity = vec3.vector
type body = {position: position,
             mass: mass,
             velocity: velocity}

type octnode = {body: body,
               children: [8]i32,
               parent: i32,
               tree_level: i32,
               is_leaf: bool}

-- | takes radix tree, assumes mortons are sorted.
let mk_octree [n] [m] (L: [n]body) (inners : [m]inner) : []octnode =
  -- (6)-(7) from karras2012. In addition, we also assign child pointers.
  let root_delta = inners[0].delta_node
  let edge delta_c delta_p : i32 = (delta_c / 3) - (delta_p / 3)
  let get_delta ptr = match ptr
                      case #leaf _i -> u64.num_bits + 2 -- u32.num_bits + 1
                      case #inner i -> inners[i].delta_node
  let edges =
    map (\n ->
           let left_edge = edge (get_delta n.left) n.delta_node
           let right_edge = edge (get_delta n.right) n.delta_node
           in {left = left_edge, right = right_edge}
        ) inners
  let prefix_sum = edges |> map (\{left, right} ->
                                  (i32.bool (left > 0)) +
                                  (i32.bool (right > 0)))
                         |> scan (+) 0
  let size = i64.i32 (last prefix_sum) + 1
  let ranges = map (+ 1) (rotate (-1) prefix_sum) with [0] = 1
  let (is, vs) = map (\x -> (i64.i32 x, 1)) (tail ranges) |> unzip
  let I_rad = reduce_by_index (replicate size 0) (+) 0 is vs
              |> scan (+) (0 : i32)

  let octree =
    map2 (\i rp ->
           let i = i32.i64 i
           -- Always relative to root.
           let radix_parent = inners[rp]
           let edges' = edges[rp]
           let beginning = ranges[rp]

           let node_level = (radix_parent.delta_node) / 3 - root_delta / 3
           let left_edge = edges'.left
           let is_left = left_edge > 0 && beginning == i
           let tree_level =
             if i == 0 then 0 else
               node_level + (if is_left then left_edge else edges'.right)
           let child = if is_left then radix_parent.left
                       else radix_parent.right
           let is_leaf = i != 0 &&
                         match child
                           case #leaf _i -> true
                           case #inner _i -> false

           let unpack ridx : i32 =
             match ridx
               case #inner ind -> ind
               case #leaf ind -> ind

           let body = if is_leaf then
                      --- if we are a leaf, child will be a pointer to
                      --- a node in the leaf array, so this is safe.
                      let b = L[unpack child]
                        in b with position = vec3.scale b.mass b.position
                      else {position = vec (0, 0, 0),
                            mass = 0.0,
                            velocity = vec (0, 0, 0)}

           let parent =
             if i == 0 then -1 -- root node
             else let ridx = rp
                  let parent' = -1
                  let (parent', _) =
                    loop (parent', ridx) while parent' == -1 do -- O(3)
                      let parent_ridx = inners[ridx].parent
                      in if parent_ridx == -1 -- octree root is parent
                        then (0, parent_ridx)
                         else
                          -- we checked parent_ridx != -1, so it must
                          -- be a pointer to a node in the radix tree,
                          -- which are assumed to be valid.
                          let parent_lchild =
                            unpack inners[parent_ridx].left
                          let parent_edges = edges[parent_ridx]

                          let is_left = parent_lchild == ridx
                          let edge_weight = if is_left then
                                              parent_edges.left
                                            else parent_edges.right
                          in if edge_weight != 0 -- find parent in
                                                 -- range
                             then
                              let start = ranges[parent_ridx]
                              in (if is_left
                                 then start
                                 else start +
                                  (1 * (i32.bool (parent_edges.left > 0)))
                                 , parent_ridx)
                             else (parent', parent_ridx)
                  in parent'

           let ne : i32 = -1
           let children =
             let arr = replicate 8 ne in
             if is_leaf then arr
             else
               let ridx = if i == 0 then 0 else unpack child
               let lvl : i32 = 0
               let to_check  = [2, 0, 0]
               let child_idx : i32 = 0
               let (arr', _, _, _, _) =
                 loop (arr, to_check, child_idx, ridx, lvl)
                 while to_check[0] > 0 do
                     --- lvl can be 0, 1 or 2.  Hence it is safe to
                     --- index into to_check with lvl.
                     --- child_idx should also stay between 0 and 7 when
                     --- loop condition is true.
                     if to_check[lvl] == 2 then -- check left edge
                       if edges[ridx].left != 0 then
                         let start = ranges[ridx]
                         let arr[child_idx] = start
                         let child_idx = child_idx + 1
                         let to_check[lvl] = 1
                         in (arr, to_check, child_idx, ridx, lvl)
                       else -- move down one level
                         let ridx = unpack inners[ridx].left
                         let lvl = lvl + 1
                         let to_check[lvl] = 2
                         in (arr, to_check, child_idx, ridx, lvl)
                     else if to_check[lvl] == 1 then -- check right edge
                       if edges[ridx].right != 0 then
                         let start = ranges[ridx]
                         let arr[child_idx] =
                           start + (i32.bool (edges[ridx].left > 0)) * 1
                         let child_idx = child_idx + 1
                         let to_check[lvl] = 0
                         in (arr, to_check, child_idx, ridx, lvl)
                       else -- move down one level
                         let ridx = unpack inners[ridx].right
                         let lvl = lvl + 1
                         let to_check[lvl] = 2
                         in (arr, to_check, child_idx, ridx, lvl)
                     else
                       let ridx = inners[ridx].parent
                       let lvl = lvl - 1
                       let to_check[lvl] = to_check[lvl] - (1 : i32)
                       in (arr, to_check, child_idx, ridx, lvl)
               in arr'
           in {body, children, parent, tree_level, is_leaf}
         ) (iota size) I_rad
   in octree

-- | Expands a 10-bit integer into 30 bits by inserting 2 zeros after
-- each bit.
let expand_bits (v: u32) : u32 =
  let v = (v * 0x00010001) & 0xFF0000FF
  let v = (v * 0x00000101) & 0x0F00F00F
  let v = (v * 0x00000011) & 0xC30C30C3
  let v = (v * 0x00000005) & 0x49249249
  in v

let morton_3D {x,y,z} : u32 =
  let x = f32.min (f32.max(x * 1024) 0) 1023
  let y = f32.min (f32.max(y * 1024) 0) 1023
  let z = f32.min (f32.max(z * 1024) 0) 1023
  let xx = expand_bits (u32.f32 x)
  let yy = expand_bits (u32.f32 y)
  let zz = expand_bits (u32.f32 z)
  in xx * 4 + yy * 2 + zz

let mk_accelerator [n] (bodies: [n]body) : ([]octnode, f32, i32, [n]body) =
  let x_max = f32.maximum (map (.position.x) bodies)
  let y_max = f32.maximum (map (.position.y) bodies)
  let z_max = f32.maximum (map (.position.z) bodies)
  let x_min = f32.minimum (map (.position.x) bodies)
  let y_min = f32.minimum (map (.position.y) bodies)
  let z_min = f32.minimum (map (.position.z) bodies)
  let normalize {x,y,z} = {x=(x-x_min)/(x_max-x_min),
                           y=(y-y_min)/(y_max-y_min),
                           z=(z-z_min)/(z_max-z_min)}
  let morton = normalize >-> morton_3D
  -- (1)-(5) from karras2012.
  let bodies =
    radix_sort_by_key (\p -> morton p.position) u32.num_bits u32.get_bit bodies
  let mortons = map (\p -> morton p.position) bodies
  let inners = mk_radix_tree mortons
  -- (6)-(7).
  let octree = mk_octree bodies inners
  let root_delta = inners[0].delta_node
  let accelerator =
    let root_leaf_delta = u64.num_bits + 2 - root_delta
    let lvl =  (root_leaf_delta + 3 - 1) / 3
    in loop octree for i in (lvl .. lvl-1 ... 0) do
         map (\(n: octnode) ->
                if n.tree_level != i then n
                else
                let j = 0
                let b = n.body
                let (body', _) =
                  loop (b, j) while j < 8 && n.children[j] != -1 do
                    let cb = octree[n.children[j]].body
                    let pos' = b.position vec3.+ cb.position
                    let mass' = b.mass + cb.mass
                    in
                    ({position = pos',
                     mass = mass',
                     velocity = b.velocity},
                     j+1)
                in n with body = body') octree
  let side_max = f32.maximum [x_max - x_min, y_max - y_min, z_max - z_min]
  in (accelerator, side_max, root_delta, bodies)
