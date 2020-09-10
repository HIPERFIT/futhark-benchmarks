let log2 x = (loop (y,c) = (x,0i32) while y > 1i32 do (y >> 1, c+1)).1

let partition2Ind [n] (cs: [n]bool) : ([n]i32, i32) =
    let tfs = map (\f -> if f then 1 else 0) cs
    let isT = scan (+) 0 tfs
    let ffs = map (\f -> if f then 0 else 1) cs
    let isF0 = scan (+) 0 ffs

    let i = isT[n-1]
    let isF = map (+ i) isF0 
    let inds = map3 (\ c iT iF ->
                        if c then iT-1 else iF-1
                    ) cs isT isF
    in (inds, i)

let sumSqrs [d] (xs: [d]f32) (ys: [d]f32) : f32 =
    map2 (\x y -> let z = x-y in z*z) xs ys |> reduce (+) 0.0f32

let sumSqrsSeq [d] (xs: [d]f32) (ys: [d]f32) : f32 =
    loop (res) = (0.0f32) for (x,y) in (zip xs ys) do
        let z = x-y in res + z*z

let gather1D 't [m] (arr1D: [m]t) (inds: [m]i32) : *[m]t =
    map (\ind -> arr1D[ind] ) inds

let gather2D 't [m][d] (arr2D: [m][d]t) (inds: [m]i32) : *[m][d]t =
    map (\ind -> map (\j -> arr2D[ind,j]) (iota d) ) inds

let scatter2D [m][k][n] 't (arr2D: *[m][k]t) (qinds: [n]i32) (vals2D: [n][k]t) : *[m][k]t =
  let nk = n*k
  let flat_qinds = map (\i -> let (d,r) = (i32.i64 i / i32.i64 k,
                                           i32.i64 i % i32.i64 k)
                              in i64.i32 (qinds[d]*i32.i64 k + r)
                       ) (iota nk)
  let res1D = scatter (flatten arr2D) flat_qinds ((flatten vals2D) :> [nk]t) 
  in  unflatten m k res1D 


let getParent (node_index: i32) = (node_index-1) / 2

let isLeaf (h: i32) (node_index: i32) =
    node_index >= ((1 << (h+1)) - 1)

-- the k'th ancestor of `node_ind` can be computed with
-- the formula: `(node_ind + 1 - (2^k)) / (2^k)`, for example
-- the parent           (k==1): `(node_ind - 1) / 2`
-- the grandparent      (k==2): `(node_ind - 3) / 4`
-- the grandgrandparent (k==3): `(node_ind - 7) / 8`
let compute_Kth_ancestor (k: i32) (node_ind: i32) =
    let tpk = 1 << k
    in  (node_ind + 1 - tpk) / tpk

let findNodeLevel (node: i32) : i32 =
 ( loop (lev, idx) = (0i32, node)
     while idx > 0i32 do
       (lev+1i32, getParent idx) ).0

-- given a tree `node1` at level `lev` and another tree leaf `leaf`,
-- this function computes the closest common ancestor of `node1` and `node2`
-- `h` is the height of the binary tree (without leaves).
let findClosestCommonAncestor (h: i32) (lev: i32) (node1: i32) (leaf: i32) : i32 =
    let node2 = compute_Kth_ancestor (h+1-lev) leaf
    let (res,_) =
      loop (node1, node2) while node1 != node2 do
        (getParent node1, getParent node2)
    in  res
