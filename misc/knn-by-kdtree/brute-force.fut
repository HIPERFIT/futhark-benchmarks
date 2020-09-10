
let sumSqrsSeq [d] (xs: [d]f32) (ys: [d]f32) : f32 =
    loop (res) = (0.0f32) for (x,y) in (zip xs ys) do
        let z = x-y in res + z*z

let bruteForce [m][d][k] (query: [d]f32) 
                         (knns0: [k](i32,f32))
                         (beg: i32, refs : [m][d]f32)
                       : [k](i32,f32) =
    --if query[0] == f32.lowest then copy knns else
    loop (knns) = (copy knns0)
      --for (i,refpt) in (zip (iota m) refs) do
      --  let dist = f32.sqrt <| sumSqrsSeq query refpt in
      for i < i32.i64 m do
        let dist = f32.sqrt <| sumSqrsSeq query (refs[i]) in
        if dist > knns[k-1].1 then knns -- early exit
        else let ref_ind = i+beg in
             let (_, _, knns') =
               loop (dist, ref_ind, knns) for j < k do
                 let cur_nn = knns[j].1  in
                 if dist >= cur_nn
                 then (dist, ref_ind, knns)
                 else let tmp_ind = knns[j].0
                      let knns[j] = (ref_ind, dist)
                      let ref_ind = tmp_ind
                      in  (cur_nn, ref_ind, knns)
             in  knns'

