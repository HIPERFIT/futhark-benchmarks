-- Code and comments based on
-- https://github.com/kkushagra/rodinia/blob/master/openmp/lavaMD/main.c
--
-- ==
-- input @ data/3_boxes.in
-- output @ data/3_boxes.out
-- compiled input @ data/10_boxes.in.gz
-- output @ data/10_boxes.out.gz

let dot ((ax,ay,az), (bx,by,bz)): f32 =
  ax*bx + ay*by + az*bz

-----------------------------------------
-- Main Computational Kernel of lavaMD --
-----------------------------------------
let main [number_boxes][par_per_box][num_neighbors]
             (alpha:  f32)
             (box_coefs_0: [number_boxes]i32)
             (box_coefs_1: [number_boxes]i32)
             (box_coefs_2: [number_boxes]i32)
             (box_coefs_3: [number_boxes]i32)
             (box_nnghs_0 : [num_neighbors][number_boxes]i32)
             (box_nnghs_1 : [num_neighbors][number_boxes]i32)
             (box_nnghs_2 : [num_neighbors][number_boxes]i32)
             (box_nnghs_3 : [num_neighbors][number_boxes]i32)
             (box_num_nghbs: [number_boxes]i32)
             (rv_0: [number_boxes][par_per_box]f32)
             (rv_1: [number_boxes][par_per_box]f32)
             (rv_2: [number_boxes][par_per_box]f32)
             (rv_3: [number_boxes][par_per_box]f32)
             (qv: [number_boxes][par_per_box]f32)
      : ([number_boxes][par_per_box]f32,
         [number_boxes][par_per_box]f32,
         [number_boxes][par_per_box]f32,
         [number_boxes][par_per_box]f32) =
  let box_coefs = zip4 (box_coefs_0) (box_coefs_1) (box_coefs_2) (box_coefs_3)
  let box_nnghs = map4 zip4 box_nnghs_0 box_nnghs_1 box_nnghs_2 box_nnghs_3
  let rv = map4 zip4 rv_0 rv_1 rv_2 rv_3
  let a2 = 2.0*alpha*alpha in
  unzip4 <| map unzip4 (map2 (\box_num_nghbs' (l: i64): [par_per_box](f32,f32,f32,f32)  ->
        let rA = rv[l]
        in
        map  (\(rA_el:  (f32,f32,f32,f32) ): (f32,f32,f32,f32)  -> --(i32 i) ->
                let (rai_v, rai_x, rai_y, rai_z) = rA_el in -- rA[i]
                let acc = (0.0,0.0,0.0,0.0)
                in loop(acc) for k < box_num_nghbs'+1 do
                          let pointer = if (k > 0)
                                        then let (_,_,_,num) = #[unsafe] box_nnghs[k-1, l] in num
                                        else i32.i64 l

                          let (_,_,_,first_j) = #[unsafe] box_coefs[pointer]
                          let rB = #[unsafe] rv[first_j]
                          let qB = #[unsafe] qv[first_j]
                          ---------------------------------------------------------
                          -- Important note: rB and qB are invariant to the      --
                          -- second map on rA -> can be blocked in shared memory --
                          ---------------------------------------------------------
                          let pres =
                            map2 (\(rbj_v,rbj_x,rbj_y,rbj_z) qbj ->
                                   let r2   = rai_v + rbj_v - dot((rai_x,rai_y,rai_z), (rbj_x,rbj_y,rbj_z))
                                   let u2   = a2*r2
                                   let vij  = f32.exp(-u2)
                                   let fs   = 2.0 * vij
                                   let d_x  = rai_x  - rbj_x
                                   let d_y  = rai_y  - rbj_y
                                   let d_z  = rai_z  - rbj_z
                                   let fxij = fs * d_x
                                   let fyij = fs * d_y
                                   let fzij = fs * d_z
                                   in (qbj*vij, qbj*fxij, qbj*fyij, qbj*fzij))
                                 rB qB

                          let (r1, r2, r3, r4) =
                            reduce (\(a1,a2,a3,a4) (b1,b2,b3,b4) ->
                                    (a1+b1, a2+b2, a3+b3, a4+b4))
                                   (0,0,0,0) pres
                          let (a1, a2, a3, a4) = acc
                          in  (a1+r1, a2+r2, a3+r3, a4+r4)
            ) rA
     ) box_num_nghbs (iota(number_boxes)))
