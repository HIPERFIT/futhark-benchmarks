-- Code and comments based on
-- https://github.com/kkushagra/rodinia/blob/master/openmp/nw
--
-- ==
--
-- compiled input @ data/large.in.gz
-- output @ data/large.out.gz

-- compiled input @ data/tiny.in
-- output @ data/tiny.out
-- compiled input @ data/small.in
-- output @ data/small.out
-- compiled input @ data/medium.in.gz
-- output @ data/medium.out.gz

let B : i64 = 32

let fInd (y:i32) (x:i32): i32 = y*(i32.i64 B+1) + x
let max3 (x:i32, y:i32, z:i32) = if x < y 
                                 then if y < z then z else y
                                 else if x < z then z else x
let mkVal [l2][l] (y:i32) (x:i32) (pen:i32) (inp_l:[l2]i32) (ref_l:[l][l]i32) : i32 = #[unsafe]
  max3( ( (inp_l[fInd (y-1) (x-1)])) + ( ref_l[y-1, x-1])
      , ( (inp_l[fInd y (x-1)])) - pen
      , ( (inp_l[fInd (y-1) x])) - pen
      )
    
let intraBlockPar [lensq][len] (penalty: i32)
                               (inputsets: [lensq]i32) 
                               (reference2: [len][len]i32) 
                               (b_y: i64) (b_x: i64)
                               : [B][B]i32 =
  -- index   =   base + cols * BLOCK_SIZE * b_index_y + BLOCK_SIZE * b_index_x + tx + ( cols + 1 );
  -- for ( int ty = 0 ; ty < BLOCK_SIZE ; ty++)
  --   REF(ty, tx) =  reference_d[index + cols * ty];

  let ref_l = replicate (B*B) 0
  let ref_l = loop ref_l for i < B do
                scatter ref_l (map (\tid->i*B+tid) (iota B))
                        (tabulate B (\j -> reference2[i+b_y*B+1,
                                                      b_x*B+1+j]))
  let ref_l = unflatten B B ref_l

  let inp_l = replicate ((B+1)*(B+1)) 0i32

  -- index_nw =  base + cols * BLOCK_SIZE * b_index_y + BLOCK_SIZE * b_index_x;
  let index_nw =  len * B * b_y + B * b_x                
  let inp_l[0] = #[unsafe] inputsets[index_nw]

  --index_w   = base + cols * BLOCK_SIZE * b_index_y + BLOCK_SIZE * b_index_x + ( cols );
  let index_w = len*B*b_y + B*b_x + len 
  -- SCORE((tx + 1), 0) = input_itemsets_d[index_w + cols * tx];
  let inp_l = scatter inp_l (map (\tx->(tx+1)*(B+1)) (iota B))
                      (map (\tx->#[unsafe] inputsets[index_w+len*tx]) (iota B))

  --index_n   = base + cols * BLOCK_SIZE * b_index_y + BLOCK_SIZE * b_index_x + tx + ( 1 );
  let index_n = len*B*b_y + B*b_x + 1
  -- SCORE(0, (tx + 1)) = input_itemsets_d[index_n];
  let inp_l = scatter inp_l (map (+1) (iota B))
                      (map (\tx->#[unsafe] inputsets[index_n+tx]) (iota B))

  let inp_l = loop inp_l for m < B do
        let (inds, vals) = unzip (
            tabulate B (\tx ->  (
                    if tx > m then (-1, 0)
                    else let ind_x = i32.i64 (tx + 1)
                         let ind_y = i32.i64 (m - tx + 1)
                         let v = mkVal ind_y ind_x penalty inp_l ref_l
                         in  (i64.i32 (fInd ind_y ind_x), v))))
        in  scatter inp_l inds vals

  let inp_l = loop inp_l for m < B-1 do
        let m = B - 2 - m
        let (inds, vals) = unzip (
            map (\tx ->  (
                    if tx > m then (-1, 0)
                    else let ind_x = i32.i64 (tx + B - m)
                         let ind_y = i32.i64 (B - tx)
                         let v = mkVal ind_y ind_x penalty inp_l ref_l
                         in  (i64.i32 (fInd ind_y ind_x), v) )
                ) (iota B) )
        in  scatter inp_l inds vals

  let inp_l2 = unflatten (B+1) (B+1) inp_l
  in  inp_l2[1:B+1,1:B+1] :> [B][B]i32


let updateBlocks [q][lensq] (len: i32) (blk: i64)
                            (mk_b_y: (i32 -> i32)) 
                            (mk_b_x: (i32 -> i32)) 
                            (block_inp: [q][B][B]i32) 
                            (inputsets:  *[lensq]i32) =
  let (inds, vals) = unzip (
    tabulate (blk*B*B) (\gid ->
                 let B2 = i32.i64 (B*B)
                 let gid = i32.i64 gid
                 let (bx, lid2) = (gid / B2, gid % B2)
                 let (ty, tx)   = (lid2 / i32.i64 B, lid2 % i32.i64 B)

                 let b_y = mk_b_y bx
                 let b_x = mk_b_x bx
                 let v = #[unsafe] block_inp[bx, ty, tx]
                 let ind = (i32.i64 B*b_y + 1 + ty) * len + (i32.i64 B*b_x + tx + 1)
                 in  (i64.i32 ind, v)))
  in  scatter inputsets inds vals


-- `len-1` should be a multiple of 16 
let main [lensq] (penalty : i32) 
                 (inputsets : *[lensq]i32) 
                 (reference : *[lensq]i32) : [lensq]i32 =
  let len = i32.f32 (f32.sqrt (f32.i64 lensq))
  let worksize = len - 1
  let block_width = worksize / i32.i64 B
  let reference2 = unflatten (i64.i32 len) (i64.i32 len) reference

  -- First loop BEGINS
  let inputsets =
    loop inputsets for blk < block_width do
        let blk = i64.i32 (blk + 1)
        let block_inp =
          tabulate blk (\b_x ->
                let b_y = blk-1-b_x
                in  intraBlockPar penalty inputsets reference2 b_y b_x
          )

        let mkBY bx = i32.i64 (blk - 1) - bx
        let mkBX bx = bx
        in  updateBlocks len blk mkBY mkBX block_inp inputsets
  -- First loop ENDS

  -- Second loop BEGINS
  let inputsets =
    loop inputsets for blk < block_width-1 do
        let blk = i64.i32 (block_width - 1 - blk)
        let block_inp =
          tabulate blk (\bx ->
                let b_y = i64.i32 block_width - 1 - bx
                let b_x = bx + i64.i32 block_width - blk
                in  intraBlockPar penalty inputsets reference2 b_y b_x
          )

        let mkBY bx = block_width - 1 - bx
        let mkBX bx = bx + block_width - i32.i64 blk
        in  updateBlocks len blk mkBY mkBX block_inp inputsets
    -- Second loop ENDS

  in  inputsets
