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

let B0: i64 = 64

let fInd (B: i64) (y:i32) (x:i32): i32 = y*(i32.i64 B+1) + x
let max3 (x:i32, y:i32, z:i32) = if x < y
                                 then if y < z then z else y
                                 else if x < z then z else x
let mkVal [l2][l] (y:i32) (x:i32) (pen:i32) (inp_l:[l2][l2]i32) (ref_l:[l][l]i32) : i32 = #[unsafe]
  max3( ( (inp_l[y - 1, x - 1])) + ( ref_l[y-1, x-1])
      , ( (inp_l[y, x - 1])) - pen
      , ( (inp_l[y - 1, x])) - pen
      )

let intraBlockPar [len] (B: i64)
                               (penalty: i32)
                               (inputsets: [len][len]i32)
                               (reference2: [len][len]i32)
                               (b_y: i64) (b_x: i64)
                               : [B][B]i32 =
  let ref_l = reference2[b_y * B + 1: b_y * B + 1 + B,
                         b_x * B + 1: b_x * B + 1 + B] :> [B][B]i32

  -- inp_l is the working memory
  let inp_l = replicate ((B+1)*(B+1)) 0i32
              |> unflatten (B+1) (B+1)

  -- Initialize inp_l with the already processed the column to the left of this
  -- block
  let inp_l[0:B+1, 0] = inputsets[b_y * B : b_y * B + B + 1, b_x * B]

  -- Initialize inp_l with the already processed the row to above this block
  let inp_l[0, 1:B+1] = inputsets[b_y * B, b_x * B + 1 : b_x * B + B + 1]

  -- Process the first half (anti-diagonally) of the block
  let inp_l = loop inp_l for m < B do
       let inds =
            -- tabulate over the m'th anti-diagonal before the middle
            tabulate B (\tx ->  (
                    if tx > m then (-1, -1)
                    else let ind_x = i32.i64 (tx + 1)
                         let ind_y = i32.i64 (m - tx + 1)
                         in  (i64.i32 ind_y, i64.i32 ind_x)))
        let vals =
            -- tabulate over the m'th anti-diagonal before the middle
            tabulate B (\tx ->  (
                    if tx > m then 0
                    else let ind_x = i32.i64 (tx + 1)
                         let ind_y = i32.i64 (m - tx + 1)
                         let v = mkVal ind_y ind_x penalty inp_l ref_l
                         in  v))
        in  scatter_2d inp_l inds vals

  -- Process the second half (anti-diagonally) of the block
  let inp_l = loop inp_l for m < B-1 do
        let m = B - 2 - m
        let inds = tabulate B (\tx ->  (
                    if tx > m then (-1, -1)
                    else let ind_x = i32.i64 (tx + B - m)
                         let ind_y = i32.i64 (B - tx)
                         in  ((i64.i32 ind_y, i64.i32 ind_x)) )
                )
        let vals =
            -- tabulate over the m'th anti-diagonal after the middle
            tabulate B (\tx ->  (
                    if tx > m then (0)
                    else let ind_x = i32.i64 (tx + B - m)
                         let ind_y = i32.i64 (B - tx)
                         let v = mkVal ind_y ind_x penalty inp_l ref_l
                         in  v ))
        in  scatter_2d inp_l inds vals

  let inp_l2 = inp_l
  in  inp_l2[1:B+1,1:B+1] :> [B][B]i32


let updateBlocks [q][len] (B: i64)
                            (blk: i64)
                            (mk_b_y: (i32 -> i32))
                            (mk_b_x: (i32 -> i32))
                            (block_inp: [q][B][B]i32)
                            (inputsets:  *[len][len]i32) =
  let (inds, vals) = unzip (
    tabulate (blk*B*B) (\gid ->
                 let B2 = i32.i64 (B*B)
                 let gid = i32.i64 gid
                 let (bx, lid2) = (gid / B2, gid % B2)
                 let (ty, tx)   = (lid2 / i32.i64 B, lid2 % i32.i64 B)

                 let b_y = mk_b_y bx
                 let b_x = mk_b_x bx
                 let v = #[unsafe] block_inp[bx, ty, tx]
                 in  ((i64.i32 (i32.i64 B*b_y + 1 + ty),
                       i64.i32 (i32.i64 B*b_x + tx + 1)),
                      v)))
  in  scatter_2d inputsets inds vals


let main [len] (penalty : i32)
                 (inputsets : *[len][len]i32)
                 (reference : *[len][len]i32) : *[len][len]i32 =
  #[unsafe]
  let worksize = len - 1
  let B = i64.min worksize B0

  -- worksize should be a multiple of B0
  let B = assert (worksize % B == 0) B

  let block_width = i32.i64 <| worksize / B

  -- First anti-diagonal half of the entire input matrix
  let inputsets =
    loop inputsets for blk < block_width do
        let blk = i64.i32 (blk + 1)
        let block_inp =
          -- Process an anti-diagonal of independent blocks
          tabulate blk (\b_x ->
                let b_y = blk-1-b_x
                in  intraBlockPar B penalty inputsets reference b_y b_x
          )

        let mkBY bx = i32.i64 (blk - 1) - bx
        let mkBX bx = bx
        in  updateBlocks B blk mkBY mkBX block_inp inputsets

  -- Second anti-diagonal half of the entire input matrix
  let inputsets =
    loop inputsets for blk < block_width-1 do
        let blk = i64.i32 (block_width - 1 - blk)
        let block_inp =
          -- Process an anti-diagonal of independent blocks
          tabulate blk (\bx ->
                let b_y = i64.i32 block_width - 1 - bx
                let b_x = bx + i64.i32 block_width - blk
                in  intraBlockPar B penalty inputsets reference b_y b_x
          )

        let mkBY bx = block_width - 1 - bx
        let mkBX bx = bx + block_width - i32.i64 blk
        in  updateBlocks B blk mkBY mkBX block_inp inputsets

  in inputsets
