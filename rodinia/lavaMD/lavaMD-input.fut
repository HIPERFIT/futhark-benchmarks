------------------------------------------
-- Util: Sobol random number generation --
------------------------------------------
let sobolDirVcts(): [30]i32 = 
  [ 536870912, 268435456, 134217728, 67108864, 33554432, 16777216, 8388608, 4194304, 2097152, 1048576, 
    524288,    262144,    131072,    65536,    32768,    16384,    8192,    4096,    2048,    1024, 
    512,       256,       128,       64,       32,       16,       8,       4,       2,       1      ] 

let sobolInd(dirVct:  [30]i32, n: i64 ): i32 =
  let n_gray = (n >> 1) ^ n
  let res = 0
  in loop (res) for (i, v) in zip (iota 30) dirVct do
       let t = 1 << i
       in if (n_gray & t) == t
          then res ^ v
          else res

-----------------------
-----------------------
let main (boxes1d: i64) (par_per_box: i64) (num_neighbors: i64):
                        (f32,
                         []i32,
                         []i32,
                         []i32,
                         []i32,

                         [][]i32,
                         [][]i32,
                         [][]i32,
                         [][]i32,

                         []i32,

                         [][]f32,
                         [][]f32,
                         [][]f32,
                         [][]f32,

                         [][]f32) =
  let number_boxes = boxes1d * boxes1d * boxes1d
  let alpha        = 0.5
  let dirVct       = sobolDirVcts()

  ----------------------------------------
  -- 1. Initialize boxs' data structure --
  ----------------------------------------
  let boxes = 
    tabulate number_boxes
        (\nh  ->
          let k = nh % boxes1d
          let nr= nh / boxes1d
          let j = nr % boxes1d
          let i = nr / boxes1d

          -- current home box
          let box_coef = ( i32.i64 k, i32.i64 j, i32.i64 i, i32.i64 nh )
          -- initialize neighbor boxes
          let (box_nngh, cur_nn) = (replicate num_neighbors (0,0,0,0), 0)
          let (box_nngh_cur_nn) = loop (box_nngh, cur_nn) for nn < num_neighbors do
            let n = (nn % 3) - 1
            let nr= nn / 3
            let m = (nr % 3) - 1
            let l = (nr / 3) - 1
            let (cur_elem, next_cur_nn) =
              -- check if (this neighbor exists) and (it is not the same as home box)
              if((( 0<=(i+l) && 0<=(j+m) && 0<=(k+n)) &&
                  ((i+l)<boxes1d && (j+m)<boxes1d && (k+n)<boxes1d)) &&
                 (!(l==0 && m==0 && n==0)))
              then  let (x, y, z) = (k+n, j+m, i+l)
                    let number = (z * boxes1d * boxes1d) + (y * boxes1d) + x
                    in ( (i32.i64 x, i32.i64 y, i32.i64 z, i32.i64 number), cur_nn+1)
              else  ( (0, 0, 0, 0     ), cur_nn  )
            in let box_nngh[cur_nn] = cur_elem
               in (box_nngh, next_cur_nn)

          let (box_nngh, cur_nn) = box_nngh_cur_nn
          in  ( box_coef, box_nngh, cur_nn )

       )

  let (box_coefs, box_nnghs0, box_num_nghbs) = unzip3 boxes
  let box_nnghs = copy(transpose(box_nnghs0))

  ----------------------------------------------
  -- 2. Initialize input distances and charge --
  ----------------------------------------------
  let rqv = map  (\i  ->
                    map (\j  ->
                           let n = (i*par_per_box + j)*5 + 1
                           let s1= sobolInd(dirVct, n  )
                           let s2= sobolInd(dirVct, n+1)
                           let s3= sobolInd(dirVct, n+2)
                           let s4= sobolInd(dirVct, n+3)
                           let s5= sobolInd(dirVct, n+4)
                           in (f32.i32(s5%10 + 1) / 10.0, 
                               ( f32.i32(s1%10 + 1) / 10.0, f32.i32(s2%10 + 1) / 10.0
                               , f32.i32(s3%10 + 1) / 10.0, f32.i32(s4%10 + 1) / 10.0 )
                              )
                       ) (iota(par_per_box))
                ) (iota(number_boxes) )

  let (qv, rv) = map unzip rqv |> unzip

  in ( alpha,

       (unzip4(box_coefs)).0,
       (unzip4(box_coefs)).1,
       (unzip4(box_coefs)).2,
       (unzip4(box_coefs)).3,

       (map (map (.0)) box_nnghs),
       (map (map (.1)) box_nnghs),
       (map (map (.2)) box_nnghs),
       (map (map (.3)) box_nnghs),

       box_num_nghbs,

       (map (map (.0)) rv),
       (map (map (.1)) rv),
       (map (map (.2)) rv),
       (map (map (.3)) rv),

       qv)
