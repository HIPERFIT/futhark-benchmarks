-- Rodinia's Back-Propagation benchmark translated to Futhark.
--
-- ==
-- compiled input @ data/small.in
-- output @ data/small.out
-- input @ data/medium.in

let eta(): f32       = 0.3
let momentum(): f32  = 0.3
let init_zero(): bool = false

let squash(x: f32): f32 = 1.0 / (1.0 + f32.exp(-x))
let fabs   (x: f32): f32 = if x < 0.0 then 0.0 - x else x

-- Computational kernels

let bpnn_output_error [n] (target: [n]f32, output: [n]f32): (f32, [n]f32) =
    let (errs, delta) = unzip (
        map ( \(t: f32, o: f32): (f32,f32)  ->
                let d = o * (1.0 - o) * (t - o)
                in ( if d < 0.0 then 0.0-d else d, d ))
            (zip target output ))
    let err = f32.sum errs
    in  ( err, delta )


let bpnn_hidden_error [no][nh] (delta_o: [no]f32, who: [nh][no]f32, hidden: [nh]f32): (f32, [nh]f32) =
    let (errs, delta_h) = unzip (
        map ( \(hidden_el: f32, who_row: []f32): (f32,f32)  ->
                let prods  = map2 (*) delta_o who_row
                let sumrow = f32.sum prods
                let new_el = hidden_el * (1.0-hidden_el) * sumrow
                in ( fabs(new_el), new_el ))
            (zip hidden who))
    let err = f32.sum errs
    in  ( err, delta_h )

let bpnn_adjust_weights [ndelta][nlym1][nly] (delta: [ndelta]f32, ly: [nlym1]f32, w: [nly][ndelta]f32, oldw: [nly][ndelta]f32): ([nly][ndelta]f32, [nly][ndelta]f32) =
  let lyext = map( \k  ->
                        if k < 1 then 1.0 else #[unsafe] ly[k-1])
                 (iota nly)
  in unzip (map (\(w_row: []f32, oldw_row: []f32, lyk: f32)  ->
                    unzip (map ( \(w_el: f32, oldw_el: f32, delta_el: f32): (f32,f32)  ->
                                   let new_dw = eta()*delta_el*lyk + momentum()*oldw_el
                                   in ( w_el+new_dw, new_dw ))
                           (zip3 (w_row) (oldw_row) delta)))
            (zip3 w oldw lyext))


let bpnn_layerforward_GOOD [n1][n2] (l1: [n1]f32, conn: [n1][n2]f32, conn_fstrow: [n2]f32): [n2]f32 =
  let connT     = transpose(conn)
  let res_tmp   = map ( \(conn_tr_row: [n1]f32): f32  ->
                            let prods = map2 (*) conn_tr_row l1
                            in f32.sum prods)
                      connT
  in map (\(pr: f32, conn0: f32): f32  -> squash(pr+conn0))
  (zip (res_tmp) (conn_fstrow))


let bpnn_layerforward [n1][n2] (l1: [n1]f32, conn: [n1][n2]f32, conn_fstrow: [n2]f32): [n2]f32 =
  let connT     = transpose(conn)
  let res_map   = map ( \(conn_tr_row: [n1]f32): [n1]f32  ->
                        map2 (*) conn_tr_row l1)
                      connT

  -- FIXME: nasty hack to avoid fusion, which presently causes the
  -- kernel extractor to sequentialise the reduction.

  let res_tmp   = map f32.sum (opaque res_map)

  in map ( \(pr: f32, conn0: f32): f32  -> squash(pr+conn0))
  (zip (res_tmp) (conn_fstrow))

--------------------------------------------------------/

let bpnn_train_kernel [n_in][n_out][n_inp1][n_hid][n_hidp1]
                     ( input_units:  [n_in]f32
                     , target: [n_out]f32
                     , input_weights: [n_inp1][n_hid]f32
                     , hidden_weights: [n_hidp1][n_out]f32
                     , input_prev_weights: [n_inp1][n_hid]f32
                     , hidden_prev_weights: [n_hidp1][n_out]f32
                     ): ( f32, f32
                        , [n_inp1][n_hid]f32
                        , [n_hidp1][n_out]f32
                        ) =
    let (inpweightsP_row0,inpweightsP) = split (1)  input_weights
    let hidden_units = bpnn_layerforward(input_units,
                                         inpweightsP :> [n_in][n_hid]f32,
                                         inpweightsP_row0[0])

    let (hidweightsP_row0,hidweightsP) = split (1) hidden_weights
    let output_units = bpnn_layerforward(hidden_units,
                                         hidweightsP :> [n_hid][n_out]f32,
                                         hidweightsP_row0[0])

    let (out_err, output_delta) = bpnn_output_error(target, output_units)
    let (hid_err, hidden_delta) = bpnn_hidden_error(output_delta,
                                                    hidweightsP :> [n_hid][n_out]f32,
                                                    hidden_units)

    let (hidden_weights, hidden_prev_weights) =
            bpnn_adjust_weights(output_delta, hidden_units, hidden_weights, hidden_prev_weights)
    let ( input_weights,  input_prev_weights) =
            bpnn_adjust_weights(hidden_delta,  input_units,  input_weights,  input_prev_weights)

    in (out_err, hid_err, input_weights, hidden_weights)


----------------------------------------------------/

let sobolIndR [num_bits] (dirVct: [num_bits]i32) (n: i32): f32 =
    -- placed norm_fact here to check that hoisting does its job!
    let norm_fact = 1.0 / ( f32.i64(1 << num_bits) + 1.0 )
    let n_gray = (n >> 1) ^ n
    let res = 0
    let res = loop (res) for (i, v) in (zip (iota num_bits) dirVct) do
        let t = i32.i64 (1 << i)
        in if (n_gray & t) == t
           then res ^ v
           else res
    in f32.i32 res * norm_fact

let bpnn_randomize_weights (m: i64) (n: i64) (offset: i32) (dirVct: []i32): ([m][n]f32,[m]f32) =
    let mat =
      tabulate_2d m n (\i j ->
                    -- (n+1) needed to create the sob num
                    -- as in the original Rodinia code.
                    let offs = i32.i64 (i*(n+1)) + offset + 1
                    in sobolIndR dirVct (offs + i32.i64 j))
    let vct =
      tabulate m (\i -> sobolIndR dirVct (offset+i32.i64 (i*(n+1))))
    in (mat, vct)

let bpnn_randomize_row (m: i64) (offset: i32) (dirVct: []i32): [m]f32 =
    map (sobolIndR(dirVct)) (map (+offset) (map i32.i64 (iota m)))

let bpnn_constant_row (m: i64) (value: f32): [m]f32 =
    replicate m value

let bpnn_zero_weights (m: i64) (n: i64): [m][n]f32 =
    replicate m (replicate n 0.0)

----------------------------------------------------/

let bpnn_create (n_in: i64) (n_inp1: i64) (n_hid: i64)
                (n_hidp1: i64) (n_out: i64) (offset: i32) (dirVct: []i32)
                : ( [n_in]f32
                  , [n_out]f32
                  ,([n_inp1][n_hid]f32, [n_inp1]f32)
                  ,([n_hidp1][n_out]f32, [n_hidp1]f32)
                  , [n_inp1][n_hid]f32
                  , [n_hidp1][n_out]f32
                  ) =
  -- [#n_out]
  let target = bpnn_constant_row n_out 0.1

  -- [#n_in][#n_hidden]f32
  let (offset, (input_weights, input_weights_fstcol)) =
    if init_zero()
    then (  offset, (bpnn_zero_weights n_inp1 n_hid, replicate n_inp1 0.0) )
    else (  offset+i32.i64 n_inp1*i32.i64 n_hidp1,
            bpnn_randomize_weights n_inp1 n_hid offset dirVct  )

  -- [#n_hidden][#n_out]f32
  let (hidden_weights, hidden_weights_fstcol) =
    bpnn_randomize_weights n_hidp1 n_out offset dirVct
  let offset = offset + i32.i64 (n_hidp1*(n_out+1))

  --[#n_in][#n_hidden]f32
  let input_prev_weights = bpnn_zero_weights n_inp1 n_hid

  --[#n_hidden][#n_out]f32
  let hidden_prev_weights= bpnn_zero_weights n_hidp1 n_out

  --[#n_in]
  let input_units = bpnn_randomize_row n_in (offset+1) dirVct

  in ( input_units, target,
       (input_weights, input_weights_fstcol),
       (hidden_weights, hidden_weights_fstcol),
       input_prev_weights, hidden_prev_weights)

let consColumn [m][n] (mat: [m][n]f32, col: [m]f32): [m][]f32 =
  let np1 = n+1
  in map (\(matrow, colelm)  ->
            map (\k -> if k < 1 then colelm else #[unsafe] matrow[k-1])
                (iota np1))
          (zip mat col)

let main [num_bits] (n_in: i32) (dirVct: [num_bits]i32): ( f32, f32, [][]f32, [][]f32 ) =
    let (n_inp1, n_hid, n_hidp1, n_out) = (i64.i32 n_in+1, 16, 16+1, 1)
    let (   input_units, target,
           ( input_weights,  input_weights_fstcol),
           (hidden_weights, hidden_weights_fstcol),
            input_prev_weights, hidden_prev_weights) =
        bpnn_create (i64.i32 n_in) n_inp1 n_hid n_hidp1 n_out 1 dirVct

    let ( out_err, hid_err, input_weights, hidden_weights ) =
        bpnn_train_kernel(  input_units, target,
                            input_weights,  hidden_weights,
                            input_prev_weights, hidden_prev_weights )

    let  input_weights' = consColumn( input_weights,  input_weights_fstcol)
    let hidden_weights' = consColumn(hidden_weights, hidden_weights_fstcol)

    in (out_err, hid_err, input_weights', hidden_weights')
