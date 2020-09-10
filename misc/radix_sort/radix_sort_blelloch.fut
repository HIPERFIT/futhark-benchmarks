-- A least significant digit radix sort to test out `scatter`; this variant
-- directly based on [1], which is apparently also the basis for one of
-- Accelerate's example programs.
--
-- [1] G. E. Blelloch. "Prefix sums and their applications." Technical Report
--     CMU-CS-90-190. Carnegie Mellon University. 1990.

let plus_scan [n] (x: [n]i32): [n]i32 =
  scan (+) 0 x

let plus_prescan [n] (x: [n]i32): [n]i32 =
  let xshifted = map (\i  -> if i == 0 then 0 else x[i - 1]) (iota n)
  in scan (+) 0 xshifted

let permute [n] (a: [n]u32, index: [n]i32): [n]u32 =
  scatter (copy a) (map i64.i32 index) a

let plus_scan_reverse_order [n] (x: [n]i32): [n]i32 =
  let xreversed = reverse x
  let x' = plus_scan xreversed
  let x'reversed = reverse x'
  in x'reversed

let split_blelloch [n] (a: [n]u32, flags: [n]i32): [n]u32 =
  let i_down = plus_prescan(map (1-) flags)
  let i_up = map (i32.i64 n-) (plus_scan_reverse_order(flags))
  let index = map3 (\flag up down -> if flag == 1 then up else down)
                   flags i_up i_down
  in permute(a, index)

let split_radix_sort [n] (a: [n]u32, number_of_bits: i32): [n]u32 =
  loop (a) for i < number_of_bits do
    let ai = map (\a -> i32.u32 (a >> u32.i32 i) & 1) a
    in split_blelloch(a, ai)
