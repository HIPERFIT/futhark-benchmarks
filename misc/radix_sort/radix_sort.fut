-- A least significant digit radix sort to test out `scatter`.

let radix_sort_step [n] (xs: [n]u32, digit_n: i32): [n]u32 =
  let bits = map (\x -> i32.u32 (x >> u32.i32 digit_n) & 1) xs
  let bits_inv = map (\(b: i32): i32  -> 1 - b) bits
  let ps0 = scan (+) 0 (bits_inv)
  let ps0_clean = map2 (*) bits_inv ps0
  let ps1 = scan (+) 0 bits
  let ps0_offset = i32.sum bits_inv
  let ps1_clean = map (+ps0_offset) ps1
  let ps1_clean' = map2 (*) bits ps1_clean
  let ps = map2 (+) ps0_clean ps1_clean'
  let ps_actual = map (\x -> x-1) ps
  in scatter (copy xs) (map i64.i32 ps_actual) xs

let radix_sort [n] (xs: [n]u32): [n]u32 =
  loop (xs) for i < 32 do radix_sort_step(xs, i)
