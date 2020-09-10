-- Generate well-conditioned matrix to avoid floating point errors.
--
-- Algorithm from Rodinia.

let main (n: i64): [n][n]f32 =
  let lambda = -0.001
  let f (p: i64) =
    (if p >= n
     then let i = -n + p + 1
          in 10.0 * f32.exp(lambda*f32.i64(i))
     else let i = n - p - 1
          in 10.0 * f32.exp(lambda*f32.i64(i)))
  in map (\i -> map (\j -> f (n-1-i+j)) (iota n))
         (iota n)
