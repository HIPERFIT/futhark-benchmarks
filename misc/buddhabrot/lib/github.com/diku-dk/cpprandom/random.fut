-- | Random number generation inspired by `<random>` in C++.
--
-- The overall idea is that you pass a low-level `rng_engine`@mtype,
-- that knows how to generate random integers, to a parametric module
-- that maps said integers to the desired distribution.  Since Futhark
-- is a pure language, the final random number function(s) will return
-- both the random number and the new state of the engine.  It is the
-- programmer's responsibility to ensure that the same state is not
-- used more than once (unless that is what is desired).  See the
-- [Examples](#examples) below.
--
-- ## Examples
--
-- This program constructs a uniform distribution of single precision
-- floats using the `minstd_rand`@term as the underlying RNG engine.
-- The `dist` module is constructed at the program top level, while we
-- use it at the expression level.  We use the `minstd_rand` module
-- for initialising the random number state using a seed, and then we
-- pass that state to the `rand` function in the generated
-- distribution module, along with a description of the distribution
-- we desire.  We get back not just the random number, but also the
-- new state of the engine.
--
-- ```
-- module dist = uniform_real_distribution f32 minstd_rand
--
-- let rng = minstd_rand.rng_from_seed [123]
-- let (rng, x) = dist.rand (1,6) rng
-- ```
--
-- The `rand` function of `uniform_real_distribution`@term simply
-- takes a pair of numbers describing the range.  In contrast,
-- `normal_distribution`@term takes a record specifying the mean and
-- standard deviation:
--
-- ```
-- module norm_dist = normal_distribution f32 minstd_rand
--
-- let (rng, y) = norm_dist.rand {mean=50, stddev=25} rng
-- ```
--
-- Since both `dist` and `norm_dist` have been initialised with the
-- same underlying `rng_engine`@mtype (`minstd_rand`@term), we can
-- re-use the same RNG state.  This is often convenient when a program
-- needs to generate random numbers from several different
-- distributions, as we still only have to manage a single RNG state.
--
-- ### Parallel random numbers
--
-- Random number generation is inherently sequential.  The `rand`
-- functions take an RNG state as input and produce a new RNG state.
-- This creates challenges when we wish to `map` a function `f` across
-- some array `xs`, and each application of the function must produce
-- some random numbers.  We generally don't want to pass the exact
-- same state to every application, as that means each element will
-- see the exact same stream of random numbers.  Common procedure is
-- to use `split_rng`, which creates any number of RNG states from
-- one, and then pass one to each application of `f`:
--
-- ```
-- let rngs = minstd_rand.split_rng n rng
-- let (rngs, ys) = unzip (map2 f rngs xs)
-- let rng = minstd.rand.join_rngs rngs
-- ```
--
-- We assume here that the function `f` returns not just the result,
-- but also the new RNG state.  Generally, all functions that accept
-- random number states should behave like this.  We subsequently use
-- `join_rngs` to combine all resulting states back into a single
-- state.  Thus, parallel programming with random numbers involves
-- frequently splitting and re-joining RNG states.  For most RNG
-- engines, these operations are generally very cheap.
--
-- ## See also
--
-- The `Sobol`@term@"sobol" module provides a very different
-- (and inherently parallel) way of generating random numbers, which
-- may be more suited for Monte Carlo applications.

-- Quick and dirty hashing to mix in something that looks like entropy.
-- From http://stackoverflow.com/a/12996028
local
let hash(x: i32): i32 =
  let x = u32.i32 x
  let x = ((x >> 16) ^ x) * 0x45d9f3b
  let x = ((x >> 16) ^ x) * 0x45d9f3b
  let x = ((x >> 16) ^ x)
  in i32.u32 x

-- | Low-level modules that act as sources of random numbers in some
-- uniform distribution.
module type rng_engine = {
  -- | A module for the type of integers generated by the engine.
  module int: integral
  -- | The state of the engine.
  type rng

  -- | Initialise an RNG state from a seed.  Even if the seed array is
  -- empty, the resulting RNG should still behave reasonably.  It is
  -- permissible for this function to process the seed array
  -- sequentially, so don't make it too large.
  val rng_from_seed: []i32 -> rng

  -- | Split an RNG state into several states.  Implementations of
  -- this function tend to be cryptographically unsound, so be
  -- careful.
  val split_rng: (n: i32) -> rng -> [n]rng

  -- | Combine several RNG states into a single state - typically done
  -- with the result of `split_rng`@term.
  val join_rng: []rng -> rng

  -- | Generate a single random element, and a new RNG state.
  val rand: rng -> (rng,int.t)

  -- | The minimum value potentially returned by the generator.
  val min: int.t

  -- | The maximum value potentially returned by the generator.
  val max: int.t
}

module type rng_distribution = {
  -- | The random number engine underlying this distribution.
  module engine: rng_engine

  -- | A module describing the type of values produced by this random
  -- distribution.
  module num: numeric

  -- | The dynamic configuration of the distribution.
  type distribution

  val rand: distribution -> engine.rng -> (engine.rng, num.t)
}

-- | A linear congruential random number generator produces numbers by
-- the recurrence relation
--
-- > X(n+1) = (a × X(n) + c) mod m
--
-- where *X* is the sequence of pseudorandom values, and
--
-- * *m, 0 < m* — "modulus"
--
-- * *a, 0 < a < m* — "multiplier"
--
-- * *c, 0 ≤ c < m* — "increment"
--
-- * *X(0), 0 ≤ X(0) < m* — "seed" or "initial value"
module linear_congruential_engine (T: integral) (P: {
  val a: T.t
  val c: T.t
  val m: T.t
}): rng_engine with int.t = T.t with rng = T.t = {
  type t = T.t
  type rng = t

  module int = T

  let rand (x: rng): (rng, t) =
    let rng' = (P.a T.* x T.+ P.c) T.%% P.m
    in (rng',rng')

  let rng_from_seed [n] (seed: [n]i32) =
    let seed' =
      loop seed' = 1 for i < n do
        u32.(((seed' >> 16) ^ seed') ^
             (i32 seed[i] ^ 0b1010101010101))
    in (rand (T.u32 seed')).1

  let split_rng (n: i32) (x: rng): [n]rng =
    map (\i -> x T.^ T.i32 (hash i)) (iota n)

  let join_rng [n] (xs: [n]rng): rng =
    reduce (T.^) (T.i32 0) xs

  let min = T.i32 0
  let max = P.m
}

-- | A random number engine that uses the *subtract with carry*
-- algorithm.  Presently quite slow.  The size of the state is
-- proportional to the long lag.
module subtract_with_carry_engine (T: integral) (P: {
  -- | Word size: number of bits in each word of the state sequence.
  -- Should be positive and less than the number of bits in `T.t`.
  val w: i32
  -- | Long lag: distance between operand values.
  val r: i32
  -- | Short lag: number of elements between advances.  Should be
  -- positive and less than `r`.
  val s: i32
}): rng_engine with int.t = T.t = {
  let long_lag = P.r
  let word_size = P.w
  let short_lag = P.s
  let modulus = T.i32 (1 << word_size)

  -- We use this one for initialisation.
  module e = linear_congruential_engine T {
    let a = T.u32 40014u32
    let c = T.u32 0u32
    let m = T.u32 2147483563u32
  }

  module int = T
  type t = T.t
  type rng = {x: [P.r]T.t,
              carry: bool,
              k: i32}

  let rand ({x, carry, k}: rng): (rng, t) =
    unsafe
    let short_index = k - short_lag
    let short_index = if short_index < 0
                      then short_index + long_lag
                      else short_index
    let (xi, carry) =
      if T.(x[short_index] >= x[k] + bool carry)
      then (T.(x[short_index] - x[k] - bool carry),
            false)
      else (T.(modulus - x[k] - bool carry + x[short_index]),
            true)
    let x = (copy x) with [k] = xi
    let k = (k + 1) % long_lag
    in ({x, carry, k}, xi)

  let rng_from_seed [n] (seed: [n]i32): rng =
    let rng = e.rng_from_seed seed
    let (x, _) = loop (x, rng) = (replicate P.r (T.i32 0), rng)
                   for i < P.r do let (v, rng) = e.rand rng
                                  in (x with [i] = T.(v % modulus),
                                      rng)
    let carry = T.(last x == i32 0)
    let k = 0
    in {x, carry, k}

  let split_rng (n: i32) ({x, carry, k}: rng): [n]rng =
    map (\i -> {x=map (T.^(T.i32 (hash i))) x, carry, k}) (iota n)

  let join_rng [n] (xs: [n]rng): rng =
    xs[0] -- FIXME

  let min = T.i32 0
  let max = T.(modulus - i32 1)
}

-- | An engine adaptor parametric module that adapts a pseudo-random
-- number generator Engine type by using only `r` elements of
-- each block of `p` elements from the sequence it produces,
-- discarding the rest.
--
-- The adaptor keeps and internal counter of how many elements have
-- been produced in the current block.
module discard_block_engine (K: {
  -- | Block size: number of elements in each block.  Must be
  -- positive.
  val p: i32
  -- | Used block: number of elements in the block that are used (not
  -- discarded). The rest `p-r` are discarded. This parameter should
  -- be greater than zero and lower than or equal to `p`.
  val r: i32}) (E: rng_engine): rng_engine with int.t = E.int.t = {
  type t = E.int.t
  module int = E.int
  type rng = (E.rng, i32)

  let min = E.min
  let max = E.max

  let rng_from_seed (xs: []i32): rng =
    (E.rng_from_seed xs, 0)

  let split_rng (n: i32) ((rng, i): rng): [n]rng =
    map (\rng' -> (rng', i)) (E.split_rng n rng)

  let join_rng (rngs: []rng): rng =
    let (rngs', is) = unzip rngs
    in (E.join_rng rngs', reduce i32.max 0 is)

  let rand ((rng,i): rng): (rng, t) =
    let (rng, i) =
      if i >= K.r then (loop rng for _j < K.r - i do (E.rand rng).1, 0)
                  else (rng, i+1)
    let (rng, x) = E.rand rng
    in ((rng, i), x)
}

-- | An engine adaptor that adapts an `rng_engine` so that the
-- elements are delivered in a different sequence.
--
-- The RNG keeps a buffer of `k` generated numbers internally, and
-- when requested, returns a randomly selected number within the
-- buffer, replacing it with a value obtained from its base engine.
module shuffle_order_engine (K: {val k: i32}) (E: rng_engine)
                          : rng_engine with int.t = E.int.t = {
  type t = E.int.t
  module int = E.int
  type rng = (E.rng, [K.k]t)

  let build_table (rng: E.rng) =
    let xs = replicate K.k (int.i32 0)
    in loop (rng,xs) for i < K.k do
         let (rng,x) = E.rand rng
         in (rng, xs with [i] = x)

  let rng_from_seed (xs: []i32) =
    build_table (E.rng_from_seed xs)

  let split_rng (n: i32) ((rng, _): rng): [n]rng =
    map build_table (E.split_rng n rng)

  let join_rng (rngs: []rng) =
    let (rngs', _) = unzip rngs
    in build_table (E.join_rng rngs')

  let rand ((rng,table): rng): (rng, int.t) =
    let (rng,x) = E.rand rng
    let i = i32.i64 (int.to_i64 x) % K.k
    let (rng,y) = E.rand rng
    in unsafe ((rng, (copy table) with [i] = y), table[i])

  let min = E.min
  let max = E.max
}

-- | A `linear_congruential_engine`@term producing `u32` values and
-- initialised with `a=48271`, `c=0` and
-- `m=2147483647`.  This is the same configuration as in C++.
module minstd_rand: rng_engine with int.t = u32 =
  linear_congruential_engine u32 {
    let a = 48271u32
    let c = 0u32
    let m = 2147483647u32
}

-- | A `linear_congruential_engine`@term producing `u32` values and
-- initialised with `a=16807`, `c=0` and
-- `m=2147483647`.  This is the same configuration as in C++.
module minstd_rand0: rng_engine with int.t = u32 =
  linear_congruential_engine u32 {
    let a = 16807u32
    let c = 0u32
    let m = 2147483647u32
}

-- | A subtract-with-carry pseudo-random generator of 24-bit numbers,
-- generally used as the base engine for the `ranlux24`@term generator.
-- It is an instantiation of `subtract_with_carry_engine`@term with
-- `w=24`, `s=10`, `r=24`.
module ranlux24_base: rng_engine with int.t = u32 =
  subtract_with_carry_engine u32 {
    let w:i32 = 24
    let s:i32 = 10
    let r:i32 = 24
  }

-- | A subtract-with-carry pseudo-random generator of 48-bit numbers,
-- generally used as the base engine for the `ranlux48`@term generator.
-- It is an instantiation of `subtract_with_carry_engine`@term with
-- `w=48`, `s=5`, `r=12`.
module ranlux48_base: rng_engine with int.t = u64 =
  subtract_with_carry_engine u64 {
    let w:i32 = 48
    let s:i32 = 5
    let r:i32 = 12
  }

-- | A subtract-with-carry pseudo-random generator of 24-bit numbers
-- with accelerated advancement.
--
-- It is an instantiation of a `discard_block_engine`@term with
-- `ranlux24_base`@term, with parameters `p=223` and `r=23`.
module ranlux24: rng_engine with int.t = u32 =
  discard_block_engine {let p:i32 = 223 let r:i32 = 23} ranlux24_base

-- | A subtract-with-carry pseudo-random generator of 48-bit numbers
-- with accelerated advancement.
--
-- It is an instantiation of a `discard_block_engine`@term with
-- `ranlux48_base`@term, with parameters `p=223` and `r=23`.
module ranlux48: rng_engine with int.t = u64 =
  discard_block_engine {let p:i32 = 389 let r:i32 = 11} ranlux48_base

-- | An engine adaptor that returns shuffled sequences generated with
-- `minstd_rand0`@term.  It is not a good idea to use this RNG in a
-- parallel setting, as the state size is fairly large.
module knuth_b: rng_engine with int.t = u32 =
  shuffle_order_engine {let k:i32 = 256} minstd_rand0

-- | The [xorshift128+](https://en.wikipedia.org/wiki/Xorshift#xorshift+) engine.  Uses
-- two 64-bit words as state.
module xorshift128plus: rng_engine with int.t = u64 = {
  module int = u64
  type rng = (u64,u64)

  let rand ((x,y): rng): (rng, u64) =
    let x = x ^ (x << 23u64)
    let new_x = y
    let new_y = x ^ y ^ (x >> 17u64) ^ (y >> 26u64)
    in ((new_x,new_y), (new_y + y))

  -- This seeding is quite a hack to ensure that we get good results
  -- even for poor seeds.  The main trick is to run a couple of rounds
  -- of the RNG after we're done.
  let rng_from_seed [n] (seed: [n]i32) =
    (loop (a,b) = (u64.i32 (hash (-n)), u64.i32 (hash n)) for i < n do
       if i % 2 == 0
       then (rand (a^u64.i32 (hash seed[i]),b)).1
       else (rand (a, b^u64.i32 (hash seed[i]))).1)
    |> rand |> (.1) |> rand |> (.1)

  let split_rng (n: i32) ((x,y): rng): [n]rng =
    map (\i -> let (a,b) = (rand (rng_from_seed [hash (i^n)])).1
               in (rand (rand (x^a,y^b)).1).1) (iota n)

  let join_rng [n] (xs: [n]rng): rng =
    reduce (\(x1,y1) (x2,y2) -> (x1^x2,y1^y2)) (0u64,0u64) xs

  let min = u64.lowest
  let max = u64.highest
}


-- | [PCG32](http://www.pcg-random.org/).  Has a state space of 128
-- bits, and produces uniformly distributed 32-bit integers.
module pcg32: rng_engine with int.t = u32 = {
  module int = u32
  type rng = {state: u64, inc: u64}

  let rand ({state, inc}: rng) =
    let oldstate = state
    let state = oldstate * 6364136223846793005u64 + (inc|1u64)
    let xorshifted = u32.u64 (((oldstate >> 18u64) ^ oldstate) >> 27u64)
    let rot = u32.u64 (oldstate >> 59u64)
    in ({state, inc},
        (xorshifted >> rot) | (xorshifted << ((-rot) & 31u32)))

  let rng_from_seed (xs: []i32) =
    let initseq = 0xda3e39cb94b95bdbu64 -- Should expose this somehow.
    let state = 0u64
    let inc = (initseq << 1u64) | 1u64
    let {state, inc} = (rand {state, inc}).1
    let state = loop state for x in xs do state + u64.i32 x
    in (rand {state, inc}).1

  let split_rng (n: i32) ({state,inc}: rng): [n]rng =
    map (\i -> (rand {state = state * u64.i32 (hash (i^n)), inc}).1) (iota n)

  let join_rng (rngs: []rng) =
    let states = map (\(x: rng) -> x.state) rngs
    let incs = map (\(x: rng) -> x.inc) rngs
    let state = reduce (*) 1u64 states
    let inc = reduce (|) 0u64 incs
    in {state, inc}

  let min = 0u32
  let max = 0xFFFFFFFFu32
}

-- | This uniform integer distribution generates integers in a given
-- range with equal probability for each, assuming the passed
-- `rng_engine`@mtype generates uniformly distributed integers.
module uniform_int_distribution (D: integral) (E: rng_engine):
  rng_distribution with num.t = D.t
                   with engine.rng = E.rng
                   with distribution = (D.t,D.t) = {

  let to_D (x: E.int.t) = D.i64 (E.int.to_i64 x)
  let to_E (x: D.t) = E.int.i64 (D.to_i64 x)

  module engine = E
  module num = D
  type distribution = (D.t,D.t) -- Lower and upper bounds.
  let uniform (min: D.t) (max: D.t) = (min,max)

  open E.int

  let rand ((min,max): distribution) (rng: E.rng) =
    let min = to_E min
    let max = to_E max
    let range = max - min + i32 1
    in if range <= i32 0
       then (rng, to_D E.min)
       else -- Avoid infinite loop if range exceeds what the RNG
            -- engine can supply.  This does not mean that we actually
            -- deliver sensible values, though.
            let secure_max = E.max - E.max %% range
            let (rng,x) = loop (rng, x) = E.rand rng
                          while x >= secure_max do E.rand rng
            in (rng, to_D (min + x / (secure_max / range)))
}

-- | This uniform integer distribution generates floats in a given
-- range with "equal" probability for each.
module uniform_real_distribution (R: real) (E: rng_engine):
  rng_distribution with num.t = R.t
                   with engine.rng = E.rng
                   with distribution = (R.t,R.t) = {
  let to_R (x: E.int.t) =
    R.u64 (u64.i64 (E.int.to_i64 x))

  module engine = E
  module num = R
  type distribution = (num.t, num.t) -- Lower and upper bounds.

  let uniform (min: num.t) (max: num.t) = (min, max)

  let rand ((min_r,max_r): distribution) (rng: E.rng) =
    let (rng', x) = E.rand rng
    let x' = R.((to_R x - to_R E.min) / (to_R E.max - to_R E.min))
    in (rng', R.(min_r + x' * (max_r - min_r)))
}

-- | Normally distributed floats.
module normal_distribution (R: real) (E: rng_engine):
  rng_distribution with num.t = R.t
                   with engine.rng = E.rng
                   with distribution = {mean:R.t,stddev:R.t} = {
  let to_R (x: E.int.t) =
    R.u64 (u64.i64 (E.int.to_i64 x))

  module engine = E
  module num = R
  type distribution = {mean: num.t, stddev: num.t}

  let normal (mean: num.t) (stddev: num.t) = {mean=mean, stddev=stddev}

  open R

  let rand ({mean,stddev}: distribution) (rng: E.rng) =
    -- Box-Muller where we only use one of the generated points.
    let (rng, u1) = E.rand rng
    let (rng, u2) = E.rand rng
    let u1 = (to_R u1 - to_R E.min) / (to_R E.max - to_R E.min)
    let u2 = (to_R u2 - to_R E.min) / (to_R E.max - to_R E.min)
    let r = sqrt (i32 (-2) * log u1)
    let theta = i32 2 * pi * u2
    in (rng, mean + stddev * (r * cos theta))
}
