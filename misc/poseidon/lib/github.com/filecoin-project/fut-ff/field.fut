import "../../athas/vector/vector"

module vector_2 = cat_vector vector_1 vector_1
module vector_3 = cat_vector vector_2 vector_1
module vector_4 = cat_vector vector_2 vector_2
module vector_5 = cat_vector vector_2 vector_3
module vector_8 = cat_vector vector_5 vector_3
module vector_9 = cat_vector vector_5 vector_4
module vector_11 = cat_vector vector_3 vector_8

-- https://github.com/filecoin-project/pairing/blob/master/src/bls12_381/fr.rs#L9
-- #[PrimeFieldModulus = "52435875175126190479447740508185965837690552500527637822603658699938581184513"]
let bls12_381_modulus  = "52435875175126190479447740508185965837690552500527637822603658699938581184513"
let r_squared_mod_p = "3294906474794265442129797520630710739278575682199800681788903916070560242797"

module type field_params = {
  module LVec: vector -- static vector to hold the limbs

  val limbs : i64  -- Number of limbs
  val p: () -> []u8 -- Size of prime field
  val r2: () -> []u8 -- R^2 mod p
}

module type fieldtype = {
  module FT: integral
  module LVec: vector
  type t

  val limbs: i64 -- number of limbs
  val double_limbs: i64 -- 2 x number of limbs
  val num_bits: i32
  val zero: t
  val one: t
  val highest: t

  -- P, the prime field modulus.
  val p_str: () -> []u8

  -- R^2 mod P
  -- It would be nice to be able to calculate this, but we need it in order to make multiplication in the field cheap.
  -- For now, we will settle for verifying it once montogmery multiplication reduction is bootstrapped. (see: field_r2_is_correct)
  val r2_str: () -> []u8

  val equal: t -> t -> bool
  val gt: t -> t -> bool
  val gte: t -> t -> bool
  val lt: t -> t -> bool
  val lte: t -> t -> bool
  val add: t -> t -> t
  val mul: t -> t -> t
  val sub: t -> t -> t
  val from_u8: u8 -> t
  val from_u16: u16 -> t
  val from_u32: u32 -> t
  val from_u64: u64 -> t
  val to_u64: t -> u64
  val mad_hi: t -> t -> t -> t

  val +: t -> t -> t
  val -: t -> t -> t
  val *: t -> t -> t

  val ==: t -> t -> bool
  val >=: t -> t -> bool
  val >: t -> t -> bool
  val <: t -> t -> bool
  val <=: t -> t -> bool

  val >>: t -> i32 -> t
  val <<: t -> i32 -> t
  val |: t -> t -> t
}

module fieldtype (T: integral) (P: field_params): fieldtype = {
  module FT = T
  module LVec = P.LVec
  type t = T.t

  let limbs =  P.limbs
  let double_limbs = 2 * P.limbs
  let num_bits = T.num_bits
  let zero = T.highest T.- T.highest
  let one = T.highest T./ T.highest
  let highest = T.highest
  let p_str = P.p
  let r2_str = P.r2

  let equal (a: t) (b: t): bool = a T.== b
  let gt (a: t) (b: t): bool = a T.> b
  let gte (a: t) (b: t): bool = a T.>= b
  let lt (a: t) (b: t): bool = a T.< b
  let lte (a: t) (b: t): bool = a T.<= b
  let add (a: t) (b: t): t = a T.+ b
  let sub (a: t) (b: t): t = a T.- b
  let mul (a: t) (b: t): t = a T.* b
  let from_u8 (n: u8): t = T.u8 n
  let from_u16 (n: u16): t = T.u16 n
  let from_u32 (n: u32): t = T.u32 n
  let from_u64 (n: u64): t = T.u64 n
  let to_u64 (n: t): u64 = u64.i64 (T.to_i64 n)
  let mad_hi (a: t) (b: t) (c: t): t = T.mad_hi a b c

  let (a: t) + (b: t) = a T.+ b
  let (a: t) - (b: t) = a T.- b
  let (a: t) * (b: t) = a T.* b
  let (a: t) == (b: t) = a T.== b
  let (a: t) >= (b: t) = a T.>= b
  let (a: t) > (b: t) = a T.> b
  let (a: t) <= (b: t) = a T.<= b
  let (a: t) < (b: t) = a T.< b
  let (x: t) >> (n: i32) = x T.>> (T.i32 n)
  let (x: t) << (n: i32) = x T.<< (T.i32 n)
  let (a: t) | (b: t) = a T.| b
}

type^ field_core 't 's = ((u8->t), t, t->t->t, t->t->t, u8->s)

-- In a perfect eventual world, this would include the integral module entirely.
module type field = {
  type t
  type s -- limb type
  type es -- element string
  type double_t

  module S: integral -- module for limb type

  val zero: t
  val one: t
  val highest: t
  val fill: s -> i32 -> t

  -- Normal operations
  val add: *t -> t -> t
  val sub: *t -> t -> t

  val ==: t -> t -> bool
  val >=: t -> t -> bool
  val >: t -> t -> bool
  val <=: t -> t -> bool
  val <: t -> t -> bool

  val big_mul: t -> t -> double_t
  val simple_reduce: double_t ->  t

  -- Field operations
  val +: t -> t -> t
  val -: t -> t -> t
  val *: t -> t -> t

  val to_mont: t -> t
  val mont_reduce: *double_t -> t
  val mont_field_mul: t -> t -> t
  val final_reduce: t -> t

  val double: t -> t
  val square: t -> t
  val pow: t -> i32 -> t
  val pow_lookup: t -> i32 -> t

  -- These from_ functions convert from the named integral type and populate the first limb.
  -- If the value must be truncated, that will happen.
  val from_u8: u8 -> t
  val from_u16: u16 -> t
  val from_u32: u32 -> t
  val from_u64: u64 -> t

  -- These mont_from_ functions are like their corresponding from_ functions
  -- but also convert the result to Montgomery representation.
  val mont_from_u8: u8 -> t
  val mont_from_u16: u16 -> t
  val mont_from_u32: u32 -> t
  val mont_from_u64: u64 -> t

  val s_from_u8: u8 -> s
  val s_from_u16: u16 -> s
  val s_from_u32: u32 -> s
  val s_from_u64: u64 -> s

  val mac_with_carry: s -> s -> s -> s -> (s, s)
  val mac_with_carry8: u8 -> u8 -> u8 -> u8 -> (s, s)
  val add_with_carry: s -> s -> (s, s)
  val add2_with_carry: s -> s -> s -> (s, s)

  val LIMBS: i64
  val FIELD_P: t
  val FIELD_R2: t
  val FIELD_P_DIFF: t
  val FIELD_INV: s

  val p_is_small: bool

  val from_string [n]: [n]u8 -> t
--  val make: (() -> []s) -> t

  val from_u64s [n]: [n]u64 -> t
  val to_u64s: t -> [LIMBS]u64
  val in_field: t -> bool

  val mont_from_u64s [n]: [n]u64 -> t
  val mont_to_u64s: t -> [LIMBS]u64

  -- Debugging
  val double_in_field: double_t -> t -> bool
  val double_sub: double_t -> double_t -> double_t
  val DOUBLE_FIELD_P: double_t
  val R_MOD_P: t
  val double_zero: double_t
  val DOUBLE_R: double_t

}

module big_field (M: fieldtype): field = {
  module LVec = M.LVec
  module LVec2 = cat_vector LVec LVec
  let DOUBLE_LIMBS = 2 * M.limbs
  let LIMBS = assert (LVec.length == M.limbs) M.limbs
  module S = M.FT

  type t = LVec.vector M.t --[LIMBS]M.t
  type s = M.t
  type double_t = LVec2.vector s --[DOUBLE_LIMBS]s

  let digit_count: i64 =
    let digits_per_byte = f64.log10(2 ** 8) in
    let bytes = LIMBS * (i64.i32 M.num_bits / 8) in
    i64.f64 (f64.ceil ((f64.i64 bytes) * digits_per_byte))

  type es = [digit_count]u8

  let limbs_are_unsigned = assert (M.lt M.zero (M.sub M.zero M.one)) true

  let zero: t = LVec.replicate M.zero --map (\_ -> M.zero) (iota LIMBS)
  let double_zero: double_t = LVec2.replicate M.zero --double_t = map (\_ -> M.zero) (iota DOUBLE_LIMBS)
  let DOUBLE_R: double_t = LVec2.map (\i -> if i == LIMBS then M.one else M.zero) LVec2.iota
  let one: t = LVec.map (\x -> if x == 0 then M.one else M.zero) LVec.iota

  let highest: t = LVec.replicate M.highest
  let fill (v: s) (count: i32) : t = LVec.map (\x -> if i32.i64 x < count then v else M.zero) LVec.iota

  let double_highest: double_t = LVec2.map (\i -> if i < LIMBS then M.highest else M.zero) LVec2.iota

  -- Do this before redefining ==.
  let no_p_str = length (M.p_str ()) == 0

  -- Primality is unchecked, but distinguish between the bignum and prime field cases. (TODO: separate the modules).
  let is_prime_field = assert (!no_p_str) true

  let (a: t) == (b: t) : bool = LVec.reduce (&&) true (LVec.map (uncurry M.equal) (LVec.zip a b))

  let bool_t (cond: bool): s = if cond then M.one else M.zero

  let (a: t) < (b: t) : bool =
    let op (a', b') (done, res) = if done then (true, res) else
                                  if M.(a' < b') then (true, true) else
                                  if M.(a' == b') then (false, false) else
                                    (true, false)
    in
    let (_, res) = LVec.foldr op (false, false) (LVec.zip a b)
    in res

  let (a: t) >= (b: t) : bool = !(a < b)
  let (a: t) > (b: t) : bool =  b < a
  let (a: t) <= (b: t) : bool = !(b < a)

  -- add is commutative because the two values (op and carry) which depend on either a or b
  -- are unaffected if a and b are reversed.
  let add (a: t) (b: t): t =
    let r: *t = copy zero in
    let carry = M.zero in
    -- op is commutative:
    let op (carry, r) (i, (a', b')) =
      -- a' and b' can be switched without affecting new because
      -- M.+ is defined to be + of a module of type integral, which is commutative.
      let new = M.(a' + b' + carry) in
      -- Note that the carry check could also have used b'.
      -- If a' + b' overflows then both a' and b' are greater than new:
      -- an overflowed sum cannot be greater than either of of its summands.
      -- If a' + b' does not overflow, then neither a' nor b' is greater than new;
      -- by definition, a non-overlowing sum is greater than either of its summands.
      -- Therefore, the choice of a' (vs b') is arbitrary and does not affect the result.
      let carry = bool_t M.(a' > new) in
      (carry, LVec.set i new r)
    -- Because op is commutative, the result is oblivious to the ordering of a and b when zipping.
    let (_, r) = LVec.foldl op (carry, r) (LVec.zip LVec.iota (LVec.zip a b))
    in r

  let sub (a: t) (b: t) : t =
    let op (borrow, r) (i, (a', b')) =
      let new = M.(a' - (b' + borrow)) in
      let borrow = bool_t M.(if borrow > zero then a' <= new else a' < new) in
      (borrow, LVec.set i new r) in
    let (_, r) = LVec.foldl op (M.zero, zero) (LVec.zip LVec.iota (LVec.zip a b))
    in r

   --  TODO: try to make this generic over size.
   let double_sub (a: double_t) (b: double_t) : double_t =
    let op (borrow, r) (i, (a', b')) =
      let new = M.(a' - (b' + borrow)) in
      let borrow = bool_t M.(if borrow > zero then a' <= new else a' < new) in
      (borrow, LVec2.set i new r) in
    let (_, r) = LVec2.foldl op (M.zero, double_zero) (LVec2.zip LVec2.iota (LVec2.zip a b))
    in r

  let mac_with_carry (a: s) (b: s) (c: s) (d: s): (s, s) =
    let lo = M.(a * b + c) in
    let hi =  M.(mad_hi a b (bool_t (lo < c))) in
    let a = lo in
    let lo = M.(lo + d) in
    let hi = M.(hi + (bool_t (lo < a))) in
    let d = hi in
    (lo, d)

    -- For testing
    let mac_with_carry8 (a: u8) (b: u8) (c: u8) (d: u8): (s, s) =
      M.(mac_with_carry (from_u8 a) (from_u8 b) (from_u8 c) (from_u8 d))

  let add_with_carry (a: s) (b: s): (s, s) =
    let lo = M.(a + b) in
    let b = bool_t M.(lo < a) in
    (lo, b)

  let add2_with_carry (a: s) (b: s) (c: s): (s, s) =
    let lo = M.(a + b) in
    let hi = bool_t M.(lo < a) in
    let a = lo in
    let lo = M.(lo + c) in
    let hi = M.(hi + (bool_t (lo < a))) in
    let c = hi in
    (lo, c)

  let double_in_field (a: double_t) (f: t): bool =
      all (\i -> M.(LVec2.get i64.(i + limbs) a == zero)) (iota LIMBS) &&
      --(f == zero) || -- Why isn't this equivalent to the next line?
      if f == zero then true else
      let (acc, _) = loop (acc, i) = (true, LIMBS - 1) while acc && (i i64.>= 0) do
                       if M.(LVec2.get i a < (LVec.get i f)) then (true, -1) else
                         (M.(LVec2.get i a == (LVec.get i f)), i - 1) in
      acc

  -- So we can read the initial FIELD_P value from a string.
  -- Naive means does not respect modulus.
  let naive_reduce (to_reduce: double_t): t =
    let in_field = loop to_reduce while !(double_in_field to_reduce zero) do
      double_sub to_reduce (copy double_highest) in
    LVec.from_array ((map (\i -> LVec2.get i in_field) (iota LIMBS)) :> [LVec.length]M.t)

  let big_mul (a: t) (b: t): LVec2.vector M.t = --[DOUBLE_LIMBS]M.t =
    let outer_op outer (i, a') =
      let (inner, carry) =
        let op (inner, carry) (j, b') =
          let (sum, carry) = mac_with_carry a' b' (LVec2.get (i+j) inner) carry
          in (LVec2.set (i+j) sum inner, carry)
        in LVec.foldl op (outer, M.zero) (LVec.zip LVec.iota b)
      in LVec2.set (i+LIMBS) carry inner
    in LVec.foldl outer_op (LVec2.replicate M.zero) (LVec.zip LVec.iota a)

  let naive_mul  (a: t) (b: t) : t =
    naive_reduce (big_mul a b)

  let from_u8 (n: u8): t = fill (M.from_u8 n) 1
  let from_u16 (n: u16): t = fill (M.from_u16 n) 1
  let from_u32 (n: u32): t = fill (M.from_u32 n) 1
  let from_u64 (n: u64): t = fill (M.from_u64 n) 1
  let ten = (from_u8 10)

  -- FIXME: Don't expose this version. Only use for reading FIELD_P.
  -- FIXME: assert valid digits
  let from_string [n] (str: [n]u8): t =
    let parse_digit (c: u8): t = from_u8 u8.(c - '0') in
    loop acc = copy zero for c in str do
      if c u8.== ' ' then
        acc else
        add (naive_mul acc (copy ten)) (parse_digit c)

  -- zero means field is the size of the underlying bits.
  let FIELD_P = if no_p_str then zero else from_string (M.p_str ())
  let FIELD_R2 = if no_p_str then zero else from_string (M.r2_str ())

  -- Could also do this by checking that the most-significant bit is 0.
  -- Is p less than half of the field size?
  let p_is_small = (add FIELD_P FIELD_P) >= FIELD_P

  let DOUBLE_FIELD_P: double_t =
    let fp = copy FIELD_P in
    LVec2.map (\i -> if i i64.< LIMBS then LVec.get i fp else M.zero) LVec2.iota

  let FIELD_P_DIFF = sub zero FIELD_P

  let calc_inv (a: s): s =
    let inv = (loop inv = M.one for _i < M.num_bits do
               let inv = M.(inv * inv) in M.(inv * a))
    in M.(zero - inv)

  let FIELD_INV = calc_inv (LVec.get 0 FIELD_P)

  let simple_reduce (to_reduce: double_t): t =
    let dfp = (copy DOUBLE_FIELD_P) in
    let in_field = loop to_reduce while !(double_in_field to_reduce FIELD_P) do
         double_sub to_reduce dfp in
    LVec.map (\i -> LVec2.get i in_field) LVec.iota

  let simple_field_mul  (a: t) (b: t) : t =
    simple_reduce (big_mul a b)

  let mont_reduce (limbs: double_t): t =
    let (FIELD_P, FIELD_INV) = (copy (FIELD_P, FIELD_INV)) in
    let (lmbs, _) =
      let outer_op (lmbs, carry2) i =
        let u = M.(FIELD_INV * (LVec2.get i lmbs)) in
        let (lmbs, carry) =
          let op (lmbs, carry) (j, fp) =
            let (x, carry) = mac_with_carry u fp (LVec2.get (i+j) lmbs) carry in
            (LVec2.set (i+j) x lmbs, carry)
          in LVec.foldl op (lmbs, M.zero) (LVec.zip LVec.iota FIELD_P) in
        let (x, carry2) = add2_with_carry (LVec2.get (i+LIMBS) lmbs) carry carry2 in
        (LVec2.set (i+LIMBS) x lmbs, carry2) in
      LVec.foldl outer_op (limbs, M.zero) LVec.iota
    in let result = LVec.map (\i -> LVec2.get (i+LIMBS) lmbs) LVec.iota
       in if result >= FIELD_P then sub result FIELD_P else result

  let final_reduce (v: t): t =
    let double = LVec2.map (\i -> if i i64.< LIMBS then LVec.get i v else M.zero) LVec2.iota in
    mont_reduce double

  let R_MOD_P = if FIELD_P == zero then zero else simple_reduce DOUBLE_R

  let to_mont(v: t): t =
    mont_reduce (big_mul (copy FIELD_R2) v)

  let mont_field_mul  (a: t) (b: t) : t =
    mont_reduce (big_mul a b)


  let (a: t) * (b: t) : t =
    mont_field_mul a b

  -- Be paranoid and fail if the provided value for r^2 mod p is incorrect.
  let field_r2_is_correct = assert (final_reduce FIELD_R2 == R_MOD_P) true

  let (a: t) - (b: t): t =
    let old = copy a in
    let res = sub a b in
    if old >= b then res else add res (copy FIELD_P)

  -- If add_expensive has been called, we can assume that the field is not 'small',
  -- which means addition can potentially overflow the underlying bignum.
  let add_expensive (a: t) (b: t): t =
    -- Both inputs must be in field.
    let fzero = FIELD_P == zero in
    -- (Can move or remove this check eventually, as long as the invariant is otherwise enforced.)
    let _ = if fzero then () else assert ((a < FIELD_P) && (b < FIELD_P)) () in
    -- res is independent of the ordering of a and b if add is commutative (see add).
    let res = add a b in
    if fzero then res else -- special case for simple bignum (make own module? or at least abstract this check)
      -- if add really does implement unchecked bignum addition, then res >= a implies res >= b,
      -- so choice of a in the test below is irrelevant.
    if res >= a then -- Can we skip this check and instead check the carry? (would need to return it)
      -- Since res >= a (or b), add did not overflow the limbs.
      (if res >= FIELD_P then
         -- At most one subtraction effects the modulus.
         sub res (copy FIELD_P) else res) else
      -- Since the field is not 'small', res could be less than a (or b).
      -- In that case, the addition was enough to actually overflow the limbs.
      -- Then we need to compensate for the difference the field modulus and the underlying bignum capacity.
      -- This amount (FIELD_P_DIFF) needs to be added back in so the overflow is effectively relative to the modulus.
      add res (copy FIELD_P_DIFF)

  -- add_cheap is commutative (and must be so + will be).
  let add_cheap (a: t) (b: t): t =
    -- Both inputs must be in field.
    let fzero = FIELD_P == zero in
    -- (Can move or remove this check eventually, as long as the invariant is otherwise enforced.)
    let _ = if fzero then () else assert ((a < FIELD_P) && (b < FIELD_P)) () in

    -- add_cheap is commutative iff add is (see add).
    let res = add a b in
    if res >= FIELD_P then sub res (copy FIELD_P) else res

  -- + should be commutative for use with reduce_comm and associative for use with reduce.
  -- Note that for bls12_381, p_is_small is true — so, for that field in particular, it suffices
  -- to show that add_cheap is commutative.
  -- Note also that if + correctly implements adition, then it will be associative and commutative
  -- by definition, so further hand-waving is mostly superfluous.
  let (a: t) + (b: t) = if p_is_small then add_cheap a b else add_expensive a b

  let square (x: t): t =
    let _ = assert false "unimplemented" in
    -- This implementation is failing with large enough values, example test case:
    -- bls12_381.((from_string "12345678912345678912345678901") * (from_string "12345678912345678912345678901") == (square (from_string "12345678912345678912345678901")))
    -- should return true, but returns false.

    -- Long multiplication (Diagonal elements are skipped)
    let (res: double_t) = LVec2.replicate M.zero in
    let res =
      loop res for i < i64.(LIMBS-1) do
      let (res, carry, _) =
        loop (res, carry, j) = (res, M.zero, i64.(i + 1)) while i64.(j < LIMBS) do
         let (sum, carry) = mac_with_carry (LVec.get i x) (LVec.get j x) (LVec2.get i64.(i+j) res) carry
         in (LVec2.set i64.(i+j) sum res, carry, i64.(j+1))
      in LVec2.set i64.(i + LIMBS) carry res in

    --  Double the result
    let res = LVec2.set i64.(LIMBS * 2 - 1) ((copy (LVec2.get i64.(LIMBS * 2 - 2) res)) M.>> i32.(M.num_bits - 1)) res in
    let (res, _) =
      loop (res, i) = (res, i64.(LIMBS * 2 - 2)) while i64.(i > 1) do
      let res = LVec2.set i M.((copy (LVec2.get i res) << 1) | ((copy (LVec2.get i64.(i - 1) res)) >> i32.(M.num_bits - 1))) res in
      (res, i64.(i - 1))
    in let res = LVec2.set 1 M.(copy (LVec2.get 1 res) << 1) res in

    --  Process diagonal elements
    let (_, _, res) =
      loop (i, carry, res) = (0, M.zero, res) while i64.(i < LIMBS) do
      let (a, carry) = mac_with_carry (LVec.get i x) (LVec.get i x) (copy (LVec2.get i64.(i * 2) res)) carry in
      let (b, carry) = add_with_carry (copy (LVec2.get i64.(i * 2 + 1) res)) carry in
      (i i64.+ 1,
       carry,
       (LVec2.set i64.(i * 2 + 1) b (LVec2.set i64.(i * 2) a res)))

    in mont_reduce res

  let double (_x: t): t = assert false (copy zero) -- TODO: implement
  let pow (_base: t) (_exp: i32): t = assert false (copy zero) -- TODO: implement
  let pow_lookup (_base: t) (_exp: i32): t = assert false (copy zero) -- FIXME: implement

  let mont_u8 (n: u8): t = to_mont (from_u8 n)
  let s_from_u8 (n: u8): s = M.from_u8 n
  let s_from_u16 (n: u16): s = M.from_u16 n
  let s_from_u32 (n: u32): s = M.from_u32 n
  let s_from_u64 (n: u64): s = M.from_u64 n

  let from_u64s [n] (limbs: [n]u64): t =
    assert i64.((length limbs) == LIMBS) (LVec.map M.from_u64 (LVec.from_array (limbs :> [LVec.length]u64)) :> t)

  -- For use in tests, return true if its argument is actually in the field.
  let in_field (x: t): bool =
    (x < FIELD_P)

  let to_u64s (limbs: t): [LIMBS]u64 =
    LVec.to_array (LVec.map M.to_u64 limbs) :> [LIMBS]u64

let mont_from_u8 (n: u8): t = to_mont (fill (M.from_u8 n) 1)
  let mont_from_u16 (n: u16): t = to_mont (fill (M.from_u16 n) 1)
  let mont_from_u32 (n: u32): t = to_mont (fill (M.from_u32 n) 1)
  let mont_from_u64 (n: u64): t = to_mont (fill (M.from_u64 n) 1)

let mont_from_u64s [n] (limbs: [n]u64): t =
  assert i64.((length limbs) == LIMBS) (to_mont (LVec.map M.from_u64 (LVec.from_array (limbs :> [LVec.length]u64)) :> t))

  let mont_to_u64s (limbs: t): [LIMBS]u64 =
    LVec.to_array (LVec.map M.to_u64 (final_reduce limbs)) :> [LIMBS]u64

}

-- Non-prime fields don't work with montgomery representation. TODO: support them separately.

-- module b32: field = big_field (fieldtype u8 { let limbs = 4i32
--                                               let p () = ""
--                                               let r2 () = "" })

-- module b64: field = big_field (fieldtype u16 { let limbs = 4i32
--                                                let p () = ""
--                                                let r2 () = "" })

-- module b24: field = big_field (fieldtype u8 { let limbs = 3i32
--                                               let p () = ""
--                                               let r2 () = "" })

-- module b256: field = big_field (fieldtype u64 { let limbs = 4i32
--                                                 let p () = ""
--                                                 let r2 () = "" })

-- module b8: field = big_field (fieldtype u8 { let limbs = 1i32
--                                              let p () = ""
--                                              let r2 () = "" })

module b8': field = big_field (fieldtype u8 { module LVec = vector_1
                                              let limbs = 1i64
                                              let p () = "251"
                                              let r2 () = "25" })
-- -- module b16: field = big_field (fieldtype u8 { let limbs = 2i32
-- --                                               let p () = ""
-- --                                               let r2 () = "" })

module bls12_381: field = big_field (fieldtype u64 {
                                                 module LVec = vector_4
                                                 let limbs = 4i64
                                                 let p () = copy bls12_381_modulus
                                                 let r2 () = copy r_squared_mod_p })
