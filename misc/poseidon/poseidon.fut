module F = import "lib/github.com/filecoin-project/fut-ff/field"

import "lib/github.com/athas/vector/vector"
module vector_2 = cat_vector vector_1 vector_1
module vector_3 = cat_vector vector_2 vector_1
module vector_4 = cat_vector vector_2 vector_2

-- For some reason, we can't access these bls12_381-related definitison from fut-ff/field, so just repeat here.
-- https://github.com/filecoin-project/pairing/blob/master/src/bls12_381/fr.rs#L9
-- #[PrimeFieldModulus = "52435875175126190479447740508185965837690552500527637822603658699938581184513"]
let bls12_381_modulus  = "52435875175126190479447740508185965837690552500527637822603658699938581184513"
let r_squared_mod_p =  "3294906474794265442129797520630710739278575682199800681788903916070560242797"

module bls12_381: F.field = F.big_field (F.fieldtype u64 {
                                                       module LVec = vector_4
                                                       let limbs = 4i64
                                                       let p () = copy bls12_381_modulus
                                                       let r2 () = copy r_squared_mod_p })

type option 't = #Some(t) | #None

module type Params = {
  val arity: i64
  val full_rounds: i64
  val partial_rounds: i64
}

type matrix 't [width] = [width][width] t

type sparse_matrix 't [width] [width_] = {
    w_hat: [width]t,
    v_rest: [width_]t
}

type sparse_matrix_in 't [n] = [n]t


type constants 't [width] [width_] [rk_count] [sparse_count]= {
    arity_tag: t,
    round_keys: [rk_count]t,
    mds_matrix: matrix t [width],
    pre_sparse_matrix: matrix t [width],
    sparse_matrixes: [sparse_count]sparse_matrix t [width] [width_]
}

let div_evenly64 (a: i64) (b: i64): i64 = assert (a % b == 0) a /b
let div_evenly (a: i32) (b: i32): i32 = assert (a % b == 0) a /b

let even_chunks 't [n]  (chunk_size: i64) (arr: [n] t) =
  let chunk_count = div_evenly64 n chunk_size in
  unflatten chunk_count chunk_size arr


module type hasher  = {
  module Field: F.field

  type state

  val width: i64
  val arity: i64
  val width_: i64
  val rk_count: i64
  val sparse_count: i64
  val sparse_matrix_size: i64

  val init: (constants: constants Field.t [width] [width_] [rk_count] [sparse_count]) -> state
  val blank_constants: constants Field.t [width] [width_] [rk_count] [sparse_count]

  val constants: state -> constants Field.t [width] [width_] [rk_count] [sparse_count]

  val set_preimage: state -> [arity]Field.t -> state
  val hash: state -> Field.t
  val hash_preimage: state -> [arity]Field.t -> Field.t
  val hash_preimages [n]: state -> [n][arity]Field.t -> [n]Field.t
  val hash_preimages_u64s [u64_count]: state -> [u64_count]u64 -> ([][Field.LIMBS]u64, state)
  val hash_preimages_monts [u64_count]: state -> [u64_count]u64 -> ([][Field.LIMBS]u64, state)

  val reset: state -> state

  val leaves_per_kib: i64 -> i64
  val leaves_per_mib: i64 -> i64
  val leaves_per_gib: i64 -> i64

  val make_constants: [Field.LIMBS]u64 -> [rk_count][Field.LIMBS]u64 -> matrix ([Field.LIMBS]u64) [width] -> matrix ([Field.LIMBS]u64) [width] -> [sparse_count][sparse_matrix_size][Field.LIMBS]u64 -> constants (Field.t) [width] [width_] [rk_count] [sparse_count]
}

module make_hasher (f: F.field) (p: Params): hasher = {
  module Field = f

  let arity = p.arity
  let width = arity + 1
  let width_ = width - 1 -- width *could* be other than 1 + arity.
  let full_rounds = p.full_rounds
  let partial_rounds = p.partial_rounds
  let rk_count = width * p.full_rounds + p.partial_rounds -- TODO: Check this. We'll find out if it's wrong when integrating, of course.
  let sparse_count = p.partial_rounds -- TODO: check
  let sparse_matrix_size = width + width_

  let full_half = div_evenly64 full_rounds 2
  let sparse_offset = full_half - 1

  type elements = [width]Field.t
  type mat = matrix Field.t [width]
  type s_mat = sparse_matrix Field.t [width] [width_]

  type state = {
      constants: constants Field.t [width] [width_] [rk_count] [sparse_count],
      elements: elements,
      current_round: i32,
      rk_offset: i32
  }

  let init (constants: constants Field.t [width] [width_] [rk_count] [sparse_count]): state = {
      constants,
      elements = map (\i -> if i == 0 then constants.arity_tag else Field.zero) (iota width),
      current_round = 0,
      rk_offset = 0
  }

  let mk_arity_tag (arity: i64): Field.t  = Field.from_u32 (u32.i64 (1 << arity))

  let blank_constants: constants Field.t [width] [width_] [rk_count] [sparse_count] =
    let dummy = Field.from_string "123" in -- dummy value for now, don't use 0 or 1, so we get hash-like results.
    let arity_tag = mk_arity_tag arity
    let round_keys = map (\_ -> dummy) (iota rk_count) :> [rk_count]Field.t in
    let mds_matrix = map (\_i -> map (\_j -> dummy) (iota width)) (iota width)
    let pre_sparse_matrix =  map (\_i -> map (\_j -> dummy) (iota width)) (iota width)
    let sparse_matrixes = map (\_ -> { w_hat = map (\_ -> dummy) (iota width),
                                       v_rest = map (\_ -> dummy) (iota width_) })
                              (iota sparse_count) in
    { arity_tag,
      round_keys,
      mds_matrix,
      pre_sparse_matrix,
      sparse_matrixes }

  let make_sparse_matrix (array: [sparse_matrix_size]Field.t): sparse_matrix Field.t [width] [width_] =
    { w_hat = take width array,
      v_rest = array |> drop width |> take width_ }
  let make_constants (arity_tag: ([Field.LIMBS]u64)) (round_keys: [rk_count]([Field.LIMBS]u64)) (mds_matrix: matrix ([Field.LIMBS]u64) [width]) (pre_sparse_matrix: matrix ([Field.LIMBS]u64) [width])
  (sparse_matrixes: [sparse_count][sparse_matrix_size][Field.LIMBS]u64): constants Field.t [width] [width_] [rk_count] [sparse_count] =
    let sparse_matrixes = map make_sparse_matrix (map (map Field.mont_from_u64s) sparse_matrixes) in
    { arity_tag = Field.mont_from_u64s arity_tag,
      round_keys = map Field.mont_from_u64s round_keys,
      mds_matrix = map (map Field.mont_from_u64s) mds_matrix,
      pre_sparse_matrix = map (map Field.mont_from_u64s) pre_sparse_matrix,
      sparse_matrixes }

  let reset (s: state): state = s with current_round = 0
                                  with rk_offset = 0
                                  with elements = map (\i -> if i == 0 then s.constants.arity_tag else Field.zero) (iota width)

  let constants (s: state): constants Field.t [width] [width_] [rk_count] [sparse_count] = s.constants

  let set_preimage (s: state) (preimage: [arity]Field.t): state =
    s with elements = (([s.constants.arity_tag] ++ preimage) :> elements)

  let result (s: state): Field.t = s.elements[1]

  -- Placeholder until Field.square is fixed.
  let naive_square (x: Field.t): Field.t = Field.(x * x)

  let quintic_s_box (x: Field.t): Field.t =
    let x2 = naive_square x in
    let x4 = naive_square x2
    in Field.(x4 * x)

  let add_round_key (s: state) (rk_offset: i32) (i: i32) =
    Field.(s.elements[i] + s.constants.round_keys[i32.(rk_offset + i)])

  -- Could be more generic, but could also use a library. Just target clarity.
  let scalar_product (a: elements) (b: elements) = Field.(reduce (+) zero (map2 (*) a b))

  let apply_matrix (m: mat) (elts: elements): elements =
    map (scalar_product elts) (transpose m)

  let apply_sparse_matrix (sm: s_mat) (elts: [width]Field.t): elements =
    let first = scalar_product elts sm.w_hat
    let rest = map2 (\x y -> Field.(x + y * elts[0]))
                    (elts[1:width] :> [width_]Field.t)
                    sm.v_rest
    in ([first] ++ rest) :> elements

  let apply_round_matrix (s: state): state =
    let elements = if s.current_round == i32.i64 sparse_offset then
                     apply_matrix s.constants.pre_sparse_matrix s.elements
                   else if (s.current_round > i32.i64 sparse_offset)
                           && (s.current_round < i32.i64 full_half + i32.i64 partial_rounds) then
                        let index = s.current_round - i32.i64 sparse_offset - 1 in
                        apply_sparse_matrix s.constants.sparse_matrixes[index] s.elements
                   else
                     apply_matrix s.constants.mds_matrix s.elements
    in s with elements = elements

  let add_full_round_keys (s: state): state =
     s with elements = map (\i -> add_round_key s s.rk_offset i)
                           (map i32.i64 (iota width))
       with rk_offset = s.rk_offset + i32.i64 width

  let add_partial_round_key (s: state): state =
    s with elements = map (\i -> if i == 0 then add_round_key s s.rk_offset i else s.elements[i])
                          (map i32.i64 (iota width))
      with rk_offset = s.rk_offset + 1

  let full_round (s: state): state =
    let elements = map quintic_s_box s.elements in
    let s = add_full_round_keys (s with elements = elements)
    in apply_round_matrix s
      with current_round = s.current_round + 1

  let last_full_round (s: state): state =
    let elements = map quintic_s_box s.elements in
    let s = s with elements = elements
              with current_round = s.current_round + 1
    in apply_round_matrix s

  let partial_round (s: state): state =
    let elements = copy s.elements with [0] = quintic_s_box s.elements[0] in
    let s = add_partial_round_key (s with elements = elements)
    in apply_round_matrix s with current_round = s.current_round + 1

  let hash (s: state): Field.t =
    let s = add_full_round_keys s in
    let s = loop s for _i < full_half do full_round s in
    let s = loop s for _i < partial_rounds do partial_round s in
    let s = loop s for _i < (full_half - 1) do full_round s in
    let s = last_full_round s
    in result s

  let hash_preimage (s: state) (preimage: [arity]Field.t) =
    hash (set_preimage s preimage)

  let hash_preimages [n] (s: state) (preimages: [n][arity]Field.t) =
    map (hash_preimage s) preimages

  let hash_preimages_monts [u64_count] (s: state) (preimages_u64s: [u64_count]u64) : ([][Field.LIMBS]u64, state) =
    let preimages = map Field.from_u64s (even_chunks Field.LIMBS preimages_u64s) in
    let hashes = (hash_preimages s) (even_chunks arity preimages) in
    let hashes_len = length hashes in
    let finalized_hashes = (map Field.to_u64s hashes) :> [hashes_len][Field.LIMBS]u64
    in (finalized_hashes, s)

  let hash_preimages_u64s [u64_count] (s: state) (preimages_u64s: [u64_count]u64) : ([][Field.LIMBS]u64, state) =
    let preimages = map Field.mont_from_u64s (even_chunks Field.LIMBS preimages_u64s) in
    let hashes = (hash_preimages s) (even_chunks arity preimages) in
    let hashes_len = length hashes in
    let finalized_hashes = (map Field.mont_to_u64s hashes) :> [hashes_len][Field.LIMBS]u64
    in (finalized_hashes, s)

  -- FIXME: this assumes 32-byte leaves, but we should actually get the element width (rounded up to nearest byte) from the hasher.
  let leaves_per_kib (kib: i64) =
    kib * 1024 / 32

  let leaves_per_mib (mib: i64) =
    mib * 1024 * 1024 / 32

  let leaves_per_gib (gib: i64) =
    gib * 1024 * 1024 * 1024 / 32
}

module type tree_builder = {
  module Hasher: hasher

  val leaves: i64
  val height: i64
  val tree_size: i64

  val build_tree:  Hasher.state -> [leaves]Hasher.Field.t -> [tree_size]Hasher.Field.t
  val compute_root: Hasher.state -> [leaves]Hasher.Field.t -> Hasher.Field.t

}

module type tree_builder_params = {
--  module Hasher: hasher
  val leaves: i64
}

module make_tree_builder (H: hasher) (P: tree_builder_params): tree_builder = {
  module Hasher = H
  module Field = Hasher.Field
  type t = Field.t

  let zero = Field.zero

  let leaves = P.leaves
  let arity = Hasher.arity

  let tree_dimensions (leaves: i64) (arity: i64): (i64, i64) =
    let (height, size, _) = loop (height, size, row_size) = (0, leaves, leaves) while row_size > 1 do
    let new_row_size = div_evenly64 row_size arity in
    (height + 1, size + new_row_size, new_row_size)

    in (height, size)

  let dimensions = tree_dimensions leaves arity
  let height = dimensions.0
  let tree_size = dimensions.1

  -- Like build_tree but does not retain intermediate results.
  let compute_root (s: Hasher.state) (base: [leaves]t): t =
    -- Base shrinks by a factor of 1/Hasher.arity on each iteration.
    let nodes = loop base while length base > 1 do
                  map  (\preimage -> (Hasher.hash_preimage s preimage))
                       (even_chunks Hasher.arity base)
    in assert (length nodes == 1) nodes[0]

  -- Like compute_root but returns an array of all rows, not just the last.
  -- This includes the original base row in the final treee, but to minimize memory, maybe we should define it not to.
  let build_tree (s: Hasher.state) (base: [leaves]t): [tree_size]t =
    -- Row shrinks by a factor of 1/Hasher.arity on each iteration.
    let (tree, _, _) =
      loop (tree, row, offset) = ((replicate tree_size Field.zero) with [0:(length base)] = base,
                                  base,
                                  length base)
      while length row > 1 do
      let preimages = (even_chunks Hasher.arity row) in
      let new_row = (Hasher.hash_preimages s) preimages in
      (tree with [offset:offset+length new_row] = new_row,
       new_row,
       offset + length new_row)
    in tree
}

module type column_tree_builder = {
  module ColumnHasher: hasher
  module TreeBuilder: tree_builder

  type state

  val column_size: i64 -- elements per column

  val init: ColumnHasher.state -> TreeBuilder.Hasher.state -> state
  val reset: state -> state
  val add_columns [elt_count]: state -> [elt_count]ColumnHasher.Field.t -> state
  val finalize: state -> ([TreeBuilder.tree_size][TreeBuilder.Hasher.Field.LIMBS]u64, state)
}

module make_column_tree_builder (ColumnHasher: hasher) (TreeBuilder: tree_builder): column_tree_builder = {
  module ColumnHasher = ColumnHasher
  module TreeBuilder = TreeBuilder

  type state = {
      column_state: ColumnHasher.state,
      tree_hasher_state: TreeBuilder.Hasher.state,
      leaves: [TreeBuilder.leaves]ColumnHasher.Field.t,
      pos: i64
  }

  let column_size = ColumnHasher.arity

  let init (column_state: ColumnHasher.state) (tree_hasher_state: TreeBuilder.Hasher.state): state =
    { column_state = copy column_state,
      tree_hasher_state = copy tree_hasher_state,
      leaves = replicate TreeBuilder.leaves ColumnHasher.Field.zero,
      pos = 0}

  let reset (s: state): state =
    s with column_state = ColumnHasher.reset s.column_state
      with tree_hasher_state = TreeBuilder.Hasher.reset s.tree_hasher_state
      with leaves = replicate TreeBuilder.leaves ColumnHasher.Field.zero
      with pos = 0

   let add_columns [elt_count] (s: state) (flat_columns: [elt_count]ColumnHasher.Field.t): state =
    let columns = even_chunks ColumnHasher.arity flat_columns in
    let column_leaves = (map (ColumnHasher.hash_preimage (copy s.column_state)) columns) in
    let new_pos = s.pos + length column_leaves in
    s with leaves = (copy s.leaves with [s.pos:new_pos] = column_leaves)
      with pos = new_pos


  let tree_size_in_limbs = TreeBuilder.tree_size * TreeBuilder.Hasher.Field.LIMBS

  let finalize (s: state): ([TreeBuilder.tree_size][TreeBuilder.Hasher.Field.LIMBS]u64, state) =
    let leaves = map (\i -> TreeBuilder.Hasher.Field.from_u64s (ColumnHasher.Field.to_u64s s.leaves[i])) (iota TreeBuilder.leaves) in
    let mont_tree = (TreeBuilder.build_tree s.tree_hasher_state leaves) in
    let output = map TreeBuilder.Hasher.Field.mont_to_u64s mont_tree
    in (output, s)
}

--------------------------------------------------------------------------------

module p2: hasher = make_hasher bls12_381 { let arity = 2i64
                                            let full_rounds = 8i64
                                            let partial_rounds = 55i64 }


module p4: hasher = make_hasher bls12_381 { let arity = 4i64
                                            let full_rounds = 8i64
                                            let partial_rounds = 56i64 }

module p8: hasher = make_hasher bls12_381 { let arity = 8i64
                                            let full_rounds = 8i64
                                            let partial_rounds = 57i64 }


-- Strengthened Poseidon: 1.25% more partial rounds.
module s2: hasher = make_hasher bls12_381 { let arity = 2i64
                                            let full_rounds = 8i64
                                            let partial_rounds = 69i64 }

module s4: hasher = make_hasher bls12_381 { let arity = 4i64
                                            let full_rounds = 8i64
                                            let partial_rounds = 70i64 }

module s8: hasher = make_hasher bls12_381 { let arity = 8i64
                                            let full_rounds = 8i64
                                            let partial_rounds = 72i64 }

module p11: hasher = make_hasher bls12_381 {
  let arity = 11i64
  let full_rounds = 8i64
  let partial_rounds = 57i64 }

module s11: hasher = make_hasher bls12_381 {
  let arity = 11i64
  let full_rounds = 8i64
  let partial_rounds = 72i64 }

module t2_2k: tree_builder =  make_tree_builder p2 { let leaves: i64 = p2.leaves_per_kib 2 }
module t4_2k: tree_builder =  make_tree_builder p4 { let leaves: i64 = p4.leaves_per_kib 2 }
module t8_2k: tree_builder =  make_tree_builder p8 { let leaves: i64 = p8.leaves_per_kib 2 }
module t8_64m: tree_builder =  make_tree_builder p8 { let leaves: i64 = p8.leaves_per_mib 64 }

module t8_512m: tree_builder =  make_tree_builder p8 { let leaves: i64 = p8.leaves_per_mib 512 }
module t8_4g: tree_builder =  make_tree_builder p8 { let leaves: i64 = p8.leaves_per_gib 4 }

type p2_state = p2.state

entry init2 (arity_tag: ([p2.Field.LIMBS]u64))
           (round_keys: [p2.rk_count]([p2.Field.LIMBS]u64))
           (mds_matrix: matrix ([p2.Field.LIMBS]u64) [p2.width])
           (pre_sparse_matrix: matrix ([p2.Field.LIMBS]u64) [p2.width])
           (sparse_matrixes: [p2.sparse_count][p2.sparse_matrix_size]([p2.Field.LIMBS]u64))
              : p2_state =
  let constants = p2.make_constants arity_tag round_keys mds_matrix pre_sparse_matrix sparse_matrixes in
  p2.init constants

type p8_state = p8.state

entry init8 (arity_tag: ([p8.Field.LIMBS]u64))
           (round_keys: [p8.rk_count]([p8.Field.LIMBS]u64))
           (mds_matrix: matrix ([p8.Field.LIMBS]u64) [p8.width])
           (pre_sparse_matrix: matrix ([p8.Field.LIMBS]u64) [p8.width])
           (sparse_matrixes: [p8.sparse_count][p8.sparse_matrix_size]([p8.Field.LIMBS]u64))
              : p8_state =
  let constants = p8.make_constants arity_tag round_keys mds_matrix pre_sparse_matrix sparse_matrixes in
  p8.init constants

entry hash8 (s: p8_state) (preimage_u64s: [32]u64) : ([p8.Field.LIMBS]u64, p8_state) =
  let preimage = map p8.Field.mont_from_u64s (even_chunks p8.Field.LIMBS preimage_u64s) in
  let hash = p8.Field.mont_to_u64s (p8.hash_preimage s (preimage :> [p8.arity]p8.Field.t)) in
  (hash, s)

type p11_state = p11.state

entry init11 (arity_tag: ([p11.Field.LIMBS]u64))
           (round_keys: [p11.rk_count]([p11.Field.LIMBS]u64))
           (mds_matrix: matrix ([p11.Field.LIMBS]u64) [p11.width])
           (pre_sparse_matrix: matrix ([p11.Field.LIMBS]u64) [p11.width])
           (sparse_matrixes: [p11.sparse_count][p11.sparse_matrix_size]([p11.Field.LIMBS]u64))
              : p11_state =
  let constants = p11.make_constants arity_tag round_keys mds_matrix pre_sparse_matrix sparse_matrixes in
  p11.init constants

type s2_state = s2.state

entry init2s (arity_tag: ([s2.Field.LIMBS]u64))
           (round_keys: [s2.rk_count]([s2.Field.LIMBS]u64))
           (mds_matrix: matrix ([s2.Field.LIMBS]u64) [s2.width])
           (pre_sparse_matrix: matrix ([s2.Field.LIMBS]u64) [s2.width])
           (sparse_matrixes: [s2.sparse_count][s2.sparse_matrix_size]([s2.Field.LIMBS]u64))
              : s2_state =
  let constants = s2.make_constants arity_tag round_keys mds_matrix pre_sparse_matrix sparse_matrixes in
  s2.init constants

type s8_state = s8.state

entry init8s (arity_tag: ([s8.Field.LIMBS]u64))
           (round_keys: [s8.rk_count]([s8.Field.LIMBS]u64))
           (mds_matrix: matrix ([s8.Field.LIMBS]u64) [s8.width])
           (pre_sparse_matrix: matrix ([s8.Field.LIMBS]u64) [s8.width])
           (sparse_matrixes: [s8.sparse_count][s8.sparse_matrix_size]([s8.Field.LIMBS]u64))
              : s8_state =
  let constants = s8.make_constants arity_tag round_keys mds_matrix pre_sparse_matrix sparse_matrixes in
  s8.init constants

type s11_state = s11.state

entry init11s (arity_tag: ([s11.Field.LIMBS]u64))
           (round_keys: [s11.rk_count]([s11.Field.LIMBS]u64))
           (mds_matrix: matrix ([s11.Field.LIMBS]u64) [s11.width])
           (pre_sparse_matrix: matrix ([s11.Field.LIMBS]u64) [s11.width])
           (sparse_matrixes: [s11.sparse_count][s11.sparse_matrix_size]([s11.Field.LIMBS]u64))
              : s11_state =
  let constants = s11.make_constants arity_tag round_keys mds_matrix pre_sparse_matrix sparse_matrixes in
  s11.init constants

entry mbatch_hash2 [u64_count] (s: p2_state) (u64s: [u64_count]u64): ([][p2.Field.LIMBS]u64, p2_state) = p2.hash_preimages_monts s u64s
entry mbatch_hash8 [u64_count] (s: p8_state) (u64s: [u64_count]u64): ([][p8.Field.LIMBS]u64, p8_state) = p8.hash_preimages_monts s u64s
entry mbatch_hash11 [u64_count] (s: p11_state) (u64s: [u64_count]u64): ([][p11.Field.LIMBS]u64, p11_state) = p11.hash_preimages_monts s u64s

-- Strengthened
entry mbatch_hash2s [u64_count] (s: s2_state) (u64s: [u64_count]u64): ([][s2.Field.LIMBS]u64, s2_state) = s2.hash_preimages_monts s u64s
entry mbatch_hash8s [u64_count] (s: s8_state) (u64s: [u64_count]u64): ([][s8.Field.LIMBS]u64, s8_state) = s8.hash_preimages_monts s u64s
entry mbatch_hash11s [u64_count] (s: s11_state) (u64s: [u64_count]u64): ([][s11.Field.LIMBS]u64, s11_state) = s11.hash_preimages_monts s u64s

module t8_64m_hasher = t8_64m.Hasher
type t8_64m_state = t8_64m_hasher.state
entry init_t8_64m (arity_tag: ([t8_64m_hasher.Field.LIMBS]u64))
           (round_keys: [t8_64m_hasher.rk_count]([t8_64m_hasher.Field.LIMBS]u64))
           (mds_matrix: matrix ([t8_64m_hasher.Field.LIMBS]u64) [t8_64m_hasher.width])
           (pre_sparse_matrix: matrix ([t8_64m_hasher.Field.LIMBS]u64) [t8_64m_hasher.width])
           (sparse_matrixes: [t8_64m_hasher.sparse_count][t8_64m_hasher.sparse_matrix_size]([t8_64m_hasher.Field.LIMBS]u64))
              : t8_64m_state =
  let constants = t8_64m_hasher.make_constants arity_tag round_keys mds_matrix pre_sparse_matrix sparse_matrixes in
  t8_64m_hasher.init constants

entry build_tree8_64m (s: t8_64m_state) (u64s: []u64): [][t8_64m_hasher.Field.LIMBS]u64 =
  map t8_64m_hasher.Field.mont_to_u64s
      (t8_64m.build_tree s ((map t8_64m_hasher.Field.mont_from_u64s (even_chunks t8_64m_hasher.Field.LIMBS u64s))
                            :> [t8_64m.leaves]t8_64m_hasher.Field.t))


  -- Used in example program
let x8 = p8.init p8.blank_constants
entry simple8 n = tabulate n (\i -> p8.Field.mont_to_u64s
                                    (p8.hash_preimage x8 ((replicate 8 (p8.Field.mont_from_u32 (u32.i64 i)) :> [p8.arity]p8.Field.t))))
