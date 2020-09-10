-- | Efficient statically sized vectors.
--
-- This file defines a module type, `vector`@mtype, that describes
-- single-dimensional vectors with an abstract interface.  The point
-- is that the interface has been designed such that it can be
-- implemented both with an ordinary Futhark array
-- (`any_vector`@term), but also with a tuple-representation for those
-- cases where the size of the vector is known statically
-- (`vector_1`@term and `cat_vector`@term).  This is useful for small
-- vectors, because they will then be kept in registers rather than
-- memory.  It also makes it explicit to the compiler that the
-- potential parallelism provided by the vector operations is not
-- worth exploiting.  This is generally the case when working with
-- many points in some space.
--
-- ## Examples
--
-- The `vector_1`@term module provides us with single-element vectors.
-- We can build those up with `cat_vector` to produce vectors of any
-- (static) size:
--
-- ```
-- module vector_2 = cat_vector vector_1 vector_1
-- module vector_3 = cat_vector vector_2 vector_1
-- module vector_5 = cat_vector vector_2 vector_3
-- module vector_8 = cat_vector vector_5 vector_3
--
-- let main (xs: [8]i32) =
--   xs
--   |> vector_8.from_array
--   |> vector_8.map (+1)
--   |> vector_8.to_array
-- ```
--
-- This is awkward for very long vectors, but these would not be
-- efficient anyway.

-- | A module type describing a one-dimensional vector (although a
-- length-n vector may well represent a point in an n-dimensional
-- space).  Note the absence of run-time vector concatenation - this
-- is omitted to facilitate statically sized representations, such as
-- tuples.
module type vector = {
  type vector 'a

  -- | Apply a function to every element of the vector.
  val map 'a 'b: (a -> b) -> vector a -> vector b

  -- | Apply a function to every pair of elements from two vectors.
  val map2 'a 'b 'c: (a -> b -> c) -> vector a -> vector b -> vector c

  -- | Perform a reduction of the vector with an associative operator
  -- that has the provided neutral element.
  val reduce 'a : (a -> a -> a) -> a -> vector a -> a

  -- | Construct a vector of pairs from two vectors.
  val zip 'a 'b : vector a -> vector b -> vector (a,b)

  -- | Turn a vector of arrays into an array of vectors.
  val vzip [n] 'a : vector ([n]a) -> [n](vector a)

  -- | Turn a array of vectors a vector of arrays.
  val vunzip [n] 'a : [n](vector a) -> vector ([n]a)

  -- | A vector with increasing elements, starting at 0.
  val iota : vector i64

  -- | A vector with filled with a replicated value
  val replicate 'a : a -> vector a

  -- | Retrieve the element at some position.
  val get 'a: i64 -> vector a -> a

  -- | Set the element at some position.
  val set 'a: i64 -> a -> vector a -> vector a

  -- | Perform a left-fold over the vector's elements.
  val foldl 'a 'b : (b -> a -> b) -> b -> vector a -> b

  -- | Perform a right-fold over the vector's elements.
  val foldr 'a 'b : (a -> b -> b) -> b -> vector a -> b

  -- | The length of vectors.
  val length : i64

  -- | Convert a vector to an array.
  val to_array 'a: vector a -> [length]a

  -- | Create a vector from an array.
  val from_array 'a : [length]a -> vector a
}

-- | An implementation of `vector`@mtype that uses an ordinary array
-- as the representation.  This is efficient for large vectors and the
-- implementation of `map` and `reduce` is parallel.  However, the
-- `set` operation is very slow for this implementation, as it
-- requires a full copy of the vector.  Unless you are interacting
-- with a parametric module that specifically requires an
-- implementation of `vector`@mtype, it is probably better to use
-- arrays directly than to go through this module.
--
-- When using this module, you need to instantiate it with another
-- module that indicates the dimensionality of the vectors you will be
-- producing.
module any_vector(P: { val length : i64 }) : vector = {
  let stdreplicate = replicate

  let length = P.length
  type vector 'a = [length]a
  let map = map
  let map2 = map2
  let reduce = reduce
  let zip = zip
  let vzip = transpose
  let vunzip = transpose
  let iota = iota length
  let replicate a = stdreplicate length a
  let get (i: i64) a = a[i]
  let set (i: i64) v a = copy a with [i] = v
  let foldl = foldl -- Prelude foldl.
  let foldr = foldr -- Prelude foldr.
  let to_array = id
  let from_array = id
}

-- | A module implementing statically-sized single-element vectors.
-- The implementation of `map` and `reduce` is sequential.
module vector_1 : vector = {
  type vector 'a = a

  let stdreplicate = replicate

  let map f a = f a
  let map2 f a b = f a b
  let reduce f ne a = f ne a
  let zip a b = (a, b)
  let vzip = id
  let vunzip = id
  let iota = 0i64
  let replicate a = a
  let get _ a = a
  let set _ x _ = x
  let length = 1i64
  let foldl f b x = f b x
  let foldr f b x = f x b
  let to_array a = stdreplicate length a
  let from_array as = as[0]
}

-- | Concatenation of statically sized vector modules.  This is used
-- to build up statically sized vectors of some size.  For example,
-- `cat_vector vector_1 vector_1` produces a module that defines
-- length-2 vectors.  The implementation of `map` and `reduce` is
-- sequential.
--
-- The `foldl` implementation is unrolled, so beware code explosion if
-- you use a large vector or a complex fold function.  You can always
-- use an ordinary `loop` and `get` instead.
module cat_vector (X: vector) (Y: vector): vector = {
  type vector 'a = (X.vector a, Y.vector a)

  let vzip (xs, ys) =
    zip (X.vzip xs) (Y.vzip ys)

  let vunzip v =
    let (xs, ys) = unzip v
    in (X.vunzip xs, Y.vunzip ys)

  let map f (xs, ys) = (X.map f xs, Y.map f ys)
  let map2 f (xs_a, ys_a) (xs_b, ys_b) = (X.map2 f xs_a xs_b, Y.map2 f ys_a ys_b)
  let reduce f ne (xs, ys) = X.reduce f ne xs `f` Y.reduce f ne ys
  let zip (xs_a, ys_a) (xs_b, ys_b) = (X.zip xs_a xs_b, Y.zip ys_a ys_b)

  let iota = (X.iota, Y.map (+X.length) Y.iota)
  let replicate a = (X.replicate a, Y.replicate a)
  let get i (xs, ys) = if i < X.length then X.get i xs else Y.get (i-X.length) ys
  let set i v (xs, ys) = if i < X.length then (X.set i v xs, ys)
                         else (xs, Y.set (i-X.length) v ys)
  let length = X.length + Y.length
  let foldl f b (xs, ys) = Y.foldl f (X.foldl f b xs) ys
  let foldr f b (xs, ys) = X.foldr f (Y.foldr f b ys) xs
  let to_array 't (xs, ys) = X.to_array xs ++ Y.to_array ys :> [length]t
  let from_array 't as = let xs = X.from_array (take X.length as)
                         let ys = Y.from_array (take Y.length (drop X.length as))
                         in (xs, ys)
}
