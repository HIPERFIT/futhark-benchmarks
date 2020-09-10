-- | ignore

import "vector"

module vector_2 = cat_vector vector_1 vector_1
module vector_3 = cat_vector vector_2 vector_1
module vector_5 = cat_vector vector_2 vector_3
module vector_8 = cat_vector vector_5 vector_3

-- ==
-- entry: main
-- input { [1,2,3,4,5,6,7,8] }
-- input { [1,3,5,7,9,11,13,15] }
let main (xs: [vector_8.length]i32) =
  xs
  |> vector_8.from_array
  |> vector_8.map i64.i32
  |> vector_8.map2 (+) vector_8.iota
  |> vector_8.map i32.i64
  |> vector_8.to_array

-- ==
-- entry: error
-- input { [1,2,3] }
-- error:

entry error (xs: []i32) =
  xs
  |> vector_8.from_array
  |> vector_8.to_array

-- ==
-- entry: foldl_add
-- input { [1,2,3,4,5,6,7,8] }
-- output { 36 }
entry foldl_add (xs: [vector_8.length]i32) =
  xs
  |> vector_8.from_array
  |> vector_8.foldl (+) 0

-- ==
-- entry: foldr_add
-- input { [1,2,3,4,5,6,7,8] }
-- output { 36 }
entry foldr_add (xs: [vector_8.length]i32) =
  xs
  |> vector_8.from_array
  |> vector_8.foldr (+) 0

-- ==
-- entry: foldl_sub
-- input { [1,2,3,4,5,6,7,8] }
-- output { -36 }
entry foldl_sub (xs: [vector_8.length]i32) =
  xs
  |> vector_8.from_array
  |> vector_8.foldl (-) 0

-- ==
-- entry: foldr_sub
-- input { [1,2,3,4,5,6,7,8] }
-- output { -4 }
entry foldr_sub (xs: [vector_8.length]i32) =
  xs
  |> vector_8.from_array
  |> vector_8.foldr (-) 0
