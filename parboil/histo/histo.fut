-- Note that this only computes one iteration of the histogram
-- (corresponding to numIterations==1 in the Parboil implementation).
-- This is because the numIterations-loop is pointless repeated
-- busywork, and the Futhark compiler would just optimise it to a
-- single iteration anyway.
-- ==
-- compiled input @ data/default.in.gz
-- output @ data/default.out.gz
-- compiled input @ data/large.in.gz
-- output @ data/large.out.gz

let sat_add_u8 (x: i32) (y: i32): i32 =
  if x + y > i32.u8 u8.highest then i32.u8 u8.highest else x + y

let main [img_width] [img_height] (histo_width: i32) (histo_height: i32)
                                  (img: [img_width][img_height]i32) =
  let histo = replicate (i64.i32 (histo_height*histo_width)) 0
  let flat = img_height*img_width
  in reduce_by_index histo sat_add_u8 0
                     (flatten_to flat img |> map i64.i32)
                     (replicate flat 1)
