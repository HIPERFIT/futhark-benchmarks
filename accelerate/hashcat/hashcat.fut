-- A port of the Hashcat benchmark from Accelerate, originally by Mads
-- Ynddal and Troels Ynddal.
--
-- The implementation has some limitations that are similar to the
-- ones in Accelerate.  Most importantly, all candidate passwords
-- assumed to take up a single MD5 block (64 bytes), and all input is
-- pre-padded to ensure this property.  Passwords beyond 64 bytes are
-- truncated.  Accelerate does not incorporate the cost of this
-- processing into the overall runtime, but we do.
--
-- The test data set is a truncated (for Git repository space reasons)
-- version of the rockyou password dump, which tries to find the tenth
-- password out of a dictionary containing two million candidates.
-- ==
-- compiled input @ rockyou.dataset
-- output { 10 }

type md5 = (u32, u32, u32, u32)

let us32 (x: i32) = u32.i32 x

let rs: [64]u32 =
  map us32
  [ 7, 12, 17, 22,  7, 12, 17, 22,  7, 12, 17, 22,  7, 12, 17, 22,
    5,  9, 14, 20,  5,  9, 14, 20,  5,  9, 14, 20,  5,  9, 14, 20,
    4, 11, 16, 23,  4, 11, 16, 23,  4, 11, 16, 23,  4, 11, 16, 23,
    6, 10, 15, 21,  6, 10, 15, 21,  6, 10, 15, 21,  6, 10, 15, 21 ]

let ks: [64]u32 =
  [ 0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee ,
    0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501 ,
    0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be ,
    0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821 ,
    0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa ,
    0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8 ,
    0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed ,
    0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a ,
    0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c ,
    0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70 ,
    0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05 ,
    0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665 ,
    0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039 ,
    0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1 ,
    0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1 ,
    0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391 ]

let rotate_left(x: u32, c: u32): u32 = (x << c) | (x >> (32u32 - c))

let bytes(x: u32): [4]u8 = [u8.u32(x),
                            u8.u32(x >> 8u32),
                            u8.u32(x >> 16u32),
                            u8.u32(x >> 24u32)]

let unbytes(bs: [4]u8): u32 =
  u32.u8(bs[0]) +
  (u32.u8(bs[1]) << 8u32) +
  (u32.u8(bs[2]) << 16u32) +
  (u32.u8(bs[3]) << 24u32)

-- Process 512 bits of the input.
let md5_chunk ((a0,b0,c0,d0): md5) (block: [16]u32): md5 =
  loop (a,b,c,d) = (a0,b0,c0,d0) for i < 64 do
    let (f,g) =
      if      i < 16 then ((b & c) | (!b & d),
                           u32.i32 i)
      else if i < 32 then ((d & b) | (!d & c),
                           (5u32*u32.i32 i + 1u32) % 16u32)
      else if i < 48 then (b ^ c ^ d,
                           (3u32*u32.i32 i + 5u32) % 16u32)
      else                (c ^ (b | !d),
                           (7u32*u32.i32 i)        % 16u32)
    in (d, b + rotate_left(a + f + ks[i] + #[unsafe] block[i32.u32 g], rs[i]), b, c)

let initial_chunk = (0x67452301_u32,
                     0xefcdab89_u32,
                     0x98badcfe_u32,
                     0x10325476_u32)

let md5s [n] (bs: [n][16]u32): [n]md5 =
  map (md5_chunk initial_chunk) bs


let first_true [n] (bools : [n]bool): i32 =
  (reduce_comm (\p (i,x) -> if x then (i,x) else p)
               (-1,false)
               (zip (map i32.i64 (iota n)) bools)).0

-- Input preprocessing.

type dict [n] = [n][16]u32

let mk_block (bs: []u8) ((i,k): (i64,i64)): [16]u32 =
  #[unsafe]
  let k = i64.min 64 k -- Truncate past first block.
  let one_bit = [0x80u8, 0u8, 0u8, 0u8]
  let block = replicate 64 0u8
  let block[0:k] = bs[i:i+k]
  let block[k:k+4] = one_bit
  let block[64-8:64-4] = bytes (u32.i64(k*8))
  in map unbytes (unflatten 16 4 block)

let md5_blocks [n][k] (bs: [k]u8) (offsets: [n]i64): [n][16]u32 =
  let lengths = map (\(i, j) -> if i > j then k-i else j-i)
                    (zip offsets (rotate 1 offsets))
  in map (mk_block bs) (zip offsets lengths)

entry mk_dict [n] (bs: []i8) (offsets: [n]i64): dict[n] =
  md5_blocks (map u8.i8 bs) offsets

-- The entry point.

entry main [n] (a: u32) (b: u32) (c: u32) (d: u32) (bs: []i8) (offsets: [n]i32): i32 =
  first_true (map (==(a,b,c,d)) (md5s (mk_dict bs (map i64.i32 offsets))))
