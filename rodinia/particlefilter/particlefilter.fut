-- | Port of the OpenMP version.  Matches original names and some
-- comments as closely as possible.  One complication is that the
-- original approach for handling random numbers is not directly
-- applicable in Futhark (due to purity).  However, the main
-- difference is that we are using proper multidimensional arrays
-- here.
--
-- The datasets we use are just noise.  With one exception, the
-- control flow of the particle filter does not depend on the concrete
-- values.  The exception is the findIndex function, and here it is
-- likely that the use of random noise actually disadvantages the
-- Futhark implementation.  The original Rodinia implementation
-- contains some crusty code for generating a simple data set, but it
-- is not parallelised or measured, so we leave it out.
--
-- The original Rodinia implemention prints a running description of
-- the 'distance' value (and a few others).  This is clearly not
-- possible in Futhark, so we instead produce an array of them.  It is
-- likely (and expected) that several of these will be infinite or
-- NaN.  This is also the case with the original Rodinia
-- implementation.
-- ==
-- compiled input @ data/128_128_10_image_10000_particles.in
-- output @ data/128_128_10_image_10000_particles.out
-- compiled input @ data/128_128_10_image_400000_particles.in
-- output @ data/128_128_10_image_400000_particles.out

import "lib/github.com/diku-dk/cpprandom/random"

module rnge = linear_congruential_engine i32 {
  let m = i32.highest
  let a: i32 = 1103515245
  let c: i32 = 12345
}

module norm_rng = normal_distribution f64 rnge

-- | Fills a radius x radius matrix representing the disk
let strelDisk (radius: i64): [][]i32 =
  let diameter = radius*2 - 1
  in tabulate_2d diameter diameter
                 (\x y -> let distance = f64.sqrt(f64.i64 (x-radius+1)**2 +
                                                  f64.i64 (y-radius+1)**2)
                          in if distance < f64.i64 radius then 1 else 0)

let getneighbors [diameter] (se: [diameter][diameter]i32) (radius: i64): [](f64,f64) =
  let center = radius - 1
  in se
     |> map2 (\x -> map2 (\y v -> (v != 0, (f64.i64 (y-center),
                                            f64.i64 (x-center))))
                         (iota diameter))
             (iota diameter)
     |> flatten
     |> filter (.0)
     |> map (.1)

-- | Finds the first element in the CDF that is greater than or equal
-- to the provided value and returns that index
let findIndex [Nparticles] (CDF: [Nparticles]f64) (value: f64): i32 =
  i32.min (i32.i64 Nparticles-1)
          (loop x = 0 while x < i32.i64 Nparticles &&
                            value < #[unsafe] CDF[x]
           do x + 1)

let particleFilter [IszX][IszY][Nfr]
                   (I: [IszX][IszY][Nfr]i32) (Nparticles: i64) =
  let seed = rnge.rng_from_seed [123] |> rnge.split_rng Nparticles
  let xe = f64.round (f64.i64 IszY / 2)
  let ye = f64.round (f64.i64 IszX / 2)
  let radius = 5i64
  let [diameter] (disk: [diameter][diameter]i32) = strelDisk radius
  let [countOnes] (objxy: [countOnes](f64, f64)) = getneighbors disk radius
  let weights = replicate Nparticles (1 / f64.i64 Nparticles)
  let arrayX = replicate Nparticles xe
  let arrayY = replicate Nparticles ye
  let seed0 = rnge.join_rng seed
  let distances = replicate Nfr 0
  let (_, _, distances, _, _) =
    -- I have no idea why this should run Nfr-1 times instead of Nfr times.
    loop (arrayX, arrayY, distances, seed, seed0) for k in 1..<Nfr do
    let (seed, x_noises) = seed
                           |> map (norm_rng.rand {stddev=0, mean=5})
                           |> unzip
    let (seed, y_noises) = seed
                           |> map (norm_rng.rand {stddev=0, mean=2})
                           |> unzip
    let arrayX = map2 (+) arrayX (map (1+) x_noises)
    let arrayY = map2 (+) arrayY (map (-2+) y_noises)

    let flikelihood x' y' =
      let ind = map
        (\(a,b) -> let indX = i64.f64 (f64.round x' + b)
                   let indY = i64.f64 (f64.round y' + a)
                   --- Weird: if either computed index is out of
                   --- bounds, then we just read element [0,0]
                   --- instead.  The original implementation did the
                   --- same.  Maybe this is a safety check that never
                   --- occurs in practice.
                   in if indX >= IszX || indY >= IszY then (0,0) else (indX, indY)
        ) objxy
      in (map2 (-) (map (\(i,j) -> (#[unsafe] I[i,j,k]-100)**2) ind)
                   (map (\(i,j) -> (#[unsafe] I[i,j,k]-228)**2) ind))
         |> map f64.i32
         |> map (/50)
         |> f64.sum
         |> (/f64.i64 countOnes)

    let likelihood = map2 flikelihood arrayX arrayY

    let weights = map2 (*) weights (map f64.exp likelihood)
    let sumWeights = f64.sum weights
    let weights = map2 (*) weights (map (/sumWeights) weights)
    let xe = f64.sum (map2 (*) arrayX weights)
    let ye = f64.sum (map2 (*) arrayY weights)
    let distance = f64.sqrt ((xe-f64.round(f64.i64 IszY/2))**2 + (ye-f64.round(f64.i64 IszX/2))**2)
    let distances[k] = distance

    let CDF = scan (+) 0 weights
    let (seed0, v) = norm_rng.rand {mean=0, stddev=1} seed0
    let u1 = (1/f64.i64 Nparticles) * v
    let u = map (f64.i64 >-> (/f64.i64 Nparticles) >-> (+u1)) (iota Nparticles)
    let (xj, yj) = u
                   |> map (\u' -> let i = findIndex CDF u'
                                  in #[unsafe] (arrayX[i], arrayY[i]))
                   |> unzip

    in (xj, yj, distances, seed, seed0)
  in distances

let main = particleFilter
