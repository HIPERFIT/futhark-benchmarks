-- An implementation of the tpacf benchmark from the parboil benchmark suite
--
-- ==
-- compiled input @ data/small.in.gz
-- output @ data/small.out
-- compiled input @ data/medium.in.gz
-- output @ data/medium.out
-- compiled input @ data/large.in.gz
-- output @ data/large.out

type vec3 = (f64, f64, f64)

let dec2rad(dec: f64): f64 = f64.pi/180.0 * dec
let rad2dec(rad: f64): f64 = 180.0f64/f64.pi * rad
let min_arcmin: f64 = 1.0
let max_arcmin: f64 = 10000.0
let bins_per_dec: f64 = 5.0
let numBins: i32 = 20

let iota32 (num: i64): [num]f64 =
  map f64.i64 (iota num)

let sumBins [numBins][numBinss] (bins: [numBinss][numBins]i32): *[numBins]i32 =
  map1 i32.sum (transpose bins)

let log10 (num: f64): f64 = f64.log(num) / f64.log(10.0)

let doCompute [numD][numBBins]
              (datapoints: [numD]vec3)
              (random: [numD]vec3)
              (numBins2: i64)
              (binb: [numBBins]f64)
             : [numBins2]i32 =
  let f k =
    let i = k / numD
    let j = k % numD
    let (xOuter, yOuter, zOuter) = datapoints[i]
    let (xInner, yInner, zInner) = random[j]
    let dot = xOuter * xInner + yOuter * yInner + zOuter * zInner
    let (min,max) =
      loop (min, max) = (0, numBins) while (min+1) < max do
        let k = (min+max) / 2
        in #[unsafe] if dot >= binb[k]
                  then (min, k)
                  else (k, max)

    let index = #[unsafe] if dot >= binb[min]
                then min
                else if dot < binb[max]
                    then max+1
                    else max

    in i64.i32 index

  let l = numD * numD
  in reduce_by_index (replicate numBins2 0) (+) 0 (tabulate l f) (replicate l 1)

let doComputeSelf [numD][numBBins]
                  (datapoints: [numD]vec3)
                  (numBins2: i64)
                  (binb: [numBBins]f64)
                 : [numBins2]i32 =
  let f k =
    let i = k / numD
    let j = k % numD
    in
    if j > i then -1
    else let (xOuter, yOuter, zOuter) = datapoints[i]
         let (xInner, yInner, zInner) = datapoints[j]
         let dot = xOuter * xInner + yOuter * yInner + zOuter * zInner
         let (min,max) = loop (min, max) = (0, numBins) while (min+1) < max do
           let k = (min+max) / 2
           in #[unsafe] if dot >= binb[k]
                     then (min, k)
                     else (k, max)

         let index = #[unsafe] if dot >= binb[min]
                     then min
                     else if dot < binb[max]
                         then max+1
                         else max

         in i64.i32 index
  let l = numD * numD
  in reduce_by_index (replicate numBins2 0) (+) 0 (tabulate l f) (replicate l 1)

let fixPoints (ra: f64) (dec: f64): vec3 =
  let rarad = dec2rad(ra)
  let decrad = dec2rad(dec)
  let cd = f64.cos(decrad)
  in (f64.cos(rarad)*cd, f64.sin(rarad)*cd, f64.sin(decrad))

let main [numD][numRs]
         (datapointsx: [numD]f64)
         (datapointsy: [numD]f64)
         (randompointsx: [numRs][numD]f64)
         (randompointsy: [numRs][numD]f64)
        : *[60]i32 =
  let numBins2 = i64.i32 numBins + 2
  let binb = map (\k ->
                    f64.cos((10.0 ** (log10(min_arcmin) + k*1.0/bins_per_dec))
                            / 60.0 * dec2rad(1.0)))
                 (iota32(i64.i32 numBins + 1))
  let datapoints = map2 fixPoints datapointsx datapointsy
  let randompoints = map2 (map2 fixPoints) randompointsx randompointsy
  let (drs, rrs) = unzip (map (\random ->
                                 (doCompute datapoints random numBins2 binb,
                                  doComputeSelf random numBins2 binb))
                               randompoints)
  let dd = doComputeSelf datapoints numBins2 binb
  let rr = sumBins rrs
  let dr = sumBins drs
  in [dd, dr, rr] |> transpose |> init |> tail |> flatten_to 60
