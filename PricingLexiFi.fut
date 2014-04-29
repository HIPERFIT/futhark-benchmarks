fun int grayCode(int x) = (x >> 1) ^ x

////////////////////////////////////////
/// Sobol Generator
////////////////////////////////////////
fun int index([int] arr, int i) = arr[i]

fun bool testBit(int n, int ind) =
    let t = (1 << ind) in (n & t) = t

/////////////////////////////////////////////////////////////////
//// DIFFICULT VERSION: filter is hoisted outside map:
////         less computation but variable array size!
/////////////////////////////////////////////////////////////////
//fun int xorInds([int] indices, [int] dir_vs ) =
//    reduce( op ^, 0, map( index(dir_vs), indices ) )
//
//fun [int] sobolIndI ( int bits_num, [[int]] dir_vs, int n ) =
//    let bits    = iota   ( bits_num ) in
//    let indices = filter ( testBit(grayCode(n)), bits )
//    in map( xorInds(indices), dir_vs )


/////////////////////////////////////////////////////////////////
//// EASY VERSION: filter is redundantly computed inside map:
////    more computation but a redofilt pattern, i.e., array
////    not instantiated!
/////////////////////////////////////////////////////////////////
fun int xorInds(int bits_num, int n, [int] dir_vs ) =
    let bits    = iota   ( bits_num )                   in
    let indices = filter ( testBit(grayCode(n)), bits ) in
    reduce( op ^, 0, map( index(dir_vs), indices ) )

fun [int] sobolIndI ( int bits_num, [[int]] dir_vs, int n ) =
    map( xorInds(bits_num, n), dir_vs )
////////////////////////////////////////////////////////////////

fun [real] sobolIndR( int bits_num, [[int]] dir_vs, int n ) =
    let divisor = 2.0 pow toReal (bits_num)        in
    let arri    = sobolIndI( bits_num, dir_vs, n ) in
        map( fn real (int x) => toReal(x) / divisor, arri )

////////////////////////////////////////
/// Inverse Gaussian
////////////////////////////////////////

//tmp_rat_evalL :: SpecReal -> [SpecReal] -> SpecReal
fun real polyAppr(      real x,
                        real a0, real a1, real a2, real a3,
                        real a4, real a5, real a6, real a7,
                        real b0, real b1, real b2, real b3,
                        real b4, real b5, real b6, real b7
                    ) =
        (x*(x*(x*(x*(x*(x*(x*a7+a6)+a5)+a4)+a3)+a2)+a1)+a0) /
        (x*(x*(x*(x*(x*(x*(x*b7+b6)+b5)+b4)+b3)+b2)+b1)+b0)

fun real smallcase(real q) =
        q * polyAppr( 0.180625 - q * q,

                      3.387132872796366608,
                      133.14166789178437745,
                      1971.5909503065514427,
                      13731.693765509461125,
                      45921.953931549871457,
                      67265.770927008700853,
                      33430.575583588128105,
                      2509.0809287301226727,

                      1.0,
                      42.313330701600911252,
                      687.1870074920579083,
                      5394.1960214247511077,
                      21213.794301586595867,
                      39307.89580009271061,
                      28729.085735721942674,
                      5226.495278852854561
                    )

fun real intermediate(real r) =
        polyAppr( r - 1.6,

                  1.42343711074968357734,
                  4.6303378461565452959,
                  5.7694972214606914055,
                  3.64784832476320460504,
                  1.27045825245236838258,
                  0.24178072517745061177,
                  0.0227238449892691845833,
                  7.7454501427834140764e-4,

                  1.0,
                  2.05319162663775882187,
                  1.6763848301838038494,
                  0.68976733498510000455,
                  0.14810397642748007459,
                  0.0151986665636164571966,
                  5.475938084995344946e-4,
                  1.05075007164441684324e-9
                )

fun real tail(real r) =
        polyAppr( r - 5.0,

                  6.6579046435011037772,
                  5.4637849111641143699,
                  1.7848265399172913358,
                  0.29656057182850489123,
                  0.026532189526576123093,
                  0.0012426609473880784386,
                  2.71155556874348757815e-5,
                  2.01033439929228813265e-7,

                  1.0,
                  0.59983220655588793769,
                  0.13692988092273580531,
                  0.0148753612908506148525,
                  7.868691311456132591e-4,
                  1.8463183175100546818e-5,
                  1.4215117583164458887e-7,
                  2.04426310338993978564e-5
                )

fun real ugaussianEl( real p ) =
    let dp = p - 0.5
    in  //if  ( fabs(dp) <= 0.425 )
        if ( ( (dp < 0.0) && (0.0 - dp <= 0.425) ) || ( (0.0 <= dp) && (dp <= 0.425) ) )
        then smallcase(dp)
        else let pp = if(dp < 0.0) then dp + 0.5 else 0.5 - dp      in
             let r  = sqrt( - log(pp) )                             in
             let x = if(r <= 5.0) then intermediate(r) else tail(r) in
                if(dp < 0.0) then 0.0 - x else x

// Transforms a uniform distribution [0,1) into a gaussian distribution (-inf, +inf)
fun [real] ugaussian([real] ps) = map(ugaussianEl, ps)


/////////////////////////////////
/// Brownian Bridge
/////////////////////////////////

fun [real] brownianBridgeDates (
                  int    num_dates,
                [[int ]] bb_inds,       // [3,  num_dates]
                [[real]] bb_data,       // [3,  num_dates]
                 [real]  gauss          // [num_dates]
            ) =
    let bi = bb_inds[0] in
    let li = bb_inds[1] in
    let ri = bb_inds[2] in
    let sd = bb_data[0] in
    let lw = bb_data[1] in
    let rw = bb_data[2] in

//  let gauss[ bi[0]-1 ] = sd[0] * gauss[0]  in
    let bbrow = replicate(num_dates, 0.0)   in
    let bbrow[ bi[0]-1 ] = sd[0] * gauss[0] in

    loop (bbrow) =
        for i < num_dates-1 do  // use i+1 since i in 1 .. num_dates-1
            let j  = li[i+1] - 1 in
            let k  = ri[i+1] - 1 in
            let l  = bi[i+1] - 1 in

            let wk = bbrow [k  ] in
            let zi = gauss [i+1] in
            let tmp= rw[i+1] * wk + sd[i+1] * zi in

            let bbrow[ l ] = if( (j + 1) = 0)   // if(j=-1)
                             then tmp
                             else tmp + lw[i+1] * bbrow[j]
            in  bbrow

        // This can be written as map-reduce, but it
        //   needs delayed arrays to be mapped nicely!
    in loop (bbrow) =
        for ii < num_dates-1 do
            let i = num_dates - (ii+1) in
            let bbrow[i] = bbrow[i] - bbrow[i-1]
            in  bbrow
       in bbrow

// [num_dates,num_paths]
fun [[real]] brownianBridge (
                  int    num_paths,
                  int    num_dates,
                [[int ]] bb_inds,       // [3,  num_dates]
                [[real]] bb_data,       // [3,  num_dates]
                 [real]  gaussian_arr   // [num_paths * num_dates]
            ) =
    let gauss2d  = reshape((num_dates,num_paths), gaussian_arr) in
    let gauss2dT = transpose(gauss2d) in
      transpose( map( brownianBridgeDates(num_dates, bb_inds, bb_data), gauss2dT ) )


/////////////////////////////////
/// Black-Scholes
/////////////////////////////////

fun real zwop(real a, real b, int j) = a * b

fun [real] take(int n, [real] a) = let {first, rest} = split(n, a) in first

fun [real] fftmp(int num_paths, [[real]] md_c, [real] zi) =
    map( fn real (int j) =>
            let x = map(zwop, zip(take(j+1,zi), take(j+1,md_c[j]), iota(j+1)))
            in  reduce(op +, 0.0, x),
         iota(num_paths)
       )

fun [[real]] correlateDeltas(int num_paths, [[real]] md_c, [[real]] zds) =
    map( fftmp(num_paths, md_c), zds )

fun [real] combineVs([real] n_row, [real] vol_row, [real] dr_row) =
    map( op +, zip(dr_row, map( op *, zip(n_row, vol_row ) )))

 fun [[real]] mkPrices ([real] md_starts, [[real]] md_vols, [[real]] md_drifts, [[real]] noises) =
    let e_rows = map( fn [real] ([real] x) => map(exp, x),
                      map(combineVs, zip(noises, md_vols, md_drifts))
                    )
        // If we use the scan(op *, e, [x1,..,xn]) = [e*x1, e*x1*x2,...,e*x1*x2*...xn]
    in  scan( fn [real] ([real] x, [real] y) => map(op *, zip(x, y)), md_starts, e_rows )

        // If we use the scan(op *, e, [x1,..,xn]) = [e, e*x1, ..., e*..*xnm1]
    //  let tmp = scan( fn [real] ([real] x, [real] y) => zipWith(op *, x, y), md_starts, e_rows )
    //  in  zipWith(zipWith(op *), tmp, e_rows)
        // If we use the scan(op *, e, [x1,...,xn]) = [e, e*x1, .., e*x1*..*xn], i.e., Haskell's scan then:
    // tail( scan( fn [real] ([real] x, [real] y) => zipWith(op *, x, y), md_starts, e_rows ) )



//[num_dates, num_paths]
fun [[real]] blackScholes(
                int      num_paths,
                [[real]] md_c,         //[num_paths, num_paths]
                [[real]] md_vols,      //[num_paths, num_dates]
                [[real]] md_drifts,    //[num_paths, num_dates]
                 [real]  md_starts,    //[num_paths]
                [[real]] bb_arr        //[num_dates,num_paths]
           ) =
    let noises = correlateDeltas(num_paths, md_c, bb_arr)
    in  mkPrices(md_starts, md_vols, md_drifts, noises)

////////////////////////////////////////
// MAIN
////////////////////////////////////////

fun real main(int num_mc_it, int num_dates, int num_und, int num_bits, [[int]] dir_vs,
             [[real]] md_c, [[real]] md_vols, [[real]] md_drifts, [real] md_st, [real] md_dv, [real] md_disc,
             [[int]] bb_inds, [[real]] bb_data) =
    let sobol_mat = map ( sobolIndR(num_bits, dir_vs), map(fn int (int x) => x + 1, iota(num_mc_it)) ) in
    //let x = write(sobol_mat) in
    let gauss_mat = map ( ugaussian, sobol_mat )                                       in
    //let x = write(gauss_mat) in
    let bb_mat    = map ( brownianBridge( num_und, num_dates, bb_inds, bb_data ), gauss_mat )    in
    let bs_mat    = map ( blackScholes( num_und, md_c, md_vols, md_drifts, md_st ), bb_mat ) in
    //let x = write(bs_mat) in

    let payoffs   = map ( payoff2(md_disc), bs_mat ) in
    let payoff    = reduce ( op +, 0.0, payoffs )       in
    payoff / toReal(num_mc_it)

////////////////////////////////////////
// PAYOFF FUNCTION
////////////////////////////////////////

fun real payoff2 ([real] md_disc, [[real]] xss) =
// invariant: length(xss) == 5, i.e., 5x3 matrix
    let divs    = [ 1.0/3758.05, 1.0/11840.0, 1.0/1200.0 ]             in
    let xss_div = map( fn [real] ([real] xs) => map(op *, zip(xs, divs)), xss     ) in
    let mins    = map( MIN, xss_div )
    in  if( 1.0 <= mins[0] ) then trajInner(1150.0, 0, md_disc)
        else if( 1.0 <= mins[1] ) then trajInner(1300.0, 1, md_disc)
             else if( 1.0 <= mins[2] ) then trajInner(1450.0, 2, md_disc)
                  else if( 1.0 <= mins[3] ) then trajInner(1600.0, 3, md_disc)
                       else if( 1.0 <= mins[4] ) then trajInner(1750.0, 4, md_disc)
                            else if( 0.75 < mins[4] ) then trajInner(1000.0, 4, md_disc)
                                 else trajInner(1000.0 * mins[4], 4, md_disc)

fun real MIN([real] arr) = reduce( fn real (real x, real y) => if(x<y) then x else y, arr[0], arr )

fun real trajInner(real amount, int ind, [real] disc) = amount * disc[ind]


////////////////////////////////////////
// INPUT DATA
////////////////////////////////////////


// Result is:   (   #Monte-Carlo Iters, #dates, #underlyings,
//                  integer bit-len for Sobol, direction vectors
//              )
fun {int, int, int, int, [[int]]} getHLdata() =
        {
            1000,           // number of Monte-Carlo Iterations
            5,              // number of dates
            3,              // number of underlyings
            30,             // integer bit-length representation for Sobol
            getDirVs()      // direction vectors
        }


// the market/model data
//     md_c     md_vols    md_drifts  md_starts md_detval md_disc
fun {[[real]], [[real]], [[real]], [real], [real], [real]} getModelData() =
        {

            [   // md_c                 [num_paths, num_paths]
                [ 1.0, 0.6, 0.8                 ],
                [ 0.6, 0.8, 0.15                ],
                [ 0.8, 0.15, 0.5809475019311124 ]
            ],
            [   // md_vols (volatility) [num_paths, num_dates]
                [ 0.19, 0.19, 0.15 ],
                    [ 0.19, 0.19, 0.15 ],
                    [ 0.19, 0.19, 0.15 ],
                    [ 0.19, 0.19, 0.15 ],
                    [ 0.19, 0.19, 0.15 ]
            ],
            [   // md_drifts            [num_paths, num_dates]
                [ -0.0283491736871803,  0.0178771081725381, 0.0043096808044729 ],
                [ -0.0183841413744211, -0.0044530897672834, 0.0024263805987983 ],
                [ -0.0172686581005089,  0.0125638544546015, 0.0094452810918001 ],
                    [ -0.0144179417871814,  0.0157411263968213, 0.0125315353728014 ],
                    [ -0.0121497422218761,  0.0182904634062437, 0.0151125070556484 ]
            ],
            // md_starts                [pc->num_paths]
            [ 3758.0500000000001819, 11840.0, 1200.0 ],

            // model deterministic values [nb_deterministic_pricers]
            [ 0.99976705777418484188956426805817 ],

            [ // model discounts        [pc->nb_cash_flows]
              0.9797862861805930, 0.9505748482484491,
              0.9214621679912968, 0.8906693055891434, 0.8588567633110704
            ]
        }

// Brownian-Bridge metadata
fun {[[int]], [[real]]} getBBdata() =
        {   // indirect accessing
            [   [ 5, 2, 1, 3, 4 ],      // bi
                [ 0, 0, 0, 2, 3 ],      // li
                [ 0, 5, 2, 5, 5 ]       // ri
            ],
            // real data
            [   [ 2.2372928847280580, 1.0960951589853829, // bb_sd, i.e., standard deviation
                  0.7075902730592357, 0.8166828043492210,
                  0.7075902730592357
                ],
                [ 0.0000000000000000, 0.5998905309250137, // bb_lw
                  0.4993160054719562, 0.6669708029197080,
                  0.5006839945280438
                ],
                [ 0.0000000000000000, 0.4001094690749863, // bb_rw
                  0.5006839945280438, 0.3330291970802919,
                  0.4993160054719562
                ]
            ]
        }


// in principle the direction vectors should be computed
//    or read from a file, but here we inline them
fun [[int]] getDirVs() =
                [
                    [
                        536870912,268435456,134217728,67108864,33554432,16777216,8388608,4194304,
                        2097152,1048576,524288,262144,131072,65536,32768,16384,8192,
                        4096,2048,1024,512,256,128,64,32,16,8,4,2,1
                    ],
                    [
                        536870912,805306368,671088640,1006632960,570425344,855638016,713031680,1069547520,
                        538968064,808452096,673710080,1010565120,572653568,858980352,715816960,1073725440,
                        536879104,805318656,671098880,1006648320,570434048,855651072,713042560,1069563840,
                        538976288,808464432,673720360,1010580540,572662306,858993459
                    ],
                    [
                        536870912,805306368,402653184,603979776,973078528,385875968,595591168,826277888,
                        438304768,657457152,999817216,358875136,538574848,807862272,406552576,605372416,
                        975183872,389033984,597170176,828646400,437926400,656873216,1002152832,357921088,
                        536885792,805312304,402662296,603992420,973085210,385885991
                    ],
                    [
                        536870912,805306368,939524096,335544320,234881024,721420288,411041792,616562688,
                        920649728,1062207488,381157376,258736128,771883008,453181440,545488896,817971200,
                        954261504,340963328,238651392,732843008,417426944,609285376,909831040,1068349120,
                        383778848,256901168,783810616,460062740,537001998,805503019
                    ],
                    [
                        536870912,805306368,402653184,1006632960,167772160,285212672,713031680,566231040,
                        853540864,489684992,952631296,208928768,316801024,758317056,550076416,813154304,
                        417505280,1009913856,172697600,297131008,704744960,553894656,847291520,499194688,
                        954376224,204607536,306915352,766893116,536972810,805552913
                    ],
                    [
                        536870912,805306368,402653184,469762048,301989888,721420288,92274688,264241152,
                        941621248,741343232,169345024,924581888,395444224,619380736,1034256384,603963392,
                        838868992,452997120,494934016,331357184,706744832,120597248,261621120,953946048,
                        800208928,148581424,935168536,350484252,630339474,1072370923
                    ],
                    [
                        536870912,805306368,134217728,1006632960,503316480,754974720,629145600,440401920,
                        94371840,711983104,229113856,374079488,330694656,996212736,907247616,557531136,
                        867573760,190918656,1041467392,490437632,766918144,643898624,462663040,125527616,
                        672545696,202454896,373006376,288845836,1000351766,930090001
                    ],
                    [
                        536870912,268435456,402653184,872415232,838860800,956301312,612368384,717225984,
                        211812352,386924544,302514176,688128000,1015414784,516751360,1051492352,773734400,
                        914432000,63877120,807741440,165200896,748683776,118489344,168296832,486802240,
                        243663648,667747216,439124552,81674924,975249610,350138737
                    ],
                    [
                        536870912,268435456,671088640,469762048,973078528,1023410176,713031680,339738624,
                        912261120,797966336,176685056,71565312,510263296,865533952,814120960,961232896,
                        887136256,668078080,116070400,382772224,1047134720,597098752,411468416,625689024,
                        249602976,449975248,745216680,43033924,134873446,201786361
                    ],
                    [
                        536870912,268435456,402653184,67108864,704643072,385875968,696254464,205520896,
                        920649728,946864128,359137280,859045888,302907392,50659328,462192640,524599296,
                        895541248,590794752,168810496,118033408,831447552,138662144,485185920,796511296,
                        1021313184,1064304752,619184920,997458052,250479054,745865975
                    ],
                    [
                        536870912,268435456,939524096,1006632960,838860800,889192448,645922816,46137344,
                        476053504,584056832,210239488,465829888,820903936,689897472,73695232,249118720,
                        110075904,315338752,610637824,517665792,1049494016,785318144,376210304,735921088,
                        402760480,738505552,168368744,151499820,344957894,936096557
                    ],
                    [
                        536870912,805306368,939524096,1006632960,503316480,922746880,41943040,423624704,
                        228589568,651165696,195559424,500957184,791019520,261292032,1040285696,118407168,
                        982065152,625250304,329533440,298984448,153690624,76845824,579619712,692987840,
                        900670432,450334832,363187112,719119956,765461306,382730781
                    ],
                    [
                        536870912,805306368,402653184,603979776,838860800,117440512,478150656,658505728,
                        752877568,1060110336,141033472,209453056,244187136,272957440,678068224,1014546432,
                        377724928,876875776,443160576,998185984,168665600,318837504,914397568,71818816,
                        40763680,527762288,939688008,335855668,705536494,587273091
                    ],
                    [
                        536870912,268435456,671088640,738197504,637534208,150994944,813694976,943718400,
                        77594624,179306496,798490624,967049216,134348800,1006698496,235044864,620937216,
                        377643008,826314752,874711040,854819840,725109248,856992512,664336768,94804544,
                        100663328,419430416,411041832,339738668,580911142,61865993
                    ],
                    [
                        536870912,805306368,939524096,603979776,100663296,452984832,998244352,188743680,
                        866123776,389021696,287834112,172228608,824836096,977731584,153714688,507854848,
                        254402560,88403968,883578880,235160576,118055424,422917888,371224704,326210368,
                        654926368,691353392,773877944,930190180,554263078,842348331
                    ]
                ]
