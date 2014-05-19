

////////////////////////////////////////////////////////////////
///// UTILITY FUNCTIONS
////////////////////////////////////////////////////////////////

fun real pi      () = 3.1415926535897932384626433832795
fun real r       () = 0.03
fun real infinity() = 1.0e49
fun real epsilon () = 1.0e-5

fun int  itMax   () = 10000


fun real MIN(real a, real b) = if(a < b) then a else b
fun real MAX(real a, real b) = if(a < b) then b else a

fun int MINI(int a, int b) = if(a < b) then a else b
fun int MAXI(int a, int b) = if(a < b) then b else a

fun real abs(real a) = if(a < 0.0) then -a else a

////////////////////////////////////////////////////////////////
///// MATH MODULE
////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------
// Cumulative Distribution Function for a standard normal distribution

fun real uGaussian_P(real x) =
    let u = x / sqrt(2.0) in
    let e = if (u < 0.0) then -erf(-u)
                         else  erf( u)
    in 0.5 * (1.0 + e)

//-------------------------------------------------------------------------
// polynomial expansion of the erf() function, with error<=1.5e-7
//   formula 7.1.26 (page 300), Handbook of Mathematical Functions, Abramowitz and Stegun
//   http://people.math.sfu.ca/~cbm/aands/frameindex.htm

fun real erf(real x) =
    let p = 0.3275911 in
    let {a1,a2,a3,a4,a5} =
        {0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429} in
    let t  = 1.0/(1.0+p*x)  in
    let t2 = t  * t         in
    let t3 = t  * t2        in
    let t4 = t2 * t2        in
    let t5 = t2 * t3        in
         1.0 - (a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t5) * exp(-(x*x))

//-------------------------------------------------------------------------
// iteration_max = 10000 (hardcoded)


////////////////////////////////////////////////////////////////////////
// if fid = 33 then: to_solve(real x) = (x+3.0)*(x-1.0)*(x-1.0)
// otherwise follows the real implementation
////////////////////////////////////////////////////////////////////////

fun real to_solve(int fid, [{real,real}] scalesbbi, real yhat) =
    if(fid == 33) then (yhat+3.0)*(yhat-1.0)*(yhat-1.0)
    else
        let tmps = map( fn real ( {real,real} scalesbbi ) =>
                            let {scales, bbi} = scalesbbi in
                                scales * exp(-bbi*yhat)
                        , scalesbbi
                      )
            in reduce(op +, 0.0, tmps) - 1.0

/////////////////////////////////////////////////////////
//// the function-parameter to rootFinding_Brent
/////////////////////////////////////////////////////////


fun {real,int,real}
rootFinding_Brent(int fid, [{real,real}] scalesbbi, real lb, real ub, real tol, int iter_max) =
    let tol      = if(tol     <= 0.0) then 1.0e-9  else tol    in
    let iter_max = if(iter_max <= 0 ) then 10000 else iter_max in
    let {a,b}    = {lb,ub}                                     in

    // `to_solve' is curried so it will need extra arguments!!!
    // said fid refers to the version of `to_solve', e.g., testing or real implem.
    // Was previous (eval(a),eval(b)), now:
    let {fa,fb}  = { to_solve(fid, scalesbbi, a), to_solve(fid, scalesbbi, b) } in

    if(0.0 <= fa*fb)
    then
        if(0.0 <= a) then { 0.0, 0,  infinity() }  // root not bracketed above
                     else { 0.0, 0, -infinity() }  // root not bracketed below
    else
    let {fa, fb} = if( abs(fa) < abs(fb) ) then {fb, fa}
                                           else {fa, fb} in
    let {c,fc}    = {a, fa} in
    let mflag     = True    in
    let it        = 0       in
    let d         = 0.0     in

    loop ({a,b,c,d,fa,fb,fc,mflag,it}) =
        for i < iter_max do

            if( fb==0.0 || abs(b-a)<tol )
            then {a,b,c,d,fa,fb,fc,mflag,it}
            else

                // the rest of the code implements the else branch!

                let s = if( fa==fc || fb == fc )
                        then    b-fb*(b-a)/(fb-fa)

                        else let s1 = (a*fb*fc)/( (fa-fb)*(fa-fc) ) in
                             let s2 = (b*fa*fc)/( (fb-fa)*(fb-fc) ) in
                             let s3 = (c*fa*fb)/( (fc-fa)*(fc-fb) ) in
                                s1 + s2 + s3
                                                                    in

                let {mflag, s} = if ( ( not ((3.0*a+b)/4.0 <= s && s <= b)    ) ||
                                      (     mflag && abs(b-c)/2.0 <= abs(s-b) ) ||
                                      ( not mflag && abs(c-d)/2.0 <= abs(s-b) ) ||
                                      (     mflag && abs(b-c)     <= abs(tol) ) ||
                                      ( not mflag && abs(c-d)     <= abs(tol) )
                                    )
                                 then {True,  (a+b)/2.0}
                                 else {False, s        }
                                                                    in
                // Was previous `eval(s)', Now:
                let fs = to_solve(fid, scalesbbi, s)                in

                // d is assigned for the first time here:
                // it's not used above because mflag is set
                let d = c in let {c,fc} = {b, fb}                   in
                let {a,b,fa,fb} = if( fa*fs < 0.0 )
                                  then {a,s,fa,fs}
                                  else {s,b,fs,fb}                  in

                let {a,b,fa,fb} = if( abs(fa) < abs(fb) )
                                  then {b,a,fb,fa} // swap args
                                  else {a,b,fa,fb}                  in

                // reporting non-convergence!
                let dummy =
                    if(i == iter_max-1)
                    then let w = trace("# ERROR: Brent method not converged, error: ") in
                         let w = trace(fb) in 0
                    else 0

                in {a,b,c,d,fa,fb,fc,mflag,i}

    // Finally, the result of function rootFinding_Brent is:
    in { b, it, fb }




//-------------------------------------------------------------------------
// Gaussian Quadrature with Hermite linear expansion: cmulative distribution function of Normal distribution

fun {[real],[real]} gauss_hermite_coefficients() =
    {
        // coefficients
        [
            0.0, 0.6568095668820999044613, -0.6568095668820997934390, -1.3265570844949334805563,
            1.3265570844949330364670,  2.0259480158257567872226, -2.0259480158257558990442,
            -2.7832900997816496513337,  2.7832900997816474308877,  3.6684708465595856630159, -3.6684708465595838866591
        ],
        // weights
        [
            0.6547592869145917315876, 0.6609604194409607336169, 0.6609604194409606225946, 0.6812118810666693002887,
            0.6812118810666689672217, 0.7219536247283847574252, 0.7219536247283852015144, 0.8025168688510405656800,
            0.8025168688510396775015, 1.0065267861723647957461, 1.0065267861723774522886
        ]
    }
//=========================================================================
fun bool equal(real x1, real x2) =
    abs(x1-x2) <= 1.0e-8

//eval=lambda x: (x+3)*(x-1)**2
//fun real eval(real x) = (x+3.0)*(x-1.0)*(x-1.0)


fun int main_test_math() =
    // Rootfinder.brent (-4.) (4./.3.) (fun x -> (x+.3.)*.(x-.1.)**2.) 1e-4 == -3
    let tmp = trace("# Brent test: ") in
    let {root, it, error} = rootFinding_Brent(33, [{0.0,0.0}], -4.0, 4.0/3.0, 0.0, 0) in
    let tmp =   if( equal(root, -3.0) ) then trace(" success!")
                                        else trace(" fails!") in

    // erf 0. == 0. ;; 100. *. erf (1./.sqrt 2.)
    let tmp = trace("Erf test: ") in
    let tmp =   if( equal(erf(0.0), 0.0) && equal( toReal(trunc(100.0*erf(1.0/sqrt(2.0)))), 68.0 ) )
                then trace(" success!")
                else trace(" fails!") in

    // ugaussian_P 0. ;; ugaussian_P 1. +. ugaussian_P (-1.)
    let tmp = trace("Gaussian test: ") in
    let tmp =   if( equal(uGaussian_P(0.0),0.5) && equal(uGaussian_P(-1.0)+uGaussian_P(1.0),1.0) )
                then trace(" success!")
                else trace(" fails!") in
        33


  //let tmp = trace(equal(brent( a=-4, b=4/3, eval=lambda x: (x+3)*(x-1)**2 )[0], -3)
  // erf 0. == 0. ;; 100. *. erf (1./.sqrt 2.)
  //print "# Erf test:", equal(erf(0),0) and equal(int(100*erf(1/N.sqrt(2))),68)
  // ugaussian_P 0. ;; ugaussian_P 1. +. ugaussian_P (-1.)
  //print "# Gaussian test: ", equal(uGaussian_P(0),0.5) and equal(uGaussian_P(-1)+uGaussian_P(1),1)



////////////////////////////////////////////////////////////////
///// DATA MODULE
////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// FILTERED VERSION:
// RESULT is of shape:
// ((matrity_in_year*swap_frequency*swap_term_in_year) * volatility)
//fun [( (real,real,real)*real )] getSwaptionQuotes() =
//    [
//        ( (1.0,  6.0, 1.0 ), 1.052   ),
//        ( (10.0, 6.0, 1.0 ), 0.2496  ),
//        ( (30.0, 6.0, 1.0 ), 0.23465 ),
//        ( (1.0,  6.0, 10.0), 0.41765 ),
//        ( (10.0, 6.0, 10.0), 0.2775  ),
//        ( (30.0, 6.0, 10.0), 0.22    ),
//        ( (1.0,  6.0, 25.0), 0.38115 ),
//        ( (10.0, 6.0, 25.0), 0.26055 ),
//        ( (30.0, 6.0, 25.0), 0.1686  ),
//        ( (30.0, 6.0, 30.0), 0.16355 )
//    ]

////////////////////////////////////////////////////////////////
///// G2PP Module
////////////////////////////////////////////////////////////////

fun real zc(real t) = exp(-r() * date_act_365(t, today()))

fun {real,real,real} accumSched({real,real,real} xx, {real,real,real} yy ) =
    let {x,d1,d2} = xx in let {y,tf,tp} = yy in
        { x + zc(tp) * date_act_365(tp, tf), MIN(d1,tf), MAX(d2,tp) }

///////////////////////////////////////////
//// the first param `swaption' is a triple of reals,
////    i.e., the swaption maturity date, frequency of payment and
////          swapt-term in years (how often to vary the condition of the contract)
////  the result is also a triple:
////    the maturity time stamp,
////    the range of time stamps of for each swap term : [(real,real)]
////    the strike price
///////////////////////////////////////////
fun {real,[{real,real}],real}
extended_swaption_of_swaption({real,real,real} swaption)  =  // swaption = (sw_mat, freq, sw_ty)
    let {sw_mat, freq, sw_ty} = swaption          in
    let maturity   = add_years( today(), sw_mat ) in
    let nschedule  = trunc(12.0 * sw_ty / freq)   in
    let swap_sched = map( fn {real,real}(int i) =>
                                let a = add_months(maturity, (toReal(i)*freq)) in
                                    {a, add_months(a,freq)}
                          , iota(nschedule)
                        ) in

    let {swap_sched1, swap_sched2} = unzip(swap_sched)                              in
    let swap_sched_new = zip(replicate(nschedule, 0.0), swap_sched1, swap_sched2)   in

    let {lvl,t0,tn}= reduce(accumSched, {0.0, max_date(), min_date()}, swap_sched_new) in
    let strike     = (zc(t0) - zc(tn)) / lvl

        in {maturity, swap_sched, strike}


fun real b_fun(real z, real tau) = (1.0-exp(-z*tau))/z

fun real t_fun(real sigma, real x, real tau) =
    let expxtau  = exp(-x*tau)     in
    let exp2xtau = expxtau*expxtau in
        sigma*sigma/(x*x)*(tau+2.0/x*expxtau-1.0/(2.0*x)*exp2xtau-3.0/(2.0*x))


///////////////////////////////////////////////////////////////////
// the first parameter `genome' is the five-genes genome used in
//     the genetic algorithms that models the interest rate
// the second parameter is the time
// the result is V in Brigo and Mercurio's book page 148.
//     \var(\int_t^T [x(u)+y(u)]du)
///////////////////////////////////////////////////////////////////
fun {real,real,real} bigv({real,real,real,real,real} genome, real tau) =
    let {g_a, g_b, g_rho, g_nu, g_sigma} = genome in

    // sanity check; this check should be hoisted higher up
    let g_sigma = if(g_sigma == 0.0) then 1.0e-10 else g_sigma in

    let ba = b_fun(g_a,        tau) in
    let bb = b_fun(g_b,        tau) in
    let t1 = t_fun(g_sigma,g_a,tau) in
    let t2 = t_fun(g_nu,   g_b,tau) in

    let t3 = 2.0 * g_rho * g_nu * g_sigma / (g_a * g_b)*
             ( tau - ba - bb + b_fun(g_a+g_b, tau) )

        in {t1+t2+t3, ba, bb}

///////////////////////////////////////////////////////////////////
// the first parameter `genome' is the five-genes genome used in
//     the genetic algorithms that models the interest rate
// the other parameter are times: today, maturity, and the
//      lower and upper bound of the considered time interval
//
// the result is: x drift term in tmat-forward measure
///////////////////////////////////////////////////////////////////

fun real bigmx( {real,real,real,real,real} genome,
                real today, real tmat, real s, real t
              ) =
    let {a, b, rho, nu, sigma} = genome   in

    let ts    = date_act_365(t,    s)     in
    let tmatt = date_act_365(tmat, t)     in

    let tmat0 = date_act_365(tmat, today) in
    let tmats = date_act_365(tmat, s)     in
    let t0    = date_act_365(t,    today) in
    let s0    = date_act_365(s,    today) in

    let tmp1  = (sigma*sigma)/(a*a)+(sigma*rho*nu)/(a*b)          in
    let tmp2  = 1.0 - exp(-a * ts)                                in
    let tmp3  = sigma * sigma / (2.0 * a * a)                     in
    let tmp4  = rho * sigma * nu / (b * (a + b))                  in
    let tmp5  = exp(-a * tmatt) - exp(-a * (tmats + ts))          in
    let tmp6  = exp(-b * tmatt) - exp(-b*tmat0 - a*t0 + (a+b)*s0)

        in tmp1 * tmp2 - ( tmp3 * tmp5 ) - ( tmp4 * tmp6 )


///////////////////////////////////////////////////////////////////
// the first parameter `genome' is the five-genes genome used in
//     the genetic algorithms that models the interest rate
// the other parameter are times: today, maturity, and the
//      lower and upper bound of the considered time interval
//
// the result is: y drift term in tmat-forward measure
///////////////////////////////////////////////////////////////////

fun real bigmy( {real,real,real,real,real} genome,
                real today, real tmat, real s, real t
              ) =
    let {a, b, rho, nu, sigma} = genome   in

    let ts    = date_act_365(t,    s)     in
    let tmatt = date_act_365(tmat, t)     in
    let tmat0 = date_act_365(tmat, today) in
    let tmats = date_act_365(tmat, s)     in
    let t0    = date_act_365(t,    today) in
    let s0    = date_act_365(s,    today) in

    let tmp1  = nu*nu/(b*b)+sigma*rho*nu/(a*b)     in
    let tmp2  = 1.0 - exp(-b * ts)                 in
    let tmp3  = nu * nu / (2.0 * b * b)            in
    let tmp4  = sigma * rho * nu / (a * (a + b))         in
    let tmp5  = exp(-b * tmatt) - exp(-b * (tmats + ts)) in
    let tmp6  = exp(-a * tmatt) - exp(-a*tmat0 - b*t0 + (a+b)*s0)

        in tmp1 * tmp2 - ( tmp3 * tmp5 ) - ( tmp4 * tmp6 )

///////////////////////////////////////////////////////////////////
// the first  parameter `today' is today (and very different from tomorrow)
// the second parameter is the swaption
// the third  parameter is the implied volability
//
// the result is: the swaption's price
///////////////////////////////////////////////////////////////////

fun real black_price(real today, {real,real,real} swaption, real vol ) =
    let {maturity, swap_sched, strike} =
                        extended_swaption_of_swaption( swaption ) in

    let sqrtt = date_act_365(maturity, today) in

    // morally equivalent to `swap_schedule2lvl(swap_schedule)' but in map-reduce form!!
    let {swap_sched1, swap_sched2} = unzip(swap_sched)                                    in
    let n = size(0, swap_sched) in
    let swap_sched_new = zip(replicate(n, 0.0), swap_sched1, swap_sched2)  in
    let {lvl,t0,tn} = reduce(accumSched, {0.0, max_date(), min_date()},  swap_sched_new)  in

    let s0 = (zc(t0) - zc(tn)) / lvl                      in
    let d1 = log(s0/strike) / (vol*sqrtt) + 0.5*vol*sqrtt in
    let d2 = d1-vol*sqrtt                                 in

        lvl * ( s0*uGaussian_P(d1) - strike*uGaussian_P(d2) )


/////////////////////////////////////////////////
/////////////////////////////////////////////////
/// Testing g2pp minus the main function,
///    i.e., pricer_of_swaption
/////////////////////////////////////////////////
/////////////////////////////////////////////////

fun int main_g2pp_header() =
    let today = 9000.0    in
    let tmat  = 18000.0   in
    let s     = 400000.0  in
    let t     = 9000000.0 in

    ///////////////////////////////////////////
    // testing b_fun, bigv, bigmx, bigmy
    ///////////////////////////////////////////

    let res_b_fun = b_fun(3.24, 1.362) in
    let res_bigv  = bigv ({0.02, 0.02, 0.0, 0.01, 0.04}, 1.12)              in
    let res_bigmx = bigmx({0.02, 0.02, 0.0, 0.01, 0.04}, today, tmat, s, t) in
    let res_bigmy = bigmy({0.02, 0.02, 0.0, 0.01, 0.04}, today, tmat, s, t) in

    let tmp = trace("b_fun test: ") in
    let tmp =   if( equal(res_b_fun, 0.30490117) )
                then trace(" SUCCESS! ")
                else trace(" fails! ") in
    let tmp = trace(res_b_fun)    in
    let tmp = trace("\n\n")       in

    let tmp = trace("bigv test: ")      in
    let {tmp1, tmp2, tmp3} = res_bigv   in
    let tmp =   if( equal(tmp1, 7.8288965347e-4) && equal(tmp2, 1.107549139) && equal(tmp3, 1.107549139)  )
                then trace(" SUCCESS! ")
                else trace(" fails! ") in
    let tmp = trace(res_bigv)    in
    let tmp = trace("\n\n")       in

    let tmp = trace("bigmx test: ") in
    let tmp =   if( equal(res_bigmx, -0.2356067470979) )
                then trace(" SUCCESS! ")
                else trace(" fails! ") in
    let tmp = trace(res_bigmx)    in
    let tmp = trace("\n\n")       in


    let tmp = trace("bigmy test: ") in
    let tmp =   if( equal(res_bigmy, -0.01472542169362) )
                then trace(" SUCCESS! ")
                else trace(" fails! ") in
    let tmp = trace(res_bigmy)    in
    let tmp = trace("\n\n")       in

    //////////////////////////////////
    // testing extended_swaption_of_swaption
    // The right value to test against is "654.142965".
    // However, because "today is not tomorrow", i.e., the date module
    // is very approximatively implemented, we test against ``655.250458''
    //////////////////////////////////

    let swaption = {10.0, 6.0, 4.0}     in
    let maturity = 22094640.0           in
    let strike   = 0.030226283149239714 in
    let swap_schedule = [{22094640.0, 22355280.0}, {22355280.0, 22620240.0}, {22620240.0, 22880880.0}, {22880880.0, 23145840.0}, {23145840.0, 23407920.0}, {23407920.0, 23672880.0}, {23672880.0, 23933520.0}, {23933520.0, 24198480.0}] in
    let {res_mat, res_swap_schd, res_strike} = extended_swaption_of_swaption(swaption) in
    let mat_ok    = equal(maturity, res_mat   ) in
    let strike_ok = equal(strike,   res_strike) in
    let sched_ok  = reduce( op &&, True,
                            map(    fn bool({{real,real},{real,real}} z) =>
                                        let {{x1,x2},{y1,y2}} = z in equal(x1,y1) && equal(x2,y2),
                                    zip(swap_schedule, res_swap_schd)
                            )
                          ) in
    let tmp = trace("Testing extended_swaption_of_swaption: ")          in
    let tmp =   if(mat_ok && strike_ok && sched_ok)
                then trace("SUCCESS! ")
                else trace("FAILS! ") in
    let tmp = trace("\n\tmaturity: ") in let tmp = trace(res_mat      ) in
    let tmp = trace("\n\tstrike: ")   in let tmp = trace(res_strike   ) in
    let tmp = trace("\n\tswapsched: ")in let tmp = trace(res_swap_schd) in

    //////////////////////////////////
    // testing black_price
    // The right value to test against is "654.142965".
    // However, because "today is not tomorrow", i.e., the date module
    // is very approximatively implemented, we test against ``654.1689526995502''
    //////////////////////////////////

    let vol      = 0.2454           in
    let swaption = {10.0, 6.0, 4.0} in
    let black_price_res = black_price(today(), swaption, vol) * 10000.0 in

    let tmp = trace("\n\nTesting Black Price: ") in
    let tmp =   if( equal(black_price_res, 654.1429648454) )
                then trace(" SUCCESS! ")
                else trace(" FAILS! ") in
    let tmp = trace(black_price_res)   in
    let tmp = trace("\n\n")            in


        33



///////////////////////////
//// test also the other ones via:
////
//// assert "%.6f" % b_fun(3.24,1.362)=="0.304901"
////# bigv [g_a=0.02; g_b=0.02; g_sigma=0.04; g_nu=0.01; g_rho=0.] 1.12
////assert "%.6f %.6f %.6f" % bigv(a=0.02, b=0.02, sigma=0.04, nu=0.01, rho=0.0, tau=1.12) == "0.000783 1.107549 1.107549"
////# bigmx [g_a=0.02; g_b=0.02; g_sigma=0.04; g_nu=0.01; g_rho=0.] 9000 18000 400000 9000000
////assert "%.6f" % bigmx(a=0.02, b=0.02, sigma=0.04, nu=0.01, rho=0.0, today=9000,tmat=18000,s=400000,t=9000000) == "-0.235607"
////# bigmy [g_a=0.02; g_b=0.02; g_sigma=0.04; g_nu=0.01; g_rho=0.] 9000 18000 400000 9000000
////assert "%.6f" % bigmy(a=0.02, b=0.02, sigma=0.04, nu=0.01, rho=0.0, today=9000,tmat=18000,s=400000,t=9000000) == "-0.014725"
////# extended_swaption_of_swaption today zc [swaption_maturity_in_year = 10; swap_term_in_year = 4; swap_frequency = Freq_6M]
////assert extended_swaption_of_swaption(today=today,zc=zc,swaption=['swaption_maturity_in_year': 10, 'swap_term_in_year': 4, 'swap_frequency': 6]) == ['maturity':22094640, 'swap_schedule':[(22094640, 22355280), (22355280, 22620240), (22620240, 22880880), (22880880, 23145840), (23145840, 23407920), (23407920, 23672880), (23672880, 23933520), (23933520, 24198480)], 'strike':0.030226283149239714]
////# black_price today zc [swaption_maturity_in_year = 10; swap_term_in_year = 4; swap_frequency = Freq_6M] 0.2454
////swaption=['swaption_maturity_in_year': 10, 'swap_term_in_year': 4, 'swap_frequency': 6]
////assert "%.6f" % (black_price(today=today,zc=zc,swaption=swaption,vol=0.2454)*10000) == "654.142965"
////////////////////////////



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
/// Main function of Module G2PP: pricer_of_swaption    ///
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

//////////////
// def a_fun(end_date,a,b,rho,nu,sigma,today,maturity,zc_mat,v0_mat):
//   # Brigo and Mercurio: defined top p. 148
//   v0_end,dummyA,dummyB=bigv(a,b,rho,nu,sigma,tau=Date.act_365(end_date,today))
//   vt_end,ba,bb=bigv(a,b,rho,nu,sigma,tau=Date.act_365(end_date,maturity))
//   res=zc(end_date)/zc_mat*N.exp(0.5*(vt_end-v0_end+v0_mat))
//   return res,ba,bb
//////////////

fun real
pricer_of_swaption( real                       today,
                    {real,real,real}           swaption,
                    {real,real,real,real,real} genome,
                    [real]                     x_quads,
                    [real]                     w_quads
                  ) =
    let swaption = extended_swaption_of_swaption(swaption) in
    let {maturity, schedulei, strike} = swaption           in

    // COSMIN: TODO: add schedulei as parameter for map!!!!!
    // COSMIN: PERFECT EXAMPLE FOR THE PAPER!!!!!
    let n_schedi = size(0, schedulei)                      in
    let ci = map(   fn real (int i) =>
                        let {d_beg,d_end} = schedulei[i]   in
                        let tau = date_act_365(d_end,d_beg)in
                        if(i == n_schedi-1)
                        then 1.0 + tau*strike
                        else       tau*strike
                    , iota(n_schedi)
                )                                       in
//
    let tmat0    = date_act_365( maturity, today() )    in
    let {v0_mat, dummyA, dummyB} = bigv( genome, tmat0) in
    let zc_mat   = zc(maturity)                         in
//
    let {a,b,rho,nu,sigma} = genome                     in
    let sigmax = sigma * sqrt( b_fun(2.0*a, tmat0) )    in
    let sigmay = nu    * sqrt( b_fun(2.0*b, tmat0) )    in
    let rhoxy  = (rho * sigma * nu) / (sigmax * sigmay)
                    * b_fun(a+b, tmat0)                 in

    let rhoxyc = 1.0 - rhoxy * rhoxy                    in
    let rhoxycs= sqrt( rhoxyc )                         in
    let t2     = rhoxy / (sigmax*rhoxycs)               in
    let sigmay_rhoxycs = sigmay * rhoxycs               in
    let t4     = (rhoxy * sigmay) / sigmax              in
//
    let mux    = -bigmx( genome, today, maturity, today, maturity ) in
    let muy    = -bigmy( genome, today, maturity, today, maturity ) in
//
    let {scheduleix, scheduleiy} = unzip(schedulei) in
//
    let {bai, bbi, aici, log_aici, t1_cst, scale} = unzip (
            map(fn {real,real,real,real,real,real} ({real,real} dc) =>
                    let {end_date, ci} = dc                                   in

                  // Begin Brigo and Mercurio: defined top p. 148
                    let {v0_end, dummyA, dummyB} =
                            bigv( genome, date_act_365(end_date, today   ) )  in

                    let {vt_end, bai, bbi} =
                            bigv( genome, date_act_365(end_date, maturity) )  in

                    let aa = zc(end_date) / zc_mat *
                                exp( 0.5 * (vt_end-v0_end+v0_mat) )           in
                  // END Brigo and Mercurio: defined top p. 148

                    let aici = ci * aa                                        in
                    let log_aici = log(aici)                                  in

                    let t3 = muy - 0.5*rhoxyc*sigmay*sigmay*bbi               in
                    let cst= bbi * (mux*t4 - t3)                              in
                    let t1_cst = aici * exp(cst)                              in
                    let scale  = -(bai + bbi*t4)                              in
                        {bai, bbi, aici, log_aici, t1_cst, scale}

                , zip(scheduleiy, ci)
            )
        )                                                               in

    let babaici = zip(bai, bbi, aici, log_aici)                         in
    let scals   = {b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, mux, muy} in

    let eps = 0.5 * sigmax                                     in
    let f   = exactYhat( n_schedi, scals, babaici, mux       ) in
    let g   = exactYhat( n_schedi, scals, babaici, mux + eps ) in
    let h   = exactYhat( n_schedi, scals, babaici, mux - eps ) in
    let df  = 0.5 * ( g - h ) / eps  in

    let sqrt2sigmax = sqrt(2.0) * sigmax                       in

    let tmps = map(
                    fn real ( {real,real} quad ) =>
                        let {x_quad, w_quad} = quad       in
                        let x = sqrt2sigmax*x_quad + mux  in

                        ///////////////////////////////////////////
                        // BEGIN function integrand(x) inlined
                        ///////////////////////////////////////////
                        let tmp = (x - mux) / sigmax      in
                        let t1  = exp( -0.5 * tmp * tmp ) in

                        let yhat_x = f + df*(x - mux)     in
                        let h1  = ( (yhat_x - muy) / sigmay_rhoxycs ) - t2*( x - mux ) in

                        let tmps= map(  fn real ( {real,real,real} bbit1cstscale ) =>
                                            let {bbii, t1_csti, scalei} = bbit1cstscale in
                                            let h2 = h1 + bbii * sigmay_rhoxycs in
                                                t1_csti * exp(scalei*x) * uGaussian_P(-h2)
                                        , zip(bbi, t1_cst, scale)
                                     ) in
                        let accum = reduce(op +, 0.0, tmps) in
                        let integrand_res = t1 * ( uGaussian_P(-h1) - accum )
                        ///////////////////////////////////////////
                        // END   function integrand(x) inlined
                        ///////////////////////////////////////////

                        in w_quad * integrand_res

                  , zip(x_quads, w_quads)
                  )                        in
    let sum = reduce(op +, 0.0, tmps)      in
            zc_mat * ( sum / sqrt( pi() ) )


//////////////////////////
// Root finder
//////////////////////////
fun real exactYhat( int n_schedi,
                    {real,real,real,real,real,real,real,real} scals,
                    [{real,real,real,real}] babaicis,
                    real x
                  ) =
    // ugaussian_Pinv(k)=1.0e-4
    let k=-3.71901648545568                                   in


    let uplos = map(  fn {real,real} ({real,real,real,real} babaici) =>
                        let {bai,bbi,aici,log_aici} = babaici in
                        let baix                    = bai * x in
                            {   aici * exp( -baix ),
                                (log_aici-baix) / bbi
                            }
                      , babaicis
                   )                                          in
    let {ups, los} = unzip(uplos)             in
    let up = reduce(op +, 0.0, ups)           in
    let lo = reduce(MAX, -infinity(), los)    in

    let {bai, bbi, aici, log_aici} = unzip(babaicis) in

    if(n_schedi == 1)
    then lo
    else
         let log_s = log(up)                  in
         let tmp   = log_s / bbi[n_schedi-1]  in
         let up    = if( tmp<= 0.0 ) then tmp
                     else
                       let tmp = log_s/bbi[0] in
                       if(0.0<= tmp) then tmp
                       else -infinity()       in

         let yl = lo - epsilon()              in
         let yu = up + epsilon()              in

         let {b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, mux, muy} = scals   in

         // y01 x = y0, y1 / phi(h_i(x, y0)) <= epsilon, 1 - phi(h_i(x, y1)) <= epsilon
         let y0 = sigmay * (rhoxy*(x-mux)/sigmax+k*rhoxycs) - rhoxyc/b + muy  in
         let y1 = sigmay * (rhoxy*(x-mux)/sigmax-k*rhoxycs) + muy             in

         if      (y1 <= yl) then y1 + 1.0  // yhat is greater than y1 => 1 - phi(h_i(x, yhat)) < epsilon
         else if (yu <= y0) then y0 - 1.0  // yhat is lower than y0 => phi(h_i(x, yhat)) < epsilon)
         else
              // `scales' is the same as `ups', however, if this branch
              // is not oftenly taken, it might make sense to duplicate computation,
              // since if the previous `ups' is in a map-reduce pattern!
              let scales  = map(  fn real ( {real,real} baiaici) =>
                                    let {bai,aici} = baiaici in
                                    aici * exp( -bai * x )
                                  , zip(bai, aici)
                               )        in

              let root_lb = MAX(yl, y0) in
              let root_ub = MIN(yu, y1) in
              let {root, iteration, error} =
                    rootFinding_Brent(1, zip(scales, bbi), root_lb, root_ub, 1.0e-4, 1000) in

              if      ( error == -infinity() ) then y0 - 1.0
              else if ( error ==  infinity() ) then y1 + 1.0
              else                                 root

////////////////////////////////////////////////////////////////
//              // WRONG: caused serious bugs...
//              if(  iteration < 1000 ) then root
//              else if( error < 0.0  ) then y0 - 1.0
//                                      else y1 + 1.0


//////////////////////////////////////////////////////
//// Test pricer_of_swaption:
////
//// params=params2dict(a = 0.02453, b = 0.98376, sigma = 0.02398, nu = 0.11830, rho = -0.82400)
//// swaption=['swaption_maturity_in_year': 10, 'swap_term_in_year': 4, 'swap_frequency': 6]
//// assert "%.3f" % (1e4*pricer_of_swaption(today=today,zc=zc,swaption=swaption,params=params)) == "657.822"
//// swaption=['swaption_maturity_in_year': 30, 'swap_term_in_year': 30, 'swap_frequency': 6]
//// assert "%.3f" % (1e4*pricer_of_swaption(today=today,zc=zc,swaption=swaption,params=params)) == "1902.976"
////
//////////////////////////////////////////////////////


fun int main_pricer_of_swaption([real] x_quads, [real] w_quads) =

    // (a,b,rho,nu,sigma) = genome
    let genome   = {0.02453, 0.98376, -0.82400, 0.11830, 0.02398}      in

    // (maturity, frequency, term) = swaption
    let swaption = {10.0, 6.0, 4.0}                                    in

    let price1   = 1.0e4*pricer_of_swaption(today(), swaption, genome, x_quads, w_quads) in

    let tmp = trace("Pricer_of_swaption test: ") in
    let tmp =   if( equal(price1, 657.82158867845) )
                then trace(" SUCCESS! ")
                else let tmp = trace(" FAILS! should be: ") in let tmp = trace(657.822) in trace(" is ") in
    let tmp = trace(price1)    in
    let tmp = trace("\n\n")    in


    // (maturity, frequency, term) = swaption
    let swaption = {30.0, 6.0, 30.0}                                   in

    let price2   = 1.0e4*pricer_of_swaption(today(), swaption, genome, x_quads, w_quads) in

    let tmp = trace("Pricer_of_swaption test: ") in
    let tmp =   if( equal(price2, 1902.97628191498) )
                then trace(" SUCCESS! ")
                else let tmp = trace(" FAILS! should be: ") in let tmp = trace(1902.976) in trace(" is ") in
    let tmp = trace(price2)    in
    let tmp = trace("\n\n")    in


    // (maturity, frequency, term) = swaption
    let swaption = {30.0, 6.0, 25.0}                                   in

    let price3   = 1.0e4*pricer_of_swaption(today(), swaption, genome, x_quads, w_quads) in

    let tmp = trace("Pricer_of_swaption test: ") in
    let tmp =   if( equal(price3, 1840.859126408099) )
                then trace(" SUCCESS! ")
                else let tmp = trace(" FAILS! should be: ") in let tmp = trace(1840.859126408099) in trace(" is ") in
    let tmp = trace(price3)    in
    let tmp = trace("\n\n")    in


        33




/////////////////////////////////////////////////////////////////////////////////
//fun real integrand(real x) =
//    let tmp = (x - mux) / sigmax      in
//    let t1  = exp( -0.5 * tmp * tmp ) in
//
//    let yhat_x = f + df*(x - mux)     in
//    let h1  = ( (yhat_x - muy) / sigmay_rhoxycs ) - t2*( x - mux )
//
//    let tmps= map(  fn real ( (real,real,real) bbit1cstscale ) =>
//                        let (bbi, t1_cst, scale) = bbit1cstscale in
//                        let h2 = h1 + bbi * sigmay_rhoxycs in
//                            t1_cst * exp(scale*x) * uGaussian_P(-h2)
//                    , zip(bbi, t1_cst, scale)
//                 ) in
//    let accum = reduce(op +, 0.0, tmps) in
//        t1 * ( uGaussian_P(-h1) - accum )
//
/////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////
////   MAIN CALIBRATION!!!
//////////////////////////////////////////////////////


fun real main([{{real,real,real} , real}] swaptionQuotes,
              // Gaussian quadrature data
              [real] x_quads,
              [real] w_quads) =
    let genome = {0.02453, 0.98376, -0.82400, 0.11830, 0.02398} in
    let prices = map( fn real ( {{real,real,real} , real} swapquote ) =>
                        let {swaption, quote} = swapquote                                in
                        let g2pp_price   = pricer_of_swaption(today(), swaption, genome, x_quads, w_quads) in
                        let market_price = black_price       (today(), swaption, quote ) in

                        // printing
                        let {mat_year, swap_freq, term_year} = swaption                in
                        //let tmp = trace("\n")                                        in
                        //let tmp = trace(trunc( mat_year)) in let tmp = trace("Y")    in
                        //let tmp = trace(trunc(term_year)) in let tmp = trace("Y")    in
                        //let tmp = trace(" Swaption Calibrated Price: ") in let tmp = trace(10000.0*g2pp_price) in
                        //let tmp = trace(" Market Price: ") in let tmp = trace(10000.0*market_price) in

                                //(g2pp_price, market_price)
                        let res = (g2pp_price - market_price) / market_price
                                in res * res

                    , swaptionQuotes
                    )    in

    let rms    = reduce(op +, 0.0, prices)      in
    let numswapts = size(0, swaptionQuotes ) in
    let rms    = 100.0 * sqrt ( rms / toReal(numswapts) ) in

    // printing the error!
    let tmp    = trace("\n\nComputed RMS is: ") in
    let tmp    = trace(rms)                     in
    let tmp    = trace("\n\n END \n\n")         in

                rms



//////////////////////////////////////////////////////
////   Date: Gregorian calendar
//////////////////////////////////////////////////////

fun int MOD(int x, int y) = x - (x/y)*y

fun int hours_in_dayI   () = 24
fun int minutes_in_dayI () = hours_in_dayI() * 60
fun int minutes_to_noonI() = (hours_in_dayI() / 2) * 60

fun real minutes_in_day  () = 24.0*60.0


//                           year*month*day*hour*mins
fun int date_of_gregorian( {int,int,int,int,int} date) =
    let {year, month, day, hour, mins} = date in
    let ym =
        if(month == 1 || month == 2)
        then    ( 1461 * ( year + 4800 - 1 ) ) / 4 +
                  ( 367 * ( month + 10 ) ) / 12 -
                  ( 3 * ( ( year + 4900 - 1 ) / 100 ) ) / 4
        else    ( 1461 * ( year + 4800 ) ) / 4 +
                  ( 367 * ( month - 2 ) ) / 12 -
                  ( 3 * ( ( year + 4900 ) / 100 ) ) / 4 in
    let tmp = ym + day - 32075 - 2444238

            in tmp * minutes_in_dayI() + hour * 60 + mins


fun {int,int,int,int,int}
gregorian_of_date ( int minutes_since_epoch ) =
    let jul = minutes_since_epoch / minutes_in_dayI() in
    let l = jul + 68569 + 2444238 in
    let n = ( 4 * l ) / 146097 in
    let l = l - ( 146097 * n + 3 ) / 4 in
    let i = ( 4000 * ( l + 1 ) ) / 1461001 in
    let l = l - ( 1461 * i ) / 4 + 31 in
    let j = ( 80 * l ) / 2447 in
    let d = l - ( 2447 * j ) / 80 in
    let l = j / 11 in
    let m = j + 2 - ( 12 * l ) in
    let y = 100 * ( n - 49 ) + i + l in

    //let daytime = minutes_since_epoch mod minutes_in_day in
    let daytime = MOD( minutes_since_epoch, minutes_in_dayI() ) in

    if ( daytime == minutes_to_noonI() )

    //then [year = y; month = m; day = d; hour = 12; minute = 0]
    then {y, m, d, 12, 0}

    //else [year = y; month = m; day = d; hour = daytime / 60; minute = daytime mod 60]
    else {y, m, d, daytime / 60, MOD(daytime, 60) }


fun bool check_date(int year, int month, int day) =
    let tmp1 = ( 1 <= day && 1 <= month && month <= 12 && 1980 <= year && year <= 2299 ) in
    let tmp2 = ( day <= 28 ) in

    let tmp3 = if      ( month == 2 )
               then let tmpmod = MOD(year, 100) in
                        ( day == 29 && MOD(year, 4) == 0 && ( year == 2000 || (not (tmpmod == 0)) ) )
               else if ( month == 4 || month == 6 || month == 9 || month == 11 )
                    then ( day <= 30 )
                    else ( day <= 31 )

        in tmp1 && (tmp2 || tmp3)


fun real days_between(real t1, real t2) =
  (t1 - t2) / minutes_in_day()

fun real date_act_365(real t1, real t2) = days_between(t1, t2) / 365.0

fun bool leap(int y) = ( MOD(y,4) == 0  && ( (not (MOD(y,100)==0)) || (MOD(y,400)==0) ) )

fun int end_of_month(int year, int month) =
    if      ( month == 2 && leap(year) )                           then 29
    else if ( month == 2)                                          then 28
    else if ( month == 4 || month == 6 || month == 9 || month == 11 ) then 30
    else                                                               31


fun real add_months ( real date, real rnbmonths ) =
    let nbmonths          = trunc(rnbmonths)                 in
    let {y, m, d, h, min} = gregorian_of_date( trunc(date) ) in
    let m = m + nbmonths                                     in
    let {y, m} = {y + (m-1) / 12, MOD(m-1, 12) + 1}          in
    let {y, m} = if (m <= 0) then {y - 1, m + 12} else {y, m} in
    let resmin = date_of_gregorian ( {y, m, MINI( d, end_of_month(y, m) ), 12, 0} ) in

            toReal(resmin)


fun real add_years(real date, real nbyears) =
        add_months(date, nbyears * 12.0)


//assert max_date==168307199, a.k.a. "2299-12-31T23:59:59"
fun real max_date () = 168307199.0

//assert min_date==3600, a.k.a., "1980-01-01T12:00:00"
fun real min_date () = 3600.0

// Date.of_string("2012-01-01")
fun real today    () = toReal( date_of_gregorian( {2012, 1, 1, 12, 0} ) )


////// Previous, approximate implementation ////////
//fun real date_act_365(real t1, real t2) = (t1 - t2) / 365.0
//fun real add_years   (real d1, real y ) = d1 + y*365.0
//fun real add_months  (real d1, real m ) = d1 + m*30.5
//fun real days_between(real t1, real t2) = t1 - t2
//fun real days_from_mins(real t) = t / (24.0*60.0)
//fun real max_date () = 128.0*365.0 //support 128 years
//fun real min_date () = -1.0


//assert add_months(min_date,1)==48240
//assert add_months(min_date,2)==90000
//assert add_years(min_date,1)==530640
//assert add_years(min_date,5)==2634480
//assert "%.6f" % days_between(max_date,min_date) == "116877.499306"
//assert "%.6f" % act_365(max_date,min_date) == "320.212327"
fun int main_dates() =
    let tmp = trace("Today: ") in let tmp = trace(trunc(today())) in let tmp = trace("\n") in
    let tmp = trace("add_months(min_date,1)==48240")            in
    let tmp = add_months(min_date(), 1.0)                       in
    let bbb = if( equal(tmp, 48240.0) ) then trace("SUCCESS ")
                                        else trace("FAILS ")    in
    let bbb = trace(tmp) in let bbb = trace("\n")               in
//
    let tmp = trace("add_months(min_date,2)==90000")            in
    let tmp = add_months(min_date(), 2.0)                       in
    let bbb = if( equal(tmp, 90000.0) ) then trace("SUCCESS ")
                                        else trace("FAILS ")    in
    let bbb = trace(tmp) in let bbb = trace("\n")               in
//
    let tmp = trace("add_years(min_date,1)==530640")            in
    let tmp = add_years(min_date(), 1.0)                        in
    let bbb = if( equal(tmp, 530640.0) )then trace("SUCCESS ")
                                        else trace("FAILS ")    in
    let bbb = trace(tmp) in let bbb = trace("\n")               in
//
    let tmp = trace("add_years(min_date,5)==2634480")           in
    let tmp = add_years(min_date(), 5.0)                        in
    let bbb = if( equal(tmp, 2634480.0) )then trace("SUCCESS ")
                                         else trace("FAILS ")   in
    let bbb = trace(tmp) in let bbb = trace("\n")               in
//
    let tmp = trace("days_between(max_date(),min_date()) == 116877.499306") in
    let tmp = days_between( max_date(), min_date() )                        in
    let bbb = if( equal(tmp, 116877.499305555) )then trace("SUCCESS ")
                                             else trace("FAILS ")           in
    let bbb = trace(tmp) in let bbb = trace("\n")                           in
//
    let tmp = trace("act_365(max_date,min_date) == 320.212327")   in
    let tmp = date_act_365(max_date(), min_date())                in
    let bbb = if( equal(tmp, 320.2123268645) ) then trace("SUCCESS ")
                                           else trace("FAILS ")   in
    let bbb = trace(tmp) in let bbb = trace("\n")                 in

            33
