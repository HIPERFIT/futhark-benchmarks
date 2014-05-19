fun {int,int,[real],[real],[real],[[real]],[[real]],[[real]],[[real]]} initGrid
    (real s0, real alpha, real nu, real t, int numX, int numY, int numT) =
    let logAlpha = log(alpha) in
    let myTimeline = map(fn real (int i) => t * toReal(i) / (toReal(numT) - 1.0), iota(numT)) in
    let stdX = 20.0 * alpha * s0 * sqrt(t) in
    let stdY = 10.0 * nu         * sqrt(t) in
    let {dx, dy} = {stdX / toReal(numX), stdY / toReal(numY)} in
    let {myXindex, myYindex} = {trunc(s0 / dx), numY / 2} in
    let myX = map(fn real (int i) => toReal(i) * dx - toReal(myXindex) * dx + s0, iota(numX)) in
    let myY = map(fn real (int i) => toReal(i) * dy - toReal(myYindex) * dy + logAlpha, iota(numY)) in
    let xXy = replicate(numX,replicate(numY, 0.0)) in
    let {myMuX, myVarX, myMuY, myVarY} = {xXy, xXy, xXy, xXy} in
    {myXindex, myYindex, myX, myY, myTimeline, myMuX, myVarX, myMuY, myVarY}

fun {[{real,real,real}],[{real,real,real}]} initOperator([real] x) =
    let n = size(0, x) in
    let dxu = x[1] - x[0] in
    let dxl = 0.0 in
    let Dxlow = [{0.0, -1.0 / dxu, 1.0 / dxu}] in
    let Dxxlow = [{0.0, 0.0, 0.0}] in
    let Dxmid = map(fn {real,real,real} (int i) =>
                        let dxl = x[i] - x[i-1] in
                        let dxu = x[i+1] - x[i] in
                            { -dxu/dxl/(dxl+dxu),
                             (dxu/dxl - dxl/dxu)/(dxl+dxu),
                             dxl/dxu/(dxl+dxu)
                            }
                    , map (op + (1), iota(n-2))) in
    let Dxxmid = map(fn {real,real,real} (int i) =>
                        let dxl = x[i] - x[i-1] in
                        let dxu = x[i+1] - x[i] in
                            { 2.0/dxl/(dxl+dxu),
                              -2.0*(1.0/dxl + 1.0/dxu)/(dxl+dxu),
                              2.0/dxu/(dxl+dxu)
                            }
                     , map (op + (1), iota(n-2))) in
    let dxl = x[n-1] - x[n-2] in
    let dxu = 0.0 in
    let Dxhigh = [{-1.0 / dxl, 1.0 / dxl, 0.0}] in
    let Dxxhigh = [{0.0, 0.0, 0.0}] in
    let Dx = concat(concat(Dxlow, Dxmid), Dxhigh) in
    let Dxx = concat(concat(Dxxlow, Dxxmid), Dxxhigh) in
    {Dx, Dxx}

fun real max(real x, real y) = if y < x then x else y
fun int maxInt(int x, int y) = if y < x then x else y

fun *[[real]] setPayoff(real strike, [real] myX, [real] myY) =
    let n = size(0, myY) in
    map(fn *[real] (real xi) => replicate(n, max(xi-strike,0.0)), myX)

// Returns new myMuX, myVarX, myMuY, myVarY.
fun {*[[real]] , *[[real]] , *[[real]] , *[[real]]} updateParams
    ([real] myX, [real] myY, [real] myTimeline, int g, real alpha, real beta, real nu) =
    unzip (map(fn {*[real],*[real],*[real],*[real]} (real xi) =>
           unzip (map (fn {real,real,real,real} (real yj) =>
                  {0.0,
                   exp(2.0*(beta*log(xi) + yj - 0.5*nu*nu*myTimeline[g])),
                   0.0,
                   nu * nu}, myY)), myX))

fun {*[real],[real]} tridag
    ([real] a, [real] b, [real] c, [real] r, int n) =
    let bet = 1.0/b[0] in
    let {u, uu} = {replicate(n,0.0), replicate(n,0.0)} in
    let u[0] = r[0] * bet in
    loop ({u, uu, bet}) =
      for j < n-1 do
        let j = j + 1 in
        let uu[j] = c[j-1]*bet in
        let bet = 1.0/(b[j] - a[j] * uu[j]) in
        let u[j] = (r[j] - a[j]*u[j-1]) * bet in
        {u, uu, bet} in
    loop (u) = for j < n - 1 do
                 let j = n - 2 - j in
                 let u[j] = u[j] - uu[j+1]*u[j+1] in
                 u
    in {u, uu}

fun [[real]] explicitX
    (int numX, int numY, real dtInv,
     [[real]] myResult, [[real]] myMuX, [{real,real,real}] myDx, [{real,real,real}] myDxx, [[real]] myVarX) =
    map(fn [real] (int j) =>
        map(fn real (int i) =>
            let res = dtInv*myResult[i,j] in
            let {{dx0, dx1, dx2}, {dxx0, dxx1, dxx2}} = {myDx[i], myDxx[i]} in
            let res = res +
                      if(i == 0) then 0.0
                      else  0.5 * (myMuX[i,j]*dx0+0.5*myVarX[i,j]*dxx0) * myResult[i-1,j] in
            let res = res + 0.5 * (myMuX[i,j]*dx1+0.5*myVarX[i,j]*dxx1) * myResult[i  ,j] in
            let res = res +
                      if(i == numX-1) then 0.0
                      else  0.5 * (myMuX[i,j]*dx2+0.5*myVarX[i,j]*dxx2) * myResult[i+1,j]
            in res
        , iota(numX))
    , iota(numY))

fun [[real]] explicitY
    (int numX, int numY, real dtInv,
     [[real]] myResult, [[real]] myMuY, [{real,real,real}] myDy, [{real,real,real}] myDyy, [[real]] myVarY) =
    map(fn [real] (int i) =>
        map(fn real (int j) =>
            let res = 0.0 in
            let {{dy0, dy1, dy2}, {dyy0, dyy1, dyy2}} = {myDy[j], myDyy[j]} in
            let res = res +
                      if(j == 0) then 0.0
                      else  (myMuY[i,j]*dy0+0.5*myVarY[i,j]*dyy0) * myResult[i,j-1] in
            let res = res + (myMuY[i,j]*dy1+0.5*myVarY[i,j]*dyy1) * myResult[i,j  ] in
            let res = res +
                      if(j == numY-1) then 0.0
                      else  (myMuY[i,j]*dy2+0.5*myVarY[i,j]*dyy2) * myResult[i,j+1]
            in res
        , iota(numY))
    , iota(numX))

fun *[[real]] rollback
    ([real] myX, [real] myY, [real] myTimeline, *[[real]] myResult,
     [[real]] myMuX, [{real,real,real}] myDx, [{real,real,real}] myDxx, [[real]] myVarX,
     [[real]] myMuY, [{real,real,real}] myDy, [{real,real,real}] myDyy, [[real]] myVarY, int g) =
    let {numX, numY} = {size(0, myX), size(0, myY)} in
    let numZ = maxInt(numX, numY) in
    let dtInv = 1.0/(myTimeline[g+1]-myTimeline[g]) in
    let u = explicitX(numX, numY, dtInv, myResult, myMuX, myDx, myDxx, myVarX) in
    let v = explicitY(numX, numY, 0.0, myResult, myMuY, myDy, myDyy, myVarY) in
    let u = map(fn [real] ([real] us, [real] vs) => map(op +, zip(us, vs)),
                zip(u, transpose(v))) in
    let u = map(fn [real] ({[real],[{real,real,real}],[{real,real,real}],[real],[real]} t) =>
                let {uj, myDx, myDxx, myMuX, myVarX} = t in
                let {a,b,c} = unzip(map(fn {real,real,real} ({{real,real,real},{real,real,real},real,real} tt) =>
                                        let {myDx, myDxx, myMuX, myVarX} = tt in
                                        let {dx0,  dx1,  dx2} = myDx  in
                                        let {dxx0, dxx1, dxx2} = myDxx in
                                        {-0.5*(myMuX*dx0 + 0.5*myVarX*dxx0),
                                         dtInv - 0.5*(myMuX*dx1 + 0.5*myVarX*dxx1),
                                         -0.5*(myMuX*dx2+0.5*myVarX*dxx2)}
                                    , zip(myDx, myDxx, myMuX, myVarX))) in
                let {uj, yy} = tridag(a,b,c,uj,numX) in uj,
            zip(u, replicate(numY, myDx), replicate(numY, myDxx), transpose(myMuX), transpose(myVarX))) in
    let myResult =
            map(fn *[real] ({[real],[real],[{real,real,real}],[{real,real,real}],[real],[real]} t) =>
                let {ui, vi, myDy, myDyy, myMuY, myVarY} = t in
                let {a,b,c} = unzip(map(fn {real,real,real} ({{real,real,real},{real,real,real},real,real} tt) =>
                                        let {myDy, myDyy, myMuY, myVarY} = tt in
                                        let {dy0,  dy1,  dy2} = myDy  in
                                        let {dyy0, dyy1, dyy2} = myDyy in
                                        {-0.5*(myMuY*dy0+0.5*myVarY*dyy0),
                                         dtInv - 0.5*(myMuY*dy1+0.5*myVarY*dyy1),
                                         -0.5*(myMuY*dy2+0.5*myVarY*dyy2)}
                                    , zip(myDy, myDyy, myMuY, myVarY))) in
                let y = map(fn real ({real,real} uv) =>
                            let {u,v} = uv in
                            dtInv * u - 0.5 * v
                        , zip(ui,vi)) in
                let {ri, yy} = tridag(a,b,c,y,numY) in
                    ri
            , zip(transpose(u), v, replicate(numX, myDy), replicate(numX, myDyy), myMuY, myVarY))
    in myResult

fun real value(real s0, real strike, real t, real alpha, real nu, real beta) =
    let {numX, numY, numT} = {32, 16, 32} in //(256, 32, 64) in
    let {myXindex, myYindex, myX, myY, myTimeline, myMuX, myVarX, myMuY, myVarY} =
        initGrid(s0, alpha, nu, t, numX, numY, numT) in
    let {myDx, myDxx} = initOperator(myX) in
    let {myDy, myDyy} = initOperator(myY) in
    let myResult = setPayoff(strike, myX, myY) in
    let indices = map(fn int (int i) => numT-2-i,  iota(numT-1)) in
    loop ({myResult, myMuX, myVarX, myMuY, myVarY}) =
        for i < numT - 1 do
            let i = numT-2-i in
            let {myMuX, myVarX, myMuY, myVarY} =
                updateParams(myX, myY, myTimeline, i, alpha, beta, nu) in
            let myResult = rollback(myX, myY, myTimeline, myResult,
                                    myMuX, myDx, myDxx, myVarX,
                                    myMuY, myDy, myDyy, myVarY, i) in
            {myResult, myMuX, myVarX, myMuY, myVarY} in
    myResult[myXindex,myYindex]

fun [real] main (int outer_loop_count, real s0, real strike, real t, real alpha, real nu, real beta) =
    let strikes = map(fn real (int i) => 0.001*toReal(i), iota(outer_loop_count)) in
    let res = map(fn real (real x) => value(s0, x, t, alpha, nu, beta), strikes) in
    res
