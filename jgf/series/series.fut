-- This benchmark uses double-precision, so note that it will likely
-- be pretty slow on consumer GPUs.
-- ==
-- nobench input { 4i64 }
-- output { [2.8729524964837996f64, 1.1161046676147885f64,
--           0.34429060398168687f64, 0.1523889870251929f64]
--          [0.0f64, -1.8819691893398025f64, -1.1645642623320955f64,
--          -0.8143461113044296f64] }
-- compiled input @ data/10000.in
-- output @ data/10000.out
-- compiled input @ data/100000.in
-- output @ data/100000.out
-- compiled input @ data/1000000.in
-- output @ data/1000000.out

let thefunction(x: f64, omegan: f64, select: i32): f64 =
  if select == 0 then (x+1.0) ** x
  else if select == 1 then (x+1.0)**x * f64.cos(omegan*x)
  else (x+1.0)**x * f64.sin(omegan*x)

let TrapezoidIntegrate(x0: f64, x1: f64, nsteps: i32, omegan: f64, select: i32): f64 =
  let x = x0
  let dx = (x1-x0) / f64.i32 nsteps
  let rvalue = thefunction(x0, omegan, select) / 2.0
  -- already done one step
  let (_,rvalue) = loop ((x,rvalue)) for _i < nsteps-2 do
                     (let x = x + dx
                      in (x, rvalue + thefunction(x,omegan,select)))
  in (rvalue + thefunction(x1,omegan,select) / 2.0) * dx

let main(array_rows: i64): ([array_rows]f64,[array_rows]f64) =
  let omega = 3.1415926535897932
  let first = [(TrapezoidIntegrate (0.0, 2.0, 1000, 0.0, 0) / 2.0,
                0.0) -- never set in reference implementation
              ]
  let rest = map (\i -> (TrapezoidIntegrate (0.0, 2.0, 1000, omega * f64.i64 i, 1),
                         TrapezoidIntegrate (0.0, 2.0, 1000, omega * f64.i64 i, 2)))
                 (map (1+) (iota (array_rows-1)))
  in unzip (first++rest :> [array_rows](f64, f64))
