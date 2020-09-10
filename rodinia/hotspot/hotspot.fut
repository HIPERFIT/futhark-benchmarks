-- Code and comments based on
-- https://github.com/kkushagra/rodinia/blob/master/openmp/hotspot/hotspot_openmp.cpp
--
-- ==
-- tags { futhark-c futhark-opencl }
-- compiled input @ data/64.in
-- output @ data/64.out
--
-- input @ data/512.in
-- output @ data/512.out
--
-- input @ data/1024.in
-- output @ data/1024.out

-- Maximum power density possible (say 300W for a 10mm x 10mm chip)
let max_pd: f32 = 3.0e6

-- Required precision in degrees
let precision: f32 = 0.001

let spec_heat_si: f32 = 1.75e6

let k_si: f32 = 100.0

-- Capacitance fitting factor
let factor_chip: f32 = 0.5

-- Chip parameters
let t_chip: f32 = 0.0005
let chip_height: f32 = 0.016
let chip_width: f32 = 0.016

-- Ambient temperature assuming no package at all
let amb_temp: f32 = 80.0

-- Single iteration of the transient solver in the grid model.
-- advances the solution of the discretized difference equations by
-- one time step
let single_iteration [row][col]
                    (temp: [row][col]f32, power: [row][col]f32,
                     cap: f32, rx: f32, ry: f32, rz: f32,
                     step: f32): [row][col]f32 =
  map (\r ->
         map (\c ->
                #[unsafe]
                let temp_el =  temp[r,c]
                let delta =
                  (step / cap) *
                  (power[r,c] +
                   (if r == 0 && c == 0 -- Corner 1
                    then (temp[r,c+1] - temp_el) / rx +
                         (temp[r+1,c] - temp_el) / ry
                    else if r == 0 && c == col-1 -- Corner 2
                    then (temp[r,c-1] - temp_el) / rx +
                         (temp[r+1,c] - temp_el) / ry
                    else if r == row-1 && c == col-1 -- Corner 3
                    then (temp[r,c-1] - temp_el) / rx +
                         (temp[r-1,c] - temp_el) / ry
                    else if r == row-1 && c == 0 -- Corner 4
                    then (temp[r,c+1] - temp_el) / rx +
                         (temp[r-1,c] - temp_el) / ry
                    else if r == 0 -- Edge 1
                    then (temp[r,c+1] + temp[r,c-1] - 2*temp_el) / rx +
                         (temp[r+1,c] - temp_el) / ry
                    else if c == col-1 -- Edge 2
                    then (temp[r,c-1] - temp_el) / rx +
                         (temp[r+1,c] + temp[r-1,c] - 2*temp_el) / ry
                    else if r == row-1 -- Edge 3
                    then (temp[r,c+1] + temp[r,c-1] - 2*temp_el) / rx +
                         (temp[r-1,c] - temp_el) / ry
                    else if c == 0 -- Edge 4
                    then (temp[r,c+1] - temp_el) / rx +
                         (temp[r+1,c] + temp[r-1,c] - 2*temp_el) / ry
                    else (temp[r,c+1] + temp[r,c-1] - 2 * temp_el) / rx +
                         (temp[r+1,c] + temp[r-1,c] - 2 * temp_el) / ry) +
                   (amb_temp - temp_el) / rz)
                in temp_el + delta)
             (iota col))
      (iota row)

-- Transient solver driver routine: simply converts the heat transfer
-- differential equations to difference equations and solves the
-- difference equations by iterating.
--
-- Returns a new 'temp' array.
entry compute_tran_temp [row][col]
                       (num_iterations: i32) (temp: [row][col]f32) (power: [row][col]f32): [row][col]f32 =
  let grid_height = chip_height / f32.i64(row)
  let grid_width = chip_width / f32.i64(col)
  let cap = factor_chip * spec_heat_si * t_chip * grid_width * grid_height
  let rx = grid_width / (2 * k_si * t_chip * grid_height)
  let ry = grid_height / (2 * k_si * t_chip * grid_width)
  let rz = t_chip / (k_si * grid_height * grid_width)
  let max_slope = max_pd / (factor_chip * t_chip * spec_heat_si)
  let step = precision / max_slope
  in loop temp for _i < num_iterations do
       single_iteration(temp, power, cap, rx, ry, rz, step)

entry ambient_temps (row: i64) (col: i64): [row][col]f32 =
  unflatten row col (replicate (row*col) (amb_temp ))

entry render_frame [row][col] (temp: [row][col]f32): [row][col][3]i8 =
  let hottest = 400f32
  let coldest = amb_temp
  in map (\(temp_r: []f32): [col][3]i8  ->
            map (\(c: f32): [3]i8  ->
                   let c' = f32.min hottest (f32.max c coldest)
                   let intensity = ((c' - coldest) / (hottest - coldest)) * 256
                   in [i8.f32(intensity),
                       i8.f32(intensity/2),
                       i8.f32(intensity/2)])
                temp_r)
         temp

let main [row][col] (num_iterations: i32) (temp: [row][col]f32) (power: [row][col]f32): [][]f32 =
  compute_tran_temp num_iterations temp power
