-- | Heston calibration based on a generic least-squares solver.

import "lib/github.com/diku-dk/date/date"
import "lib/github.com/diku-dk/cpprandom/random"
import "distance"
import "least_squares"
import "price_european_calls"

module heston (real: real)
              (rand: rng_engine): {
  type real = real.t
  val heston [num_quotes]:
    i32 -> i32 -> i32 -> i32 ->
    [num_quotes]i32 ->
    [num_quotes]real ->
    [num_quotes]real ->
    (real,i32,real,real,real,real,real)
} = {

type real = real.t
let int (x: i32) = real.i32 x

let heston_parameters_from_vector (x: [5]real) =
  { initial_variance = x[0]
  , long_term_variance = x[1]
  , correlation = x[2]
  , mean_reversion = x[3]
  , variance_volatility = x[4] }

module pricer = price_european_calls real

module real_least_squares = mk_least_squares real rand
module real_distance = absolute_distance real

type quote = {maturity: date, strike: real, quote: real}

let distinct_maturities [n] (dates: [n]date): ([]date, [n]i32) =
  let switched (x: date) (i: i32) =
    i == 0 || #[unsafe] !(same_date x dates[i - 1])
  let switches = map2 switched dates (map i32.i64 (iota n))
  in ((unzip (filter (.0) (zip switches dates))).1,
      map (\x -> x-1) (scan (+) 0 (map i32.bool switches)))

let run_calibration({today,
                     quotes,
                     max_global,
                     np,
                     strike_weight_bandwidth,
                     maturity_weight_x0,
                     maturity_weight_gamma,
                     integral_iterations,
                     variables}) =
  let price_and_vega_of_quote (strike: real) (maturity: date) (quote: real) =
    (let (price, vega) = pricer.bs_call true today (int 1) strike maturity quote
     in (price, real.max (real.f64 1e-1) vega))
  let strike_weight (p: real) (x: real) = real.(exp (p * (log x + int 1 - x)))
  let maturity_weight (x0: real) (mgamma: real) (x: real) =
      (let k = real.(int 1 / (exp(mgamma * x0) - int 1))
       in real.(if x <= x0 then k * (exp(mgamma * x) - int 1) else int 1))
  let weight (strike: real) (mat: date) =
    real.(maturity_weight maturity_weight_x0 maturity_weight_gamma (f64 (diff_dates today mat)) *
          strike_weight strike_weight_bandwidth strike)


  let (maturity_dates, quotes_to_maturities) =
    distinct_maturities (map (.maturity) quotes)
  let weights = map (\{maturity, strike, quote=_} -> weight strike maturity) quotes
  let prices_and_vegas = map (\{maturity, strike, quote} ->
                              price_and_vega_of_quote strike maturity quote) quotes

  let day_count_fractions =
    map real.f64 (map (diff_dates today) maturity_dates)
  let gauss_laguerre_coefficients =
    pricer.gauss_laguerre_coefficients integral_iterations

  let objective (x: []real): real =
    let heston_parameters = heston_parameters_from_vector x
    let x_prices = pricer.price_european_calls
                   gauss_laguerre_coefficients
                   false (int 1) (int 1) (int 1)
                   heston_parameters
                   day_count_fractions
                   (map2 (\i q -> {maturity=i, strike=q.strike}) quotes_to_maturities quotes)
    in real_distance.distance
         (map2 (\(p,v) w -> real.(w * p / v)) prices_and_vegas weights)
         (map3 (\w (_, v) p -> real.(w * p / v)) weights prices_and_vegas x_prices)

  in real_least_squares.least_squares variables objective max_global np
                                      (i32.i64 (length quotes))

let date_of_int(x: i32) =
  let d = x%100
  let m = (x/100)%100
  let y = x/10000
  in date_of_triple (y, m, d)

let default_variables: [5]real_least_squares.optimization_variable =
  [real_least_squares.optimize_value
   {lower_bound =  real.f64 1e-6, initial_value = real.f64 4e-2, upper_bound = int 1},
   real_least_squares.optimize_value
   {lower_bound =  real.f64 1e-6, initial_value = real.f64 4e-2, upper_bound = int 1},
   real_least_squares.optimize_value
   {lower_bound =  int (-1), initial_value = real.f64 (-0.5), upper_bound = int 0},
   real_least_squares.optimize_value
   {lower_bound =  real.f64 1e-4, initial_value = real.f64 1e-2, upper_bound = int 4},
   real_least_squares.optimize_value
   {lower_bound =  real.f64 1e-4, initial_value = real.f64 0.4, upper_bound = int 2}
  ]

let heston [num_quotes]
           (max_global: i32)
           (num_points: i32)
           (np: i32)
           (today: i32)
           (quotes_maturity: [num_quotes]i32)
           (quotes_strike: [num_quotes]real)
           (quotes_quote: [num_quotes]real) =
  let result =
    run_calibration { today = date_of_int today
                    , quotes = map3 (\m k q -> {maturity = date_of_int m, strike = k, quote = q})
                      quotes_maturity quotes_strike quotes_quote
                    , max_global = max_global
                    , np = np
                    , strike_weight_bandwidth = int 0
                    , maturity_weight_x0 = int 0
                    , maturity_weight_gamma = int 1
                    , integral_iterations = if num_points == 10 then pricer.ten else pricer.twenty
                    , variables = default_variables
                    }
  let { initial_variance,
        long_term_variance,
        correlation,
        mean_reversion,
        variance_volatility} = heston_parameters_from_vector result.parameters
  in (result.root_mean_squared_error,
      result.num_feval,
      initial_variance,
      long_term_variance,
      mean_reversion,
      variance_volatility,
      correlation
      )
}
