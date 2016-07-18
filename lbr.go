package lbr

import (
	"math"
)

const TWO_PI float64 = 6.283185307179586476925286766559005768394338798750
const SQRT_PI_OVER_TWO float64 = 1.253314137315500251207882642405522626503493370305
const SQRT_THREE float64 = 1.732050807568877293527446341505872366942805253810
const SQRT_ONE_OVER_THREE float64 = 0.577350269189625764509148780501957455647601751270
const TWO_PI_OVER_SQRT_TWENTY_SEVEN float64 = 1.209199576156145233729385505094770488189377498728
const PI_OVER_SIX float64 = 0.523598775598298873077107230546583814032861566563
const ONE_OVER_SQRT_TWO float64 = 0.7071067811865475244008443621048490392848359376887
const ONE_OVER_SQRT_TWO_PI float64 = 0.3989422804014326779399460599343818684758586311649
const SQRT_TWO_PI float64 = 2.506628274631000502415765284811045253006986740610

var DBL_EPSILON float64 = math.Nextafter(1, 2) - 1
var SQRT_DBL_EPSILON float64 = math.Sqrt(DBL_EPSILON)
var FOURTH_ROOT_DBL_EPSILON float64 = math.Sqrt(SQRT_DBL_EPSILON)
var EIGHTH_ROOT_DBL_EPSILON float64 = math.Sqrt(FOURTH_ROOT_DBL_EPSILON)
var SIXTEENTH_ROOT_DBL_EPSILON float64 = math.Sqrt(EIGHTH_ROOT_DBL_EPSILON)
var SQRT_DBL_MIN float64 = math.Sqrt(math.SmallestNonzeroFloat64)
var SQRT_DBL_MAX float64 = math.Sqrt(math.MaxFloat64)

const DENORMALIZATION_CUTOFF float64 = 0
const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC float64 = -math.MaxFloat64
const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM float64 = math.MaxFloat64

const implied_volatility_maximum_iterations = 2

const asymptotic_expansion_accuracy_threshold = -10

func is_below_horizon(x float64) bool {
	return math.Abs(x) < DENORMALIZATION_CUTOFF
}

func householder_factor(newton, halley, hh3 float64) float64 {
	return (1.0 + 0.5*halley*newton) / (1.0 + newton*(halley+hh3*newton/6.0))
}

func normalised_intrinsic(x, q float64) float64 {
	if q*x <= 0 {
		return 0
	}
	var x2 float64 = x * x
	var q_sign float64 = 1
	if q < 0 {
		q_sign = -1
	}
	if x2 < 98.0*FOURTH_ROOT_DBL_EPSILON {
		return math.Abs(math.Max(q_sign*x*(1+x2*((1.0/24.0)+x2*((1.0/1920.0)+x2*((1.0/322560.0)+(1.0/92897280.0)*x2)))), 0.0))
	}
	var b_max float64 = math.Exp(0.5 * x)
	var one_over_b_max float64 = 1.0 / b_max
	return math.Abs(math.Max(q_sign*(b_max-one_over_b_max), 0.0))
}

func normalised_intrinsic_call(x float64) float64 {
	return normalised_intrinsic(x, 1.0)
}

func asymptotic_expansion_of_normalized_black_call(h, t float64) float64 {
	var e float64 = (t / h) * (t / h)
	var r float64 = ((h + t) * (h - t))
	var q float64 = (h / r) * (h / r)
	var asymptotic_expansion_sum float64 = (2.0 + q*(-6.0E0-2.0*e+3.0*q*(1.0E1+e*(2.0E1+2.0*e)+5.0*q*(-1.4E1+e*(-7.0E1+e*(-4.2E1-2.0*e))+7.0*q*(1.8E1+e*(1.68E2+e*(2.52E2+e*(7.2E1+2.0*e)))+9.0*q*(-2.2E1+e*(-3.3E2+e*(-9.24E2+e*(-6.6E2+e*(-1.1E2-2.0*e))))+1.1E1*q*(2.6E1+e*(5.72E2+e*(2.574E3+e*(3.432E3+e*(1.43E3+e*(1.56E2+2.0*e)))))+1.3E1*q*(-3.0E1+e*(-9.1E2+e*(-6.006E3+e*(-1.287E4+e*(-1.001E4+e*(-2.73E3+e*(-2.1E2-2.0*e))))))+1.5E1*q*(3.4E1+e*(1.36E3+e*(1.2376E4+e*(3.8896E4+e*(4.862E4+e*(2.4752E4+e*(4.76E3+e*(2.72E2+2.0*e)))))))+1.7E1*q*(-3.8E1+e*(-1.938E3+e*(-2.3256E4+e*(-1.00776E5+e*(-1.84756E5+e*(-1.51164E5+e*(-5.4264E4+e*(-7.752E3+e*(-3.42E2-2.0*e))))))))+1.9E1*q*(4.2E1+e*(2.66E3+e*(4.0698E4+e*(2.3256E5+e*(5.8786E5+e*(7.05432E5+e*(4.0698E5+e*(1.08528E5+e*(1.197E4+e*(4.2E2+2.0*e)))))))))+2.1E1*q*(-4.6E1+e*(-3.542E3+e*(-6.7298E4+e*(-4.90314E5+e*(-1.63438E6+e*(-2.704156E6+e*(-2.288132E6+e*(-9.80628E5+e*(-2.01894E5+e*(-1.771E4+e*(-5.06E2-2.0*e))))))))))+2.3E1*q*(5.0E1+e*(4.6E3+e*(1.0626E5+e*(9.614E5+e*(4.08595E6+e*(8.9148E6+e*(1.04006E7+e*(6.53752E6+e*(2.16315E6+e*(3.542E5+e*(2.53E4+e*(6.0E2+2.0*e)))))))))))+2.5E1*q*(-5.4E1+e*(-5.85E3+e*(-1.6146E5+e*(-1.77606E6+e*(-9.37365E6+e*(-2.607579E7+e*(-4.01166E7+e*(-3.476772E7+e*(-1.687257E7+e*(-4.44015E6+e*(-5.9202E5+e*(-3.51E4+e*(-7.02E2-2.0*e))))))))))))+2.7E1*q*(5.8E1+e*(7.308E3+e*(2.3751E5+e*(3.12156E6+e*(2.003001E7+e*(6.919458E7+e*(1.3572783E8+e*(1.5511752E8+e*(1.0379187E8+e*(4.006002E7+e*(8.58429E6+e*(9.5004E5+e*(4.7502E4+e*(8.12E2+2.0*e)))))))))))))+2.9E1*q*(-6.2E1+e*(-8.99E3+e*(-3.39822E5+e*(-5.25915E6+e*(-4.032015E7+e*(-1.6934463E8+e*(-4.1250615E8+e*(-6.0108039E8+e*(-5.3036505E8+e*(-2.8224105E8+e*(-8.870433E7+e*(-1.577745E7+e*(-1.472562E6+e*(-6.293E4+e*(-9.3E2-2.0*e))))))))))))))+3.1E1*q*(6.6E1+e*(1.0912E4+e*(4.74672E5+e*(8.544096E6+e*(7.71342E7+e*(3.8707344E8+e*(1.14633288E9+e*(2.07431664E9+e*(2.33360622E9+e*(1.6376184E9+e*(7.0963464E8+e*(1.8512208E8+e*(2.7768312E7+e*(2.215136E6+e*(8.184E4+e*(1.056E3+2.0*e)))))))))))))))+3.3E1*(-7.0E1+e*(-1.309E4+e*(-6.49264E5+e*(-1.344904E7+e*(-1.4121492E8+e*(-8.344518E8+e*(-2.9526756E9+e*(-6.49588632E9+e*(-9.0751353E9+e*(-8.1198579E9+e*(-4.6399188E9+e*(-1.6689036E9+e*(-3.67158792E8+e*(-4.707164E7+e*(-3.24632E6+e*(-1.0472E5+e*(-1.19E3-2.0*e)))))))))))))))))*q)))))))))))))))))
	var b float64 = ONE_OVER_SQRT_TWO_PI * math.Exp((-0.5 * (h*h + t*t))) * (t / r) * asymptotic_expansion_sum
	return math.Abs(math.Max(b, 0.0))
}

func erfcx_cody(x float64) float64 {
	return math.Exp(x*x) * math.Erfc(x)
}

func normalised_black_call_using_erfcx(h, t float64) float64 {
	var b float64 = 0.5 * math.Exp(-0.5*(h*h+t*t)) * (erfcx_cody(-ONE_OVER_SQRT_TWO*(h+t)) - erfcx_cody(-ONE_OVER_SQRT_TWO*(h-t)))
	return math.Abs(math.Max(b, 0.0))
}

func small_t_expansion_of_normalized_black_call(h, t float64) float64 {
	var a float64 = 1 + h*(0.5*SQRT_TWO_PI)*erfcx_cody(-ONE_OVER_SQRT_TWO*h)
	var w float64 = t * t
	var h2 float64 = h * h
	var expansion float64 = 2 * t * (a + w*((-1+3*a+a*h2)/6+w*((-7+15*a+h2*(-1+10*a+a*h2))/120+w*((-57+105*a+h2*(-18+105*a+h2*(-1+21*a+a*h2)))/5040+w*((-561+945*a+h2*(-285+1260*a+h2*(-33+378*a+h2*(-1+36*a+a*h2))))/362880+w*((-6555+10395*a+h2*(-4680+17325*a+h2*(-840+6930*a+h2*(-52+990*a+h2*(-1+55*a+a*h2)))))/39916800+((-89055+135135*a+h2*(-82845+270270*a+h2*(-20370+135135*a+h2*(-1926+25740*a+h2*(-75+2145*a+h2*(-1+78*a+a*h2))))))*w)/6227020800.0))))))
	var b float64 = ONE_OVER_SQRT_TWO_PI * math.Exp((-0.5 * (h*h + t*t))) * expansion
	return math.Abs(math.Max(b, 0.0))
}

var small_t_expansion_of_normalized_black_threshold float64 = 2 * SIXTEENTH_ROOT_DBL_EPSILON

const norm_cdf_asymptotic_expansion_first_threshold float64 = -10.0

var norm_cdf_asymptotic_expansion_second_threshold float64 = -1 / SQRT_DBL_EPSILON

func norm_pdf(x float64) float64 {
	return ONE_OVER_SQRT_TWO_PI * math.Exp(-.5*x*x)
}

func norm_cdf(z float64) float64 {
	if z <= norm_cdf_asymptotic_expansion_first_threshold {
		var sum float64 = 1.0
		if z >= norm_cdf_asymptotic_expansion_second_threshold {
			var zsqr float64 = z * z
			var i float64 = 1
			var g float64 = 1
			var x, y float64
			var a float64 = math.MaxFloat64
			var lasta float64
			for {
				lasta = a
				x = (4*i - 3) / zsqr
				y = x * ((4*i - 1) / zsqr)
				a = g * (x - y)
				sum -= a
				g *= y
				i++
				a = math.Abs(a)
				if !(lasta > a && a >= math.Abs(sum*DBL_EPSILON)) {
					break
				}
			}
		}
		return -norm_pdf(z) * sum / z
	}
	return 0.5 * math.Erfc(-z*ONE_OVER_SQRT_TWO)
}

func inverse_norm_cdf(u float64) float64 {
	const split1 float64 = 0.425
	const split2 float64 = 5.0
	const const1 float64 = 0.180625
	const const2 float64 = 1.6

	// Coefficients for P close to 0.5
	const A0 float64 = 3.3871328727963666080E0
	const A1 float64 = 1.3314166789178437745E+2
	const A2 float64 = 1.9715909503065514427E+3
	const A3 float64 = 1.3731693765509461125E+4
	const A4 float64 = 4.5921953931549871457E+4
	const A5 float64 = 6.7265770927008700853E+4
	const A6 float64 = 3.3430575583588128105E+4
	const A7 float64 = 2.5090809287301226727E+3
	const B1 float64 = 4.2313330701600911252E+1
	const B2 float64 = 6.8718700749205790830E+2
	const B3 float64 = 5.3941960214247511077E+3
	const B4 float64 = 2.1213794301586595867E+4
	const B5 float64 = 3.9307895800092710610E+4
	const B6 float64 = 2.8729085735721942674E+4
	const B7 float64 = 5.2264952788528545610E+3
	// Coefficients for P not close to 0, 0.5 or 1.
	const C0 float64 = 1.42343711074968357734E0
	const C1 float64 = 4.63033784615654529590E0
	const C2 float64 = 5.76949722146069140550E0
	const C3 float64 = 3.64784832476320460504E0
	const C4 float64 = 1.27045825245236838258E0
	const C5 float64 = 2.41780725177450611770E-1
	const C6 float64 = 2.27238449892691845833E-2
	const C7 float64 = 7.74545014278341407640E-4
	const D1 float64 = 2.05319162663775882187E0
	const D2 float64 = 1.67638483018380384940E0
	const D3 float64 = 6.89767334985100004550E-1
	const D4 float64 = 1.48103976427480074590E-1
	const D5 float64 = 1.51986665636164571966E-2
	const D6 float64 = 5.47593808499534494600E-4
	const D7 float64 = 1.05075007164441684324E-9
	// Coefficients for P very close to 0 or 1
	const E0 float64 = 6.65790464350110377720E0
	const E1 float64 = 5.46378491116411436990E0
	const E2 float64 = 1.78482653991729133580E0
	const E3 float64 = 2.96560571828504891230E-1
	const E4 float64 = 2.65321895265761230930E-2
	const E5 float64 = 1.24266094738807843860E-3
	const E6 float64 = 2.71155556874348757815E-5
	const E7 float64 = 2.01033439929228813265E-7
	const F1 float64 = 5.99832206555887937690E-1
	const F2 float64 = 1.36929880922735805310E-1
	const F3 float64 = 1.48753612908506148525E-2
	const F4 float64 = 7.86869131145613259100E-4
	const F5 float64 = 1.84631831751005468180E-5
	const F6 float64 = 1.42151175831644588870E-7
	const F7 float64 = 2.04426310338993978564E-15

	if u <= 0 {
		return math.Log(u)
	}
	if u >= 1 {
		return math.Log(1 - u)
	}

	var q float64 = u - 0.5
	if math.Abs(q) <= split1 {
		var r float64 = const1 - q*q
		return q * (((((((A7*r+A6)*r+A5)*r+A4)*r+A3)*r+A2)*r+A1)*r + A0) /
			(((((((B7*r+B6)*r+B5)*r+B4)*r+B3)*r+B2)*r+B1)*r + 1.0)
	} else {
		var r float64
		if q < 0 {
			r = u
		} else {
			r = 1.0 - u
		}
		r = math.Sqrt(-math.Log(r))
		var ret float64
		if r < split2 {
			r = r - const2
			ret = (((((((C7*r+C6)*r+C5)*r+C4)*r+C3)*r+C2)*r+C1)*r + C0) /
				(((((((D7*r+D6)*r+D5)*r+D4)*r+D3)*r+D2)*r+D1)*r + 1.0)
		} else {
			r = r - split2
			ret = (((((((E7*r+E6)*r+E5)*r+E4)*r+E3)*r+E2)*r+E1)*r + E0) /
				(((((((F7*r+F6)*r+F5)*r+F4)*r+F3)*r+F2)*r+F1)*r + 1.0)
		}
		if q < 0 {
			return -ret
		} else {
			return ret
		}

	}
}

func normalized_black_call_using_norm_cdf(x, s float64) float64 {
	var h float64 = x / s
	var t float64 = 0.5 * s
	var b_max float64 = math.Exp(0.5 * x)
	var b float64 = norm_cdf(h+t)*b_max - norm_cdf(h-t)/b_max
	return math.Abs(math.Max(b, 0.0))
}

func normalised_black_call(x, s float64) float64 {
	if x > 0 {
		return normalised_intrinsic_call(x) + normalised_black_call(-x, s)
	}
	var ax float64 = math.Abs(x)
	if s <= ax*DENORMALIZATION_CUTOFF {
		return normalised_intrinsic_call(x)
	}
	if x < s*asymptotic_expansion_accuracy_threshold && 0.5*s*s+x < s*(small_t_expansion_of_normalized_black_threshold+asymptotic_expansion_accuracy_threshold) {
		return asymptotic_expansion_of_normalized_black_call(x/s, 0.5*s)
	}
	if 0.5*s < small_t_expansion_of_normalized_black_threshold {
		return small_t_expansion_of_normalized_black_call(x/s, 0.5*s)
	}
	if x+0.5*s*s > s*0.85 {
		return normalized_black_call_using_norm_cdf(x, s)
	}
	return normalised_black_call_using_erfcx(x/s, 0.5*s)
}

func square(x float64) float64 {
	return x * x
}

func normalised_vega(x, s float64) float64 {
	var ax float64 = math.Abs(x)
	if ax <= 0 {
		return ONE_OVER_SQRT_TWO_PI * math.Exp(-0.125*s*s)
	} else if s <= 0 || s <= ax*SQRT_DBL_MIN {
		return 0
	} else {
		return ONE_OVER_SQRT_TWO_PI * math.Exp(-0.5*(square(x/s)+square(0.5*s)))
	}
}

func normalised_black(x, s, q float64) float64 {
	if q < 0 {
		return normalised_black_call(-x, s)
	} else {
		return normalised_black_call(x, s)
	}
}

func Black(F, K, sigma, T, q float64) float64 {
	var intrinsic float64
	if q < 0 {
		intrinsic = math.Abs(math.Max(K-F, 0.0))
	} else {
		intrinsic = math.Abs(math.Max(F-K, 0.0))
	}

	if q*(F-K) > 0 {
		return intrinsic + Black(F, K, sigma, T, -q)
	}
	return math.Max(intrinsic, (math.Sqrt(F)*math.Sqrt(K))*normalised_black(math.Log(F/K), sigma*math.Sqrt(T), q))
}

var minimum_rational_cubic_control_parameter_value float64 = -(1 - math.Sqrt(DBL_EPSILON))
var maximum_rational_cubic_control_parameter_value float64 = 2 / (DBL_EPSILON * DBL_EPSILON)

func is_zero(x float64) bool {
	return math.Abs(x) < math.SmallestNonzeroFloat64
}

func rational_cubic_interpolation(x, x_l, x_r, y_l, y_r, d_l, d_r, r float64) float64 {
	var h float64 = (x_r - x_l)
	if math.Abs(h) <= 0 {
		return 0.5 * (y_l + y_r)
	}
	var t float64 = (x - x_l) / h
	if !(r >= maximum_rational_cubic_control_parameter_value) {
		var t float64 = (x - x_l) / h
		var omt float64 = 1 - t
		var t2 float64 = t * t
		var omt2 float64 = omt * omt

		return (y_r*t2*t + (r*y_r-h*d_r)*t2*omt + (r*y_l+h*d_l)*t*omt2 + y_l*omt2*omt) / (1 + (r-3)*t*omt)
	}
	return y_r*t + y_l*(1-t)
}

func rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_l float64) float64 {
	var h float64 = (x_r - x_l)
	var numerator float64 = 0.5*h*second_derivative_l + (d_r - d_l)

	if is_zero(numerator) {
		return 0
	}
	var denominator float64 = (y_r-y_l)/h - d_l
	if is_zero(denominator) {
		if numerator > 0 {
			return maximum_rational_cubic_control_parameter_value
		} else {
			return minimum_rational_cubic_control_parameter_value
		}
	}
	return numerator / denominator
}

func rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_r float64) float64 {
	var h float64 = (x_r - x_l)
	var numerator float64 = 0.5*h*second_derivative_r + (d_r - d_l)
	if is_zero(numerator) {
		return 0
	}
	var denominator float64 = d_r - (y_r-y_l)/h
	if is_zero(denominator) {
		if numerator > 0 {
			return maximum_rational_cubic_control_parameter_value
		} else {
			return minimum_rational_cubic_control_parameter_value
		}
	}
	return numerator / denominator
}

func minimum_rational_cubic_control_parameter(d_l, d_r, s float64, preferShapePreservationOverSmoothness bool) float64 {
	var monotonic bool = d_l*s >= 0 && d_r*s >= 0
	var convex bool = d_l <= s && s <= d_r
	var concave bool = d_l >= s && s >= d_r
	if !monotonic && !convex && !concave {
		return minimum_rational_cubic_control_parameter_value
	}
	var d_r_m_d_l float64 = d_r - d_l
	var d_r_m_s float64 = d_r - s
	var s_m_d_l float64 = s - d_l
	var r1 float64 = -math.MaxFloat64
	var r2 float64 = r1
	if monotonic {
		if !is_zero(s) {
			r1 = (d_r + d_l) / s
		} else if preferShapePreservationOverSmoothness {

			r1 = maximum_rational_cubic_control_parameter_value
		}
	}
	if convex || concave {
		if !(is_zero(s_m_d_l) || is_zero(d_r_m_s)) {
			r2 = math.Max(math.Abs(d_r_m_d_l/d_r_m_s), math.Abs(d_r_m_d_l/s_m_d_l))
		} else if preferShapePreservationOverSmoothness {
			r2 = maximum_rational_cubic_control_parameter_value
		}
	} else if monotonic && preferShapePreservationOverSmoothness {
		r2 = maximum_rational_cubic_control_parameter_value
	}
	return math.Max(minimum_rational_cubic_control_parameter_value, math.Max(r1, r2))
}

func convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_l float64, preferShapePreservationOverSmoothness bool) float64 {
	var r float64 = rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_l)
	var r_min float64 = minimum_rational_cubic_control_parameter(d_l, d_r, (y_r-y_l)/(x_r-x_l), preferShapePreservationOverSmoothness)
	return math.Max(r, r_min)
}

func convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_r float64, preferShapePreservationOverSmoothness bool) float64 {
	var r float64 = rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_r)
	var r_min float64 = minimum_rational_cubic_control_parameter(d_l, d_r, (y_r-y_l)/(x_r-x_l), preferShapePreservationOverSmoothness)
	return math.Max(r, r_min)
}

func compute_f_lower_map_and_first_two_derivatives(x, s float64) (f, fp, fpp float64) {
	var ax float64 = math.Abs(x)
	var z float64 = SQRT_ONE_OVER_THREE * ax / s
	var y float64 = z * z
	var s2 float64 = s * s
	var Phi float64 = norm_cdf(-z)
	var phi float64 = norm_pdf(z)

	fpp = PI_OVER_SIX * y / (s2 * s) * Phi * (8*SQRT_THREE*s*ax + (3*s2*(s2-8)-8*x*x)*Phi/phi) * math.Exp(2*y+0.25*s2)

	if is_below_horizon(s) {
		fp = 1
		f = 0
	} else {
		var Phi2 float64 = Phi * Phi
		fp = TWO_PI * y * Phi2 * math.Exp(y+0.125*s*s)
		if is_below_horizon(x) {
			f = 0
		} else {
			f = TWO_PI_OVER_SQRT_TWENTY_SEVEN * ax * (Phi2 * Phi)
		}
	}
	return
}

func inverse_f_lower_map(x, f float64) float64 {
	if is_below_horizon(f) {
		return 0
	}
	return math.Abs(x / (SQRT_THREE * inverse_norm_cdf(math.Pow(f/(TWO_PI_OVER_SQRT_TWENTY_SEVEN*math.Abs(x)), 1./3.))))
}

func compute_f_upper_map_and_first_two_derivatives(x, s float64) (f, fp, fpp float64) {
	f = norm_cdf(-0.5 * s)
	if is_below_horizon(x) {
		fp = -0.5
		fpp = 0
	} else {
		var w float64 = square(x / s)
		fp = -0.5 * math.Exp(0.5*w)
		fpp = SQRT_PI_OVER_TWO * math.Exp(w+0.125*s*s) * w / s
	}
	return
}

func inverse_f_upper_map(f float64) float64 {
	return -2. * inverse_norm_cdf(f)
}

func unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(beta, x, q float64, N int) float64 {
	if q*x > 0 {
		beta = math.Abs(math.Max(beta-normalised_intrinsic(x, q), 0.0))
		q = -q
	}
	// Map puts to calls
	if q < 0 {
		x = -x
		q = -q
	}
	if beta <= 0 { // For negative or zero prices we return 0.
		return 0
	}
	if beta < DENORMALIZATION_CUTOFF { // For positive but denormalized (a.k.a. 'subnormal') prices, we return 0 since it would be impossible to converge to full machine accuracy anyway.
		return 0
	}
	var b_max float64 = math.Exp(0.5 * x)
	if beta >= b_max {
		return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM
	}
	var iterations int = 0
	var direction_reversal_count int = 0
	var f float64 = -math.MaxFloat64
	var s float64 = -math.MaxFloat64
	var ds float64 = s
	var ds_previous float64 = 0
	var s_left float64 = math.SmallestNonzeroFloat64
	var s_right float64 = math.MaxFloat64
	var s_c float64 = math.Sqrt(math.Abs(2 * x))
	var b_c float64 = normalised_black_call(x, s_c)
	var v_c float64 = normalised_vega(x, s_c)

	if beta < b_c {
		var s_l float64 = s_c - b_c/v_c
		var b_l float64 = normalised_black_call(x, s_l)
		if beta < b_l {
			var f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2 float64
			f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2 = compute_f_lower_map_and_first_two_derivatives(x, s_l)
			var r_ll float64 = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(0., b_l, 0., f_lower_map_l, 1., d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2, true)
			f = rational_cubic_interpolation(beta, 0., b_l, 0., f_lower_map_l, 1., d_f_lower_map_l_d_beta, r_ll)
			if !(f > 0) {
				var t float64 = beta / b_l
				f = (f_lower_map_l*t + b_l*(1-t)) * t
			}
			s = inverse_f_lower_map(x, f)
			s_right = s_l
			for ; iterations < N && math.Abs(ds) > DBL_EPSILON*s; iterations++ {
				if ds*ds_previous < 0 {
					direction_reversal_count++
				}
				if iterations > 0 && (3 == direction_reversal_count || !(s > s_left && s < s_right)) {
					s = 0.5 * (s_left + s_right)
					if s_right-s_left <= DBL_EPSILON*s {
						break
					}
					direction_reversal_count = 0
					ds = 0
				}
				ds_previous = ds
				var b float64 = normalised_black_call(x, s)
				var bp float64 = normalised_vega(x, s)
				if b > beta && s < s_right {
					s_right = s
				} else if b < beta && s > s_left {
					s_left = s
				}
				if b <= 0 || bp <= 0 {
					ds = 0.5*(s_left+s_right) - s
				} else {
					var ln_b float64 = math.Log(b)
					var ln_beta float64 = math.Log(beta)
					var bpob float64 = bp / b
					var h float64 = x / s
					var b_halley float64 = h*h/s - s/4
					var newton float64 = (ln_beta - ln_b) * ln_b / ln_beta / bpob
					var halley float64 = b_halley - bpob*(1+2/ln_b)

					var b_hh3 float64 = b_halley*b_halley - 3*square(h/s) - 0.25
					var hh3 float64 = b_hh3 + 2*square(bpob)*(1+3/ln_b*(1+1/ln_b)) - 3*b_halley*bpob*(1+2/ln_b)
					ds = newton * householder_factor(newton, halley, hh3)
				}
				ds = math.Max(-0.5*s, ds)
				s += ds
			}
			return s
		} else {
			var v_l float64 = normalised_vega(x, s_l)
			var r_lm float64 = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(b_l, b_c, s_l, s_c, 1/v_l, 1/v_c, 0.0, false)
			s = rational_cubic_interpolation(beta, b_l, b_c, s_l, s_c, 1/v_l, 1/v_c, r_lm)
			s_left = s_l
			s_right = s_c
		}
	} else {
		var s_h float64
		if v_c > math.SmallestNonzeroFloat64 {
			s_h = s_c + (b_max-b_c)/v_c
		} else {
			s_h = s_c
		}
		var b_h float64 = normalised_black_call(x, s_h)
		if beta <= b_h {
			var v_h float64 = normalised_vega(x, s_h)
			var r_hm float64 = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_c, b_h, s_c, s_h, 1/v_c, 1/v_h, 0.0, false)
			s = rational_cubic_interpolation(beta, b_c, b_h, s_c, s_h, 1/v_c, 1/v_h, r_hm)
			s_left = s_c
			s_right = s_h
		} else {
			var f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2 float64
			f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2 = compute_f_upper_map_and_first_two_derivatives(x, s_h)
			if d2_f_upper_map_h_d_beta2 > -SQRT_DBL_MAX && d2_f_upper_map_h_d_beta2 < SQRT_DBL_MAX {
				var r_hh float64 = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_h, b_max, f_upper_map_h, 0., d_f_upper_map_h_d_beta, -0.5, d2_f_upper_map_h_d_beta2, true)
				f = rational_cubic_interpolation(beta, b_h, b_max, f_upper_map_h, 0., d_f_upper_map_h_d_beta, -0.5, r_hh)
			}
			if f <= 0 {
				var h float64 = b_max - b_h
				var t float64 = (beta - b_h) / h
				f = (f_upper_map_h*(1-t) + 0.5*h*t) * (1 - t)
			}
			s = inverse_f_upper_map(f)
			s_left = s_h
			if beta > 0.5*b_max {
				for ; iterations < N && math.Abs(ds) > DBL_EPSILON*s; iterations++ {
					if ds*ds_previous < 0 {
						direction_reversal_count++
					}
					if iterations > 0 && (3 == direction_reversal_count || !(s > s_left && s < s_right)) {
						s = 0.5 * (s_left + s_right)
						if s_right-s_left <= DBL_EPSILON*s {
							break
						}
						direction_reversal_count = 0
						ds = 0
					}
					ds_previous = ds
					var b float64 = normalised_black_call(x, s)
					var bp float64 = normalised_vega(x, s)
					if b > beta && s < s_right {
						s_right = s
					} else if b < beta && s > s_left {
						s_left = s
					}
					if b >= b_max || bp <= math.SmallestNonzeroFloat64 {
						ds = 0.5*(s_left+s_right) - s
					} else {
						var b_max_minus_b float64 = b_max - b
						var g float64 = math.Log((b_max - beta) / b_max_minus_b)
						var gp float64 = bp / b_max_minus_b
						var b_halley float64 = square(x/s)/s - s/4
						var b_hh3 float64 = b_halley*b_halley - 3*square(x/(s*s)) - 0.25
						var newton float64 = -g / gp
						var halley float64 = b_halley + gp
						var hh3 float64 = b_hh3 + gp*(2*gp+3*b_halley)
						ds = newton * householder_factor(newton, halley, hh3)
					}
					ds = math.Max(-0.5*s, ds)
					s += ds
				}
				return s
			}
		}
	}
	for ; iterations < N && math.Abs(ds) > DBL_EPSILON*s; iterations++ {
		if ds*ds_previous < 0 {
			direction_reversal_count++
		}
		if iterations > 0 && (3 == direction_reversal_count || !(s > s_left && s < s_right)) {
			s = 0.5 * (s_left + s_right)
			if s_right-s_left <= DBL_EPSILON*s {
				break
			}
			direction_reversal_count = 0
			ds = 0
		}
		ds_previous = ds
		var b float64 = normalised_black_call(x, s)
		var bp float64 = normalised_vega(x, s)
		if b > beta && s < s_right {
			s_right = s
		} else if b < beta && s > s_left {
			s_left = s
		}
		var newton float64 = (beta - b) / bp
		var halley float64 = square(x/s)/s - s/4
		var hh3 float64 = halley*halley - 3*square(x/(s*s)) - 0.25
		ds = math.Max(-0.5*s, newton*householder_factor(newton, halley, hh3))
		s += ds
	}
	return s
}

func implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(price, F, K, T, q float64, N int) float64 {
	var intrinsic float64
	if q < 0 {
		intrinsic = math.Abs(math.Max(K-F, 0.0))
	} else {
		intrinsic = math.Abs(math.Max(F-K, 0.0))
	}
	if price < intrinsic {
		return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC
	}
	var max_price float64
	if q < 0 {
		max_price = K
	} else {
		max_price = F
	}
	if price >= max_price {
		return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM
	}
	var x float64 = math.Log(F / K)
	// Map in-the-money to out-of-the-money
	if q*x > 0 {
		price = math.Abs(math.Max(price-intrinsic, 0.0))
		q = -q
	}
	return unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(price/(math.Sqrt(F)*math.Sqrt(K)), x, q, N) / math.Sqrt(T)
}

func implied_volatility_from_a_transformed_rational_guess(price, F, K, T, q float64) float64 {
	return implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(price, F, K, T, q, implied_volatility_maximum_iterations)
}

func normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(beta, x, q float64, N int) float64 {
	// Map in-the-money to out-of-the-money
	if q*x > 0 {
		beta -= normalised_intrinsic(x, q)
		q = -q
	}
	if beta < 0 {
		return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC
	}
	return unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(beta, x, q, N)
}

func Normalised_implied_volatility_from_a_transformed_rational_guess(beta, x, q float64) float64 {
	return normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(beta, x, q, implied_volatility_maximum_iterations)
}
