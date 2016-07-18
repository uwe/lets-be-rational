package lbr

import (
	"math"
	"testing"
)

type Test struct {
	CallPut      float64
	Strike       float64
	Maturity     float64
	RiskFreeRate float64
	DividendRate float64
	Underlying   float64
	Price        float64
	IV           float64
}

var tests = []Test{
	Test{1, 1220, 6.0 / 365, 0.0001, 0.0127, 1240.41, 22.1, 13.46},
	Test{-1, 1220, 6.0 / 365, 0.0001, 0.0127, 1240.41, 3.1, 16.13},
	Test{1, 1240, 6.0 / 365, 0.0001, 0.0127, 1240.41, 8.1, 12.65},
	Test{-1, 1240, 6.0 / 365, 0.0001, 0.0127, 1240.41, 9.4, 14.94},

	Test{1, 1220, 4.0 / 252, 0.0001, 0.0127, 1240.41, 22.1, 13.68},
	Test{-1, 1220, 4.0 / 252, 0.0001, 0.0127, 1240.41, 3.1, 16.42},
	Test{1, 1240, 4.0 / 252, 0.0001, 0.0127, 1240.41, 8.1, 12.87},
	Test{-1, 1240, 4.0 / 252, 0.0001, 0.0127, 1240.41, 9.4, 15.21},
}

func TestIV(t *testing.T) {
	for _, test := range tests {
		F := test.Underlying * math.Exp((test.RiskFreeRate-test.DividendRate)*test.Maturity)
		iv := 100 * implied_volatility_from_a_transformed_rational_guess(
			test.Price,
			F,
			test.Strike,
			test.Maturity,
			test.CallPut,
		)
		if math.Abs(iv-test.IV) > 0.01 {
			t.Error("Expected", test.IV, "Got", iv)
		}
	}
}
