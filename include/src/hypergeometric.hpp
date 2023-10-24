#ifndef HYPERGEO_H
#define HYPERGEO_H

#include <cmath>
#include <cassert>
#include <iostream>

#define HYP1F1_RELTOL 1e-8
#define HYP1F1_MAXITT 1000000
#define BETAIN_RELTOL 1e-8
#define BETAIN_MAXITT 100

#define TINY 1e-30
#define PI 3.14159265358979323846
#define EULER 2.71828182845904523536


/* ------------ POLYGAMMA -------------- */

double digamma(double x)
{
  /* Used and modified under GNU LGPL license 
     Original source: https://people.math.sc.edu/Burkardt/cpp_src/asa103
     Copyright (c) 2008 John Burkardt */
  if (x <= 0.0) return digamma(1.0 - x) - PI / std::tan(PI * x); // analytic continuation
  if (x <= 1e-5) return -EULER - 1.0 / x;
  if (x < 8.5) return digamma(1.0 + x) - 1.0 / x;
  double xpm2 = 1.0 / std::pow(x, 2);
  return log(x) - 0.5 / x // de Moivre's expansion
      - 0.083333333333333333 * xpm2
      + 0.008333333333333333 * std::pow(xpm2, 2)
      - 0.003968253968253968 * std::pow(xpm2, 3)
      + 0.004166666666666667 * std::pow(xpm2, 4)
      - 0.007575757575757576 * std::pow(xpm2, 5)
      + 0.021092796092796094 * std::pow(xpm2, 6);
}


double trigamma(double x)
{
  /* Used and modified under GNU LGPL license 
     Original source: https://people.math.sc.edu/Burkardt/cpp_src/asa121
     Copyright (c) 2008 John Burkardt */
  if (x <= 0.0) return -trigamma(1.0 - x) + std::pow(PI / std::sin(PI * x), 2); // analytic continuation
  if (x <= 1e-4) return 1.0 / std::pow(x, 2);
  if (x < 5) return trigamma(1.0 + x) + 1.0 / std::pow(x, 2);
  double xpm1 = 1.0 / x, xpm2 = 1.0 / std::pow(x, 2);
  return xpm1 * ( // de Moivre's expansion
      1.000000000000000000
      + 0.500000000000000000 * xpm1
      + 0.166666666666666667 * xpm2
      - 0.033333333333333333 * std::pow(xpm2, 2)
      + 0.023809523809523808 * std::pow(xpm2, 3)
      - 0.033333333333333333 * std::pow(xpm2, 4)
      + 0.075757575757575756 * std::pow(xpm2, 5)
      - 0.253113553113553102 * std::pow(xpm2, 6)
      + 1.166666666666666741 * std::pow(xpm2, 7)
  );
}


/* -------------- INCOMPLETE BETA ---------------- */

struct betainc_t 
{
  static
  double lentz(double a, double b, double x)
  {
    /* Used and modified under zlib license 
       Original source: https://github.com/codeplea/incbeta
       Copyright (c) 2016, 2017 Lewis Van Winkle */
  
    if (x < 0.0 || x > 1.0) return std::numeric_limits<double>::quiet_NaN();
    if (x > (a + 1.0) / (a + b + 2.0)) {
        return 1.0 - betainc_t::lentz(b, a, 1.0 - x);
    }
    
    /* find the first part before the continued fraction */
    const double lbeta_ab = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
    const double front = exp(log(x) * a + log(1.0 - x) * b - lbeta_ab) / a;
    
    /* use Lentz's algorithm to evaluate the continued fraction */
    double f = 1.0, c = 1.0, d = 0.0;
    
    int i, m;
    for (i = 0; i < BETAIN_MAXITT; ++i) {
      m = i / 2;
      
      double numerator;
      if (i == 0) {
        numerator = 1.0;
      } else if (i % 2 == 0) {
        numerator = (m * (b - m) * x)/((a + 2.0 * m - 1.0) * (a + 2.0 * m));
      } else {
        numerator = -((a + m) * (a + b + m) * x)/((a + 2.0 * m) * (a + 2.0 * m + 1));
      }
      
      /* do an iteration of Lentz's algorithm */
      d = 1.0 + numerator * d;
      if (fabs(d) < TINY) d = TINY;
      d = 1.0 / d;
      c = 1.0 + numerator / c;
      if (fabs(c) < TINY) c = TINY;
      
      const double cd = c * d;
      f *= cd;
      
      /* check for stop */
      if (fabs(1.0 - cd) < BETAIN_RELTOL) return front * (f - 1.0);
    }
    
    return std::numeric_limits<double>::quiet_NaN(); // did not converge
  }
};


/* --------------- CONFLUENT HYPERGEOMETRIC --------------- */

struct hyp1f1_t
{
  //static
  //double is_valid_1f1(f1, f2, a, b, x)
  //{
  //  /* ... */

  //}

  //static
  //double poinare_expansion (double a, double b, double x, double& da, double& db, double& dx, double& d2x)
  //{
  //  /* ... */
  //}

  static
  double taylor_series (double a, double b, double x, double& da, double& db, double& dx, double& d2x) 
  {
    /* Evaluate via a series expansion around the origin. Returns log function
       value, and derivatives scaled by inverse of function value */
    assert (a > 0 && b > 0);
    double val;
    if (x == 0.0) {
      val = 0.0;
      da = 0.0;
      db = 0.0;
      dx = a / b;
      d2x = dx * (a + 1) / (b + 1);
      return val;
    } else if (x < 0.0) {
      /* Kummer's transformation, 1F1(a, b, x) = exp(x) 1F1(b - a, b, -x) */
      val = x + hyp1f1_t::taylor_series(b - a, b, -x, da, db, dx, d2x);
      db += da;
      da *= -1.0;
      d2x += 1.0 - 2*dx;
      dx = 1.0 - dx;
      return val;
    }
    //else if (x > 100.0) {
    //  /* Poincare-type asymptotic expansion */
    //}

    double bi, ai, xk, bk, ak, xd, bd, ad, offset, weight, norm;
    std::size_t k = 1;
    const double ltol = log(HYP1F1_RELTOL);
    ai = a, bi = b; // state (increment) of args
    xk = bk = ak = 0.0; // multiplicative update for args
    xd = bd = ad = 0.0; // additive update for args
    offset = 0.0;
    da = db = dx = d2x = 0.0;
    val = 1.0;
    do {
    	ak += log(ai + k - 1);
    	bk += log(bi + k - 1);
      ad += 1.0 / (ai + k - 1);
      bd += 1.0 / (bi + k - 1);
      xk += log(x) - log(k);
      xd += 1.0 / x;
      weight = ak - bk + xk;
      norm = exp(weight - offset);
      if (weight < offset) {
        val += norm;
        da += ad * norm;
        db += bd * norm;
        dx += xd * norm;
        d2x += xd * (xd - 1.0 / x) * norm;
      } else {
        val = 1.0 + val / norm;
        da = ad + da / norm;
        db = bd + db / norm;
        dx = xd + dx / norm;
        d2x = xd * (xd - 1.0 / x) + d2x / norm;
        offset = weight;
      }
      if (!std::isfinite(val) || k > HYP1F1_MAXITT) break;
      //if (weight < ltol + log(val) && hyp1f1_t::is_valid_1f1(...)) {
      if (weight < ltol + log(val)) {
        da /= val;
        db /= -val;
        dx /= val;
        d2x /= val;
        val = log(val) + offset;
        return val;
      }
      k++;
    } while (true);
    /* failed to converge */
    val = da = db = dx = d2x = std::numeric_limits<double>::quiet_NaN();
    return val;
  }
};

#endif
