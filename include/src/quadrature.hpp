#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

#include "hypergeometric.hpp"

/* Gauss-Laguerre and Gauss-Jacobi quadrature rules using Golub-Welsch */

/* ----------- BASE CLASSES ------------ */

struct quadrature_node_t
{
  /* Container for quadrature nodes with positive weights */
  constexpr static double inf = std::numeric_limits<double>::infinity();
  double abscissa, logweight;
  quadrature_node_t (void) : abscissa(0.0), logweight(-inf) {};
  quadrature_node_t (double x, double w) : abscissa(x), logweight(w) {};
};

struct gauss_quadrature_t
{
  /* Machinery for the Golub-Welsch algorithm */

  protected:

    void diagonalize_matrix(const std::size_t n, std::vector<double>& d, std::vector<double>& e, std::vector<double>& z, const std::size_t max_iter = 30)
    {
      /* Diagonalize a symmetric, tridiagonal matrix. 
       * The diagonal `d` and offdiagonal `e` are overwritten.
       * The vector `z` is multiplied in-place by the eigenvectors.
       *
       * Used and modified under the GNU LGPL license. 
       * Source: https://people.math.sc.edu/burkardt/cpp_src/jacobi_rule/jacobi_rule.cpp
       * Modified from C++ translation of FORTRAN77 by John Burkardt
       * Original FORTRAN77 by Slyvan Elhay and Jaroslav Kautsky
       */
      double b, c, f, g, p, r, s;
      const double prec = 2.220446049250313E-016;
      int i, ii, j, k, l, m, mml;
      if (n == 1) return;
      e[n-1] = 0.0;
      for (l = 1; l <= n; ++l) {
        j = 0;
        do {
            for (m = l; m <= n; ++m) {
              if (m == n) break;
              if (fabs(e[m-1]) <= prec*(fabs(d[m-1]) + fabs(d[m]))) break;
            }
            p = d[l-1];
            if (m == l) break;
            if (j >= max_iter) {
              /* failed to converge */
              std::fill(d.begin(), d.end(), std::numeric_limits<double>::quiet_NaN());
              std::fill(e.begin(), e.end(), std::numeric_limits<double>::quiet_NaN());
              std::fill(z.begin(), z.end(), std::numeric_limits<double>::quiet_NaN());
              return;
            }
            j += 1;
            g = (d[l] - p) / (2.0*e[l-1]);
            r =  std::sqrt(std::pow(g, 2) + 1.0);
            g = d[m-1] - p + e[l-1] / (g + fabs(r)*std::copysign(1.0, g));
            s = 1.0;
            c = 1.0;
            p = 0.0;
            mml = m - l;
            for (ii = 1; ii <= mml; ++ii) {
              i = m - ii;
              f = s*e[i-1];
              b = c*e[i-1];
              if (fabs(g) <= fabs(f)) {
                c = g / f;
                r =  std::sqrt(std::pow(c, 2) + 1.0);
                e[i] = f*r;
                s = 1.0 / r;
                c = c*s;
              } else {
                s = f / g;
                r =  std::sqrt(std::pow(s, 2) + 1.0);
                e[i] = g*r;
                c = 1.0 / r;
                s = s*c;
              }
              g = d[i] - p;
              r = (d[i-1] - g)*s + 2.0*c*b;
              p = s*r;
              d[i] = g + p;
              g = c*r - b;
              f = z[i];
              z[i] = s*z[i-1] + c*f;
              z[i-1] = c*z[i-1] - s*f;
            }
            d[l-1] = d[l-1] - p;
            e[l-1] = g;
            e[m-1] = 0.0;
        } while (true);
      } 
      /* sorting */
      for (ii = 2; ii <= m; ++ii) {
        i = ii - 1;
        k = i;
        p = d[i-1];
        for (j = ii; j <= n; ++j) {
          if (d[j-1] < p) {
             k = j;
             p = d[j-1];
          }
        }
        if (k != i) {
          d[k-1] = d[i-1];
          d[i-1] = p;
          p = z[i-1];
          z[i-1] = z[k-1];
          z[k-1] = p;
        }
      }
      return;
    }

    void golub_welsch (const std::size_t order, std::vector<double>& aj, std::vector<double>& bj, std::vector<double>& wj)
    {
      /* Calculate knots and weights using Golub-Welsch algorithm */
      assert (wj.size() == order);
      assert (aj.size() == order);
      assert (bj.size() == order);
      std::fill(wj.begin(), wj.end(), 0.0); wj[0] = 1.0; // normalized weights
      diagonalize_matrix(order, aj, bj, wj);
      std::for_each(wj.begin(), wj.end(), [](double& x){x = 2*log(fabs(x));});
    }


    void recurrence_from_moments(const std::vector<double>& moments, std::vector<double>& aj, std::vector<double>& bj)
    {
      /* Upper-triangular Cholesky factorization of the Hankel matrix of moments  
         See Golub & Welsch 1969 for motivation */
      assert (moments.size() % 2 == 0);
      const std::size_t n = moments.size() / 2;
      assert (aj.size() == n - 1);
      assert (bj.size() == n - 1);
      std::vector<double> upper (n*n, 0.0); // flattened in row-major
      double s;
      for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < i + 1; ++j) {
          s = 0.0;
          if (j == i) {
            for (std::size_t k = 0; k < j; ++k) s += std::pow(upper[k*n+j], 2);
            upper[j*n+j] = sqrt(moments[j+j] - s);
          } else {
            for (std::size_t k = 0; k < j; ++k) s += upper[k*n+i] * upper[k*n+j];
            if (upper[j*n+j] > 0.0) upper[j*n+i] = (moments[i+j] - s) / upper[j*n+j];
          }
        }
      }

      /* Recursion coefficients for orthogonal polynomial 
         Eq. 4.3 in Golub & Welsch 1969 */
			aj.resize(n - 1);
    	bj.resize(n - 1);
      for (std::size_t j = 0; j < n - 1; ++j) {
        aj[j] = upper[n*j+j+1]/upper[n*j+j]; //R[j,j+1]/R[j,j]
        if (j > 0.0) aj[j] -= upper[n*(j-1)+j]/upper[n*(j-1)+j-1]; //R[j-1,j]/R[j-1, j-1]
        bj[j] = upper[n*(j+1)+j+1]/upper[n*j+j]; //R[j+1,j+1]/R[j,j]
      } 
    }
};


/* ----------- SPECIFIC RULES ------------ */

struct gauss_laguerre_t : protected gauss_quadrature_t
{
  /* Gauss-Laguerre rule for `f(x) x^alpha exp(-x)` on [0, infty) */

  private:
  double recurrence (const std::size_t m, const double alpha, std::vector<double>& aj, std::vector<double>& bj) 
  {
    /* Calculate Jacobi matrix */ 
    assert (bj.size() == m);
    assert (aj.size() == m);
    double logconst = std::lgamma(alpha + 1.0);
    for (std::size_t i = 1; i <= m; ++i) {
        aj[i-1] = 2.0*i - 1.0 + alpha;
        bj[i-1] = std::sqrt(i*(i + alpha));
    }
    return logconst;
  }

  public:
  std::vector<quadrature_node_t> grid;

  gauss_laguerre_t (const std::size_t order, const double alpha)
  {
	  std::vector<double> offdiag(order), knots(order), logweights(order);
    double logconst = recurrence(order, alpha, knots, offdiag);
    golub_welsch(order, knots, offdiag, logweights);
    std::for_each(logweights.begin(), logweights.end(), 
        [&logconst](double& x){ x += logconst; }
    );
    /* fill grid */
    grid.reserve(order);
    for (std::size_t i = 0; i < order; ++i) grid.emplace_back(knots[i], logweights[i]);
  }
};


struct gauss_jacobi_t : protected gauss_quadrature_t
{
  /* Gauss-Jacobi quadrature for `f(x) (1 - x)^(alpha - beta) x^(beta - 1)` on [0, 1] */

  private:
  double recurrence (const std::size_t m, const double alpha, const double beta, std::vector<double>& aj, std::vector<double>& bj) 
  {
    /* Calculate Jacobi matrix */ 
    assert (bj.size() == m);
    assert (aj.size() == m);
    double ab = alpha + beta, abi = 2.0 + ab, a2b2 = std::pow(beta, 2) - std::pow(alpha, 2);
    double logconst = log(2) * (ab + 1.0) + std::lgamma(alpha + 1.0) + 
      std::lgamma(beta + 1.0) - std::lgamma(abi);
    aj[0] = (beta - alpha) / abi;
    bj[0] = std::sqrt(4*(1.0 + alpha)*(1.0 + beta)/((abi + 1.0)*std::pow(abi, 2)));
    for (std::size_t i = 2; i <= m; ++i)
    {
      abi = 2.0*i + ab;
      aj[i-1] = a2b2 / ((abi - 2.0)*abi);
      abi = abi * abi;
      bj[i-1] = std::sqrt(4.0*i*(i + alpha)*(i + beta)*(i + ab) / ((abi - 1.0) * abi));
    }
    return logconst;
  }

  public:
  std::vector<quadrature_node_t> grid;

  gauss_jacobi_t (const std::size_t order, const double alpha, const double beta)
  {
	  std::vector<double> offdiag(order), knots(order), logweights(order);
    auto logconst = recurrence(order, alpha - beta, beta - 1.0, knots, offdiag);
    golub_welsch(order, knots, offdiag, logweights);
    std::for_each(logweights.begin(), logweights.end(), 
        [&logconst](double& x){ x += logconst; }
    );
    /* rescale to [0, 1] */
    std::for_each(knots.begin(), knots.end(), [](double& x){ x = (x + 1.0) / 2.0; });
    double scale = log(2.0) * alpha;
    std::for_each(logweights.begin(), logweights.end(), [&scale](double& x){ x -= scale; });
    /* fill grid */
    grid.reserve(order);
    for (std::size_t i = 0; i < order; ++i) grid.emplace_back(knots[i], logweights[i]);
  }
};


struct gauss_2f1_t : protected gauss_quadrature_t
{
  /* Gauss quadrature for `f(x) x^alpha (1 - zeta x)^beta (alpha + 1)` on [0, 1] */

  private:
  double recurrence(const std::size_t order, const double alpha, const double beta, const double zeta, std::vector<double>& aj, std::vector<double>& bj) {
    /* n-th moment is 2F1(-beta, alpha + n + 1; alpha + n + 2; zeta)` */
    std::vector<double> moments(2*order + 2);
    for (std::size_t n = 0; n < 2*order + 2; ++n) {
      moments[n] = std::pow(1.0/zeta, alpha + n + 1) *
        betainc_t::lentz(alpha + n + 1, beta + 1, zeta); // rewrite 2F1 as incomplete beta
    }
    // TODO: this could be done using a recurrence relation:
    // double I, B;
    // B = exp(a * log(x) + b * log(1.0 - x) - log(a) - std::lgamma(a) - std::lgamma(b) + std::lgamma(a + b));
    // I = betainc_t::lentz(a, b, x);
    // for (std::size_t i = 0; i < order; ++i) {
    //   I -= B;
    //   B *= (a + b) / (a + 1.0) * x;
    // }
    // (however, `I` is regularized)

    /* recurrence coefficients from moments */
    recurrence_from_moments(moments, aj, bj);

    return moments[0];
  }

  public:
  std::vector<quadrature_node_t> grid;

  gauss_2f1_t (const std::size_t order, const double alpha, const double beta, const double zeta)
  {
	  std::vector<double> offdiag(order), knots(order), logweights(order);
    auto logconst = recurrence(order, alpha, beta, zeta, knots, offdiag);
    golub_welsch(order, knots, offdiag, logweights);
    std::for_each(logweights.begin(), logweights.end(), 
        [&logconst](double& x){ x += logconst; }
    );
    /* fill grid */
    grid.reserve(order);
    for (std::size_t i = 0; i < order; ++i) grid.emplace_back(knots[i], logweights[i]);
  }
};


/* ----------- STRATIFIED RULES ------------ */

//template <class T>
//struct antigauss_quadrature : protected gauss_quadrature_t
//{
//  // take an existing rule, use it's class matrix to make a new class matrix
//};


struct nested_gauss_laguerre_t : protected gauss_quadrature_t
{
  /* Stratified Gauss-Laguerre rule for `f(x) x^alpha exp(-x) / gamma(alpha + 1)` on [0, infty) 
     The stratification method in XXXX gives rise to two nested quadrature rules */

  private:
  double recurrence (const std::size_t m, const double alpha, const double beta, std::vector<double>& aj, std::vector<double>& bj) 
  {
    /* Calculate Jacobi matrix */ 
    std::size_t M = 2*m + 1;
    assert (bj.size() == M);
    assert (aj.size() == M);
    double logconst = std::lgamma(alpha + 1.0);
    // (m + 1)-order rule
    for (std::size_t i = 1; i <= m + 1; ++i) {
        aj[i-1] = 2.0*i - 1.0 + alpha;
        bj[i-1] = std::sqrt(i*(i + alpha));
    }
    // m-order rule in reverse
    std::size_t k = m + 1;
    for (std::size_t i = m; i >= 1; --i) {
        aj[k] = aj[i-1];
        if (i != 1) bj[k] = bj[i-2];
        k++;
    }
    return logconst;
  }

  public:
  std::vector<quadrature_node_t> grid;

  nested_gauss_laguerre_t (const std::size_t order, const double alpha)
  {
	  std::vector<double> offdiag(2*order + 1), knots(2*order + 1), logweights(2*order + 1);
    double logconst = recurrence(order, alpha, 1.0, knots, offdiag);
    golub_welsch(2*order + 1, knots, offdiag, logweights);
    std::for_each(logweights.begin(), logweights.end(), 
        [&logconst](double& x){ x += logconst; }
    );
    /* fill grid */
    grid.reserve(2*order + 1);
    for (std::size_t i = 0; i < 2*order + 1; ++i) grid.emplace_back(knots[i], logweights[i]);
  }
};

struct extended_gauss_laguerre_t : protected gauss_quadrature_t
{
  /* Stratified Gauss-Laguerre rule for `f(x) x^alpha exp(-x) / gamma(alpha + 1)` on [0, infty) 
     The stratification method in XXXX gives rise to two nested quadrature rules */

  private:
  double recurrence (const std::size_t m, const double alpha, const double beta, std::vector<double>& aj, std::vector<double>& bj) 
  {
    /* Calculate Jacobi matrix */ 
    assert (bj.size() == m + 1);
    assert (aj.size() == m + 1);
    double logconst = std::lgamma(alpha + 1.0);
    // (m + 1)-order rule
    for (std::size_t i = 1; i <= m + 1; ++i) {
        aj[i-1] = 2.0*i - 1.0 + alpha;
        bj[i-1] = std::sqrt(i*(i + alpha));
    }
    bj[m-1] = std::sqrt(bj[m-1]*bj[m-1] + bj[m]*bj[m]);
    bj[m] = 0.0;
    return logconst;
  }

  public:
  std::vector<quadrature_node_t> grid;

  extended_gauss_laguerre_t (const std::size_t order, const double alpha)
  {
	  std::vector<double> offdiag(order + 1), knots(order + 1), logweights(order + 1);
    double logconst = recurrence(order, alpha, 1.0, knots, offdiag);
    golub_welsch(order + 1, knots, offdiag, logweights);
    std::for_each(logweights.begin(), logweights.end(), 
        [&logconst](double& x){ x += logconst; }
    );
    /* fill grid */
    grid.reserve(order + 1);
    for (std::size_t i = 0; i < order + 1; ++i) grid.emplace_back(knots[i], logweights[i]);
  }
};


#endif
