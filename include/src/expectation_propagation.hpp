#ifndef EXPECTPROP_H
#define EXPECTPROP_H

/* Use expectation propagation to infer node ages in marginal trees */

#include <iostream>
#include <random>
#include <cassert>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

#include "anc.hpp"
#include "data.hpp"

#define OPTIM_MAXITT 100
#define OPTIM_RELTOL 1e-8
#define IS_MINGAP 1e-12
#define HYP1F1_RELTOL 1e-8
#define HYP1F1_MAXITT 1000000
#define BETAIN_RELTOL 1e-8
#define BETAIN_MAXITT 100
#define TINY 1e-30
#define PI 3.14159265358979323846
#define EULER 2.71828182845904523536


// ------------- SPECIAL FUNCTIONS -------------- //

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

double betainc(double a, double b, double x)
{
  /* Used and modified under zlib license 
     Original source: https://github.com/codeplea/incbeta
     Copyright (c) 2016, 2017 Lewis Van Winkle */

  if (x < 0.0 || x > 1.0) return std::numeric_limits<double>::quiet_NaN();
  if (x > (a + 1.0) / (a + b + 2.0)) {
      return 1.0 - betainc(b, a, 1.0 - x);
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


// ---------- GAUSSIAN QUADRATURE -------------- //

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
    assert (alpha > -1.0);
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
  /* Gauss-Jacobi quadrature for `f(x) (1 - x)^alpha x^beta` on [0, 1] */

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
    assert (beta > -1.0);
    assert (alpha > -1.0);
	  std::vector<double> offdiag(order), knots(order), logweights(order);
    auto logconst = recurrence(order, alpha, beta, knots, offdiag);
    golub_welsch(order, knots, offdiag, logweights);
    std::for_each(logweights.begin(), logweights.end(), 
        [&logconst](double& x){ x += logconst; }
    );
    /* rescale to [0, 1] */
    std::for_each(knots.begin(), knots.end(), [](double& x){ x = (x + 1.0) / 2.0; });
    double scale = log(2.0) * (alpha + beta);
    std::for_each(logweights.begin(), logweights.end(), [&scale](double& x){ x -= scale; });
    /* fill grid */
    grid.reserve(order);
    for (std::size_t i = 0; i < order; ++i) grid.emplace_back(knots[i], logweights[i]);
  }
};


// --------------- GAMMA DISTRIBUTION ------------------- //

struct gamma_distr_t
{
  /* Gamma distribution in natural parameterization */

  double alpha, beta;

  gamma_distr_t(void)
    : alpha(0.0)
    , beta(0.0)
  {}

  gamma_distr_t(const double& alpha, const double& beta) 
    : alpha(alpha)
    , beta(beta)
  {}

  /* alternative parameterizations */
  double shape(void) const {
    return alpha + 1;
  }

  double rate(void) const {
    return -beta;
  }

  double scale(void) const {
    return beta < 0 ? -1.0 / beta : 0.0;
  }

  double mean(void) const {
    return beta < 0 ? (alpha + 1) / -beta : 0.0;
  }

  double square(void) const {
    return beta < 0 ? variance() * (alpha + 2) : 0.0;
  }

  double logmean(void) const {
    return beta < 0 ? digamma(alpha + 1) - log(-beta) : std::numeric_limits<double>::infinity();
  }

  double variance(void) const {
    return beta < 0 ? (alpha + 1) / std::pow(beta, 2) : 0.0;
  }

  /* random number generator */
  std::gamma_distribution<double> random(void) const {
    return std::gamma_distribution<double>(shape(), scale());
  }

  /* log likelihood */
  double loglik(const double& x) const {
    return log(x) * alpha + x * beta;
  }

  /* instantiate from alternative parameterizations */
  static gamma_distr_t from_moments(const double& mn, const double& va)
  {
    return va > 0 && mn > 0 ? 
      gamma_distr_t(std::pow(mn, 2) / va - 1, -mn / va) :
      gamma_distr_t();
  }

  static gamma_distr_t from_canonical(const double& shape, const double& rate)
  {
    return rate > 0 && shape > 0 ? 
      gamma_distr_t(shape - 1, -rate) :
      gamma_distr_t();
  }

  static gamma_distr_t from_sufficient(const double& mn, const double& ln)
  {
    /* use Newton root finding to solve: E[log x] = psi(shape) - log(rate), E[x] = shape / rate */
    if (mn <= 0.0 || std::isinf(ln)) return gamma_distr_t();
    assert (log(mn) > ln);
    double shape = 0.5 / (log(mn) - ln);
    if (1.0 / shape < 1e-4) return gamma_distr_t::from_canonical(shape, shape / mn);
    double delta = std::numeric_limits<double>::infinity();
    for (std::size_t itt=0; itt < OPTIM_MAXITT; itt++)
    {
      if (fabs(delta) < shape * OPTIM_RELTOL) break;
      delta = digamma(shape) - log(shape) + log(mn) - ln;
      delta /= trigamma(shape) - 1 / shape;
      shape -= delta;
    }
    assert (shape > 0);
    return gamma_distr_t::from_canonical(shape, shape / mn);
  }

  /* operator overloads */
  void operator=(const gamma_distr_t& other)
  {
    alpha = other.alpha;
    beta = other.beta;
  }

  gamma_distr_t operator+(const gamma_distr_t& other) const  
  {
    return gamma_distr_t(alpha + other.alpha, beta + other.beta);
  }

  void operator+=(const gamma_distr_t& other)
  {
    alpha += other.alpha;
    beta += other.beta;
  }

  gamma_distr_t operator-(const gamma_distr_t& other) const
  {
    return gamma_distr_t(alpha - other.alpha, beta - other.beta);
  }

  void operator-=(const gamma_distr_t& other)
  {
    alpha -= other.alpha;
    beta -= other.beta;
  }

  friend std::ostream& operator<<(std::ostream& os, const gamma_distr_t& self)
  {
    os << "Gamma(" << self.shape() << ", " << self.rate() << ")";
    return os;
  }

  explicit operator bool() const
  {
    /* is integrable */
    return beta < 0.0 && alpha > -1.0 ? true : false;
  }
};


// ---------------- PRIORS ----------------- //

struct trio_t
{
  /* Container for parent-child trios */
  std::size_t parent, lchild, rchild;
  trio_t (std::size_t parent, std::size_t lchild, std::size_t rchild)
    : parent(parent), lchild(lchild), rchild(rchild) 
  {}
};


struct timescale_t
{
  /* Convert node ages between coalescent and generational scales under a piecewise-constant coalescence rate */

  const std::size_t dim;

  private:
    std::vector<double> gens_breaks, coal_breaks, gens_step, coal_step, rate;

    static
    double average_rate (const std::vector<double>& epoch, const std::vector<double>& coal) {
				/* expected pairwise coalescence rate under the piecewise constant history 
           (integrate PDF for pairwise coalescence times under time rescaling) */
        std::vector<double> duration (epoch.size() - 1);
        std::vector<double> weight (epoch.size(), 1.0);
        std::vector<double> integrand (epoch.size(), 1.0);
        for (std::size_t i = 1; i < epoch.size(); ++i) {
          duration[i-1] = exp(-(epoch[i] - epoch[i-1]) * coal[i-1]);
          weight[i] = weight[i-1] * duration[i-1];
          integrand[i-1] = (1.0 - duration[i-1]) / coal[i-1];
        }
        integrand.back() /= coal.back();
        double expected_rate = 
          std::inner_product(weight.begin(), weight.end(), integrand.begin(), 0.0);
        return expected_rate;
    }

  public:
    double pairwise_tmrca;

    timescale_t(const std::vector<double>& epoch_breaks, const std::vector<double>& epoch_rate)
      : dim (epoch_rate.size())
      , rate (epoch_rate)
      , gens_breaks (dim, 0.0)
      , coal_breaks (dim, 0.0)
      , gens_step (dim, 0.0)
      , coal_step (dim, 0.0)
    {
			assert (epoch_breaks.size() == dim);
      assert (epoch_breaks[0] == 0.0);
      for (std::size_t i = 1; i < dim; ++i) {
        assert (epoch_breaks[i] - epoch_breaks[i - 1] > 0);
        gens_breaks[i] = epoch_breaks[i];
        gens_step[i] = gens_step[i - 1] + gens_breaks[i] * (rate[i - 1] - rate[i]);
        coal_breaks[i] = gens_breaks[i] * rate[i] + gens_step[i];
        coal_step[i] = coal_step[i - 1] + coal_breaks[i] * (1./rate[i - 1] - 1./rate[i]);
      }
      pairwise_tmrca = average_rate(epoch_breaks, epoch_rate);
    }

    double to_generations(const double& coalescents) const
    {
      if (coalescents == 0.0) return 0.0;
      auto index = 
        std::lower_bound(coal_breaks.begin(), coal_breaks.end(), coalescents) - 
        coal_breaks.begin() - 1; // O(log dim)
      return coalescents / rate[index] + coal_step[index];
    }

    double to_coalescents(const double& generations) const
    {
      if (generations == 0.0) return 0.0;
      auto index = 
        std::lower_bound(gens_breaks.begin(), gens_breaks.end(), generations) - 
        gens_breaks.begin() - 1; // O(log dim)
      return generations * rate[index] + gens_step[index];
    }
};


struct sum_decomposition_t
{
  /* Decomposition of node marginals into differences between parent (the node) and the maximum of its children */

  private:
    gamma_distr_t maximum_moments (const gamma_distr_t& a, const gamma_distr_t& b)
    {
      /* calculate moments of maximum of two gamma distributions */
      if (a && b) { // both children non-fixed
        double mn, sq, u, v, x;
        x = a.rate() / (a.rate() + b.rate());
        u = betainc(b.shape(), a.shape() + 1, 1.0 - x);
        v = betainc(a.shape(), b.shape() + 1, x);
        mn = a.mean() * u + b.mean() * v;
        u = betainc(b.shape(), a.shape() + 2, 1.0 - x);
        v = betainc(a.shape(), b.shape() + 2, x);
        sq = a.square() * u + b.square() * v;
        return gamma_distr_t::from_moments(mn, sq - std::pow(mn, 2));
      } else if (a) {
        return a;
      } else if (b) {
        return b;
      } else {
        return gamma_distr_t();
      }
    }

    gamma_distr_t gamma_difference (const gamma_distr_t& a, const gamma_distr_t& b)
    {
      /* gamma distribution where E[c] = E[a] - E[b] and V[c] = V[a] - V[b] */
      return gamma_distr_t::from_moments(a.mean() - b.mean(), a.variance() - b.variance());
    }

  public:
    std::vector<gamma_distr_t> factors;

    sum_decomposition_t (const std::vector<trio_t>& trios, const std::vector<gamma_distr_t>& marginals) 
      : factors(2 * trios.size() + 1)
    {
      assert (marginals.size() == factors.size());
      gamma_distr_t eldest;
      for (auto& i : trios) {
        eldest = maximum_moments(marginals[i.lchild], marginals[i.rchild]);
        factors[i.parent] = gamma_difference(marginals[i.parent], eldest);
      }
    }
};


struct conditional_coalescent_prior_t
{
  /* Topology-conditioned coalescent prior, on the coalescent timescale */

  private:

    std::vector<double> pr_node_is_event (std::vector<trio_t> trios)
    {
      /* u[i][k] is the proportion of permissible permutations where node i is event k.
       * a permutation is permissible if it satisfies the partial ordering imposed by the tree.
       * assumes that "trios" is in order (children before parents).
       * returned in flattened row-major form.
       */

      double a, b;
      std::size_t E = trios.size(), N = 2 * E + 1, root = trios.back().parent;

      /* find number of internal descendants */
      std::vector<std::size_t> n (N, 0);
      for (auto& i : trios) {
        n[i.parent] = 1 + n[i.lchild] + n[i.rchild];
      }
      assert (n[root] == E);

      /* recursively calculate state probabilities */
      std::reverse(trios.begin(), trios.end());
      std::vector<double> u (N*E, 0.0); // row-major order
      u[root*E + 0] = 1.0;
      for (auto& i : trios) {
        for (auto& c : {i.lchild, i.rchild}) {
          for (std::size_t t = 0; t < E - 1; ++t) {
            a = (E - t - n[c]) / (E - t - 1.0);
            b = n[c] / (E - t - 1.0);
            u[E*c+t+1] = a * u[E*c+t] + b * u[E*i.parent+t];
          }
        }
      }

      return u;
    }

  public:
    std::vector<gamma_distr_t> edges, nodes;

    conditional_coalescent_prior_t (const std::vector<trio_t>& trios)
      : edges (2 * trios.size() + 1)
      , nodes (2 * trios.size() + 1)
    {
      std::size_t E = trios.size(), N = 2 * E + 1;

      /* gamma projections for hypoexponential event times */
      std::vector<double> mean, logmean;
      mean.reserve(E); logmean.reserve(E);
      std::size_t k = E + 1;
      double la, ln, mn = 0.0, va = 0.0;
      for (std::size_t i = 0; i < E; ++i, --k) {
        la = 2.0 / (k * (k - 1.0));
        mn += la; va += std::pow(la, 2);
        ln = log(mn) - 0.5 * va / std::pow(mn, 2); // taylor approximation
        mean.push_back(mn); logmean.push_back(ln);
      }

      /* marginalize over node-event probabilities */
      auto P = pr_node_is_event(trios);
      for (auto& i : trios) {
        auto row = P.begin() + E*i.parent;
        nodes[i.parent] = gamma_distr_t::from_sufficient(
            std::inner_product(row, row + E, mean.rbegin(), 0.0),
            std::inner_product(row, row + E, logmean.rbegin(), 0.0)
        );
      }

      /* convert to prior on edge from eldest child to parent */
      edges = sum_decomposition_t(trios, nodes).factors;
    }
};


// -------------- MAIN ALGORITHM ------------------ //

template <class RNG>
struct expectation_propagation_t
{
  /* Date a tree (array of parent/child trios) using expectation propagation */

  private:
    std::vector<gamma_distr_t> parent_message, lchild_message, rchild_message;

  public:

    /* --- QUADRATURE SCHEME (VARIABLE NE) --- */

    static
    double quadrature_conditional(
        gamma_distr_t& p_i, // posterior (output)
        const gamma_distr_t& q_i, // cavity (input)
        const gamma_distr_t& prior, const gamma_distr_t& l_j, const gamma_distr_t& l_k, // likelihood and prior
        const timescale_t& timescale, const std::size_t order
    ) {
      /* Estimate moments of surrogate via Gaussian quadrature, conditioned on
         t_i > t_j = t_k = 0 */

      auto loglik = [&l_j,&l_k,&q_i] (const double t_i, const double t_j, const double t_k) {
        return l_j.loglik(t_i - t_j) + l_k.loglik(t_i - t_k) + q_i.loglik(t_i);
      };

      double c_i, t_i, t_j, t_k, lp, alpha, beta, scale;
      const std::size_t dim = order;
      std::vector<double> weight, age_i, log_i;
      weight.reserve(dim); age_i.reserve(dim); log_i.reserve(dim);

      assert (prior.alpha > -1.0);
      assert (prior.beta < 0.0);
      scale = log(-prior.beta) * (-prior.alpha - 1.0);
      gauss_laguerre_t quad_i (order, prior.alpha);
      for (auto& i : quad_i.grid) {
        c_i = -i.abscissa / prior.beta;
        t_i = timescale.to_generations(c_i);
        t_j = 0.0; // TODO use sample ages
        t_k = 0.0; // TODO use sample ages
        lp = loglik(t_i, t_j, t_k) + i.logweight + scale;
        weight.push_back(lp); age_i.push_back(t_i); log_i.push_back(log(t_i));
      }

      /* normalizing constant */
      double maxlp, norm, logconst;
      maxlp = *std::max_element(weight.begin(), weight.end());
      std::for_each(weight.begin(), weight.end(), [&maxlp](double& x){ x = exp(x - maxlp); });
      norm = std::accumulate(weight.begin(), weight.end(), 0.0);
      logconst = log(norm) + maxlp;

      /* project onto gamma */
      double mn_i, ln_i;
      mn_i = std::inner_product(weight.begin(), weight.end(), age_i.begin(), 0.0) / norm;
      ln_i = std::inner_product(weight.begin(), weight.end(), log_i.begin(), 0.0) / norm;
      p_i = gamma_distr_t::from_sufficient(mn_i, ln_i);

      return logconst;
    }


    static
    double quadrature_conditional(
        gamma_distr_t& p_i, gamma_distr_t& p_j, // posterior (output)
        const gamma_distr_t& q_i, const gamma_distr_t& q_j, // cavities (input)
        const gamma_distr_t& prior, const gamma_distr_t& l_j, const gamma_distr_t& l_k, // likelihood and prior
        const timescale_t& timescale, const std::size_t order
    ) {
      /* Estimate moments of surrogate via Gaussian quadrature, conditioned on
         t_i > t_j and t_k = 0. */

      /* loglikelihood on natural timescale */
      auto loglik = [&l_j,&l_k,&q_i] (const double t_i, const double t_j, const double t_k) {
        return l_j.loglik(t_i - t_j) + l_k.loglik(t_i - t_k) + q_i.loglik(t_i);
      };

      /* absorb as much as possible into quadrature weight */
      auto proposal = [&l_j,&l_k,&q_i] (const double t_i, const double t_j, const double t_k) {
        return log(t_j)*(q_i.alpha + l_k.alpha) + t_j*(q_i.beta + l_k.beta);
      };

      /* quadrature over ages of parent and older child */
      double c_ij, c_j, t_i, t_j, t_k, lp, alpha, beta, scale;
      const std::size_t dim = std::pow(order, 2);
      std::vector<double> weight, age_i, log_i, age_j, log_j;
      weight.reserve(dim); age_i.reserve(dim); log_i.reserve(dim);
      age_j.reserve(dim); log_j.reserve(dim);

      gauss_laguerre_t quad_ij (order, prior.alpha);
      alpha = q_j.alpha + q_i.alpha + l_k.alpha;
      beta = q_j.beta + q_i.beta + l_k.beta;
      gauss_laguerre_t quad_j (order, alpha);
      assert (beta < 0.0 && prior.beta < 0.0);
      scale = log(-beta) * (-alpha - 1.0) + log(-prior.beta) * (-prior.alpha - 1.0);
      for (auto& j : quad_j.grid) {
        t_j = -j.abscissa / beta;
        c_j = timescale.to_coalescents(t_j);
        for (auto& i : quad_ij.grid) {
          c_ij = -i.abscissa / prior.beta;
          t_i = timescale.to_generations(c_j + c_ij);
          t_k = 0.0; //TODO use sample age
          lp = loglik(t_i, t_j, t_k) - proposal(t_i, t_j, t_k) + i.logweight + j.logweight + scale;
          weight.push_back(lp);
          age_i.push_back(t_i); log_i.push_back(log(t_i));
          age_j.push_back(t_j); log_j.push_back(log(t_j));
        }
      }

      /* normalizing constant */
      double maxlp, norm, logconst;
      maxlp = *std::max_element(weight.begin(), weight.end());
      std::for_each(weight.begin(), weight.end(), [&maxlp](double& x){ x = exp(x - maxlp); });
      norm = std::accumulate(weight.begin(), weight.end(), 0.0);
      logconst = log(norm) + maxlp;

      /* project onto gamma */
      double mn_i, ln_i, mn_j, ln_j;
      mn_i = std::inner_product(weight.begin(), weight.end(), age_i.begin(), 0.0) / norm;
      ln_i = std::inner_product(weight.begin(), weight.end(), log_i.begin(), 0.0) / norm;
      mn_j = std::inner_product(weight.begin(), weight.end(), age_j.begin(), 0.0) / norm;
      ln_j = std::inner_product(weight.begin(), weight.end(), log_j.begin(), 0.0) / norm;
      p_i = gamma_distr_t::from_sufficient(mn_i, ln_i);
      p_j = gamma_distr_t::from_sufficient(mn_j, ln_j);

      return logconst;
    }

    static
    double quadrature_conditional(
        gamma_distr_t& p_i, gamma_distr_t& p_j, gamma_distr_t& p_k, // posterior (output)
        const gamma_distr_t& q_i, const gamma_distr_t& q_j, const gamma_distr_t& q_k,  // cavities (input)
        const gamma_distr_t& prior, const gamma_distr_t& l_j, const gamma_distr_t& l_k, // likelihood and prior
        const timescale_t& timescale, const std::size_t order
    ) {
      /* Estimate moments of surrogate via Gaussian quadrature, conditioned on
         t_i > t_j > t_k. */

      /* loglikelihood on natural timescale */
      auto loglik = [&l_j,&l_k,&q_i] (const double t_i, const double t_j, const double t_k) {
        return l_j.loglik(t_i - t_j) + l_k.loglik(t_i - t_k) + q_i.loglik(t_i);
      };

      /* absorb as much as possible into quadrature weight */
      auto proposal = [&l_j,&l_k,&q_i] (const double t_i, const double t_j, const double t_k) {
        return log(t_j)*(q_i.alpha + l_k.alpha) + t_j*q_i.beta + (t_j - t_k)*l_k.beta;
      };

      /* quadrature over ages of parent and older child */
      double c_ij, c_j, t_i, t_j, t_k, lp, alpha, beta, scale;
      const std::size_t dim = std::pow(order, 3);
      std::vector<double> weight, age_i, log_i, age_j, log_j, age_k, log_k; 
      weight.reserve(dim); age_i.reserve(dim); log_i.reserve(dim);
      age_j.reserve(dim); log_j.reserve(dim); age_k.reserve(dim); log_k.reserve(dim);

      gauss_laguerre_t quad_ij (order, prior.alpha);
      alpha = q_i.alpha + q_j.alpha + q_k.alpha + l_k.alpha + 1.0;
      gauss_laguerre_t quad_j (order, alpha);
      gauss_jacobi_t quad_k (order, 0, q_k.alpha);
      for (auto& k : quad_k.grid) {
        for (auto& j : quad_j.grid) {
          beta = q_i.beta + q_j.beta + q_k.beta*k.abscissa + l_k.beta*(1.0 - k.abscissa);
          assert (prior.beta < 0.0);
          assert (beta < 0.0);
          scale = log(-beta) * (-alpha - 1.0) + log(-prior.beta) * (-prior.alpha - 1);
          t_j = -j.abscissa / beta;
          c_j = timescale.to_coalescents(t_j);
          t_k = k.abscissa * t_j;
          for (auto& ij : quad_ij.grid) {
            c_ij = -ij.abscissa / prior.beta;
            t_i = timescale.to_generations(c_j + c_ij);
            lp = loglik(t_i, t_j, t_k) - proposal(t_i, t_j, t_k) + ij.logweight + j.logweight + k.logweight + scale;
            weight.push_back(lp);
            age_i.push_back(t_i); log_i.push_back(log(t_i));
            age_j.push_back(t_j); log_j.push_back(log(t_j));
            age_k.push_back(t_k); log_k.push_back(log(t_k));
          }
        }
      }

      /* normalizing constant */
      double maxlp, norm, logconst;
      maxlp = *std::max_element(weight.begin(), weight.end());
      std::for_each(weight.begin(), weight.end(), [&maxlp](double& x){ x = exp(x - maxlp); });
      norm = std::accumulate(weight.begin(), weight.end(), 0.0);
      logconst = log(norm) + maxlp;

      /* project onto gamma */
      double mn_i, ln_i, mn_j, ln_j, mn_k, ln_k;
      mn_i = std::inner_product(weight.begin(), weight.end(), age_i.begin(), 0.0) / norm;
      ln_i = std::inner_product(weight.begin(), weight.end(), log_i.begin(), 0.0) / norm;
      mn_j = std::inner_product(weight.begin(), weight.end(), age_j.begin(), 0.0) / norm;
      ln_j = std::inner_product(weight.begin(), weight.end(), log_j.begin(), 0.0) / norm;
      mn_k = std::inner_product(weight.begin(), weight.end(), age_k.begin(), 0.0) / norm;
      ln_k = std::inner_product(weight.begin(), weight.end(), log_k.begin(), 0.0) / norm;
      p_i = gamma_distr_t::from_sufficient(mn_i, ln_i);
      p_j = gamma_distr_t::from_sufficient(mn_j, ln_j);
      p_k = gamma_distr_t::from_sufficient(mn_k, ln_k);

      return logconst;
    }


    /* --- QUADRATURE SCHEME (CONSTANT NE) --- */

    static
    double quadrature_conditional(
        gamma_distr_t& p_i, // posterior (output)
        const gamma_distr_t& q_i, // cavity (input)
        const gamma_distr_t& prior, const gamma_distr_t& l_j, const gamma_distr_t& l_k, // likelihood and prior
        const double coal_rate, const std::size_t order
    ) {
      /* Estimate moments of surrogate via Gaussian quadrature, conditioned on
         t_i > t_j = t_k = 0. The integral (omitting the prior) */

      /* NB: no need for quadrature when both children are contemporary leaves, but
         will be needed for ancient samples */
      double t_i, lp, alpha, beta, scale;
      const std::size_t dim = order;
      std::vector<double> weight, age_i, log_i;
      weight.reserve(dim); age_i.reserve(dim); log_i.reserve(dim);

      alpha = q_i.alpha + l_k.alpha + l_j.alpha + prior.alpha;
      beta = q_i.beta + l_k.beta + l_j.beta + prior.beta*coal_rate;
      assert (alpha > -1.0);
      assert (beta < 0.0);
      scale = log(-beta) * (-alpha - 1.0);
      gauss_laguerre_t quad_i (order, alpha);
      for (auto& i : quad_i.grid) {
        t_i = -i.abscissa / beta;
        lp = i.logweight + scale;
        weight.push_back(lp); age_i.push_back(t_i); log_i.push_back(log(t_i));
      }

      /* normalizing constant */
      double maxlp, norm, logconst;
      maxlp = *std::max_element(weight.begin(), weight.end());
      std::for_each(weight.begin(), weight.end(), [&maxlp](double& x){ x = exp(x - maxlp); });
      norm = std::accumulate(weight.begin(), weight.end(), 0.0);
      logconst = log(norm) + maxlp;

      /* project onto gamma */
      double mn_i, ln_i;
      mn_i = std::inner_product(weight.begin(), weight.end(), age_i.begin(), 0.0) / norm;
      ln_i = std::inner_product(weight.begin(), weight.end(), log_i.begin(), 0.0) / norm;
      p_i = gamma_distr_t::from_sufficient(mn_i, ln_i);

      return logconst;
    }

    static
    double quadrature_conditional(
        gamma_distr_t& p_i, gamma_distr_t& p_j, // posterior (output)
        const gamma_distr_t& q_i, const gamma_distr_t& q_j, // cavities (input)
        const gamma_distr_t& prior, const gamma_distr_t& l_j, const gamma_distr_t& l_k, // likelihood and prior
        const double coal_rate, const std::size_t order
    ) {
      /* Estimate moments of surrogate via Gaussian quadrature, conditioned on
         t_i > t_j and t_k = 0. */

      double t_i, t_j, lp, alpha, beta, scale;
      const std::size_t dim = std::pow(order, 2);
      std::vector<double> weight, age_i, log_i, age_j, log_j;
      weight.reserve(dim); age_i.reserve(dim); log_i.reserve(dim);
      age_j.reserve(dim); log_j.reserve(dim);

      alpha = q_i.alpha + q_j.alpha + l_k.alpha + l_j.alpha + prior.alpha + 1.0;
      gauss_laguerre_t quad_i (order, alpha);
      gauss_jacobi_t quad_j (order, l_j.alpha + prior.alpha, q_j.alpha);
      for (auto& j : quad_j.grid) {
        beta = q_i.beta + l_k.beta + q_j.beta*j.abscissa +
          (l_j.beta + prior.beta*coal_rate)*(1.0 - j.abscissa);
        assert (beta < 0.0);
        scale = log(-beta) * (-alpha - 1.0);
        for (auto& i : quad_i.grid) {
          t_i = -i.abscissa / beta;
          t_j = j.abscissa * t_i;
          lp = i.logweight + j.logweight + scale;
          weight.push_back(lp);
          age_i.push_back(t_i); log_i.push_back(log(t_i));
          age_j.push_back(t_j); log_j.push_back(log(t_j));
        }
      }

      /* normalizing constant */
      double maxlp, norm, logconst;
      maxlp = *std::max_element(weight.begin(), weight.end());
      std::for_each(weight.begin(), weight.end(), [&maxlp](double& x){ x = exp(x - maxlp); });
      norm = std::accumulate(weight.begin(), weight.end(), 0.0);
      logconst = log(norm) + maxlp;

      /* project onto gamma */
      double mn_i, ln_i, mn_j, ln_j;
      mn_i = std::inner_product(weight.begin(), weight.end(), age_i.begin(), 0.0) / norm;
      ln_i = std::inner_product(weight.begin(), weight.end(), log_i.begin(), 0.0) / norm;
      mn_j = std::inner_product(weight.begin(), weight.end(), age_j.begin(), 0.0) / norm;
      ln_j = std::inner_product(weight.begin(), weight.end(), log_j.begin(), 0.0) / norm;
      p_i = gamma_distr_t::from_sufficient(mn_i, ln_i);
      p_j = gamma_distr_t::from_sufficient(mn_j, ln_j);

      return logconst;
    }

    static
    double quadrature_conditional(
        gamma_distr_t& p_i, gamma_distr_t& p_j, gamma_distr_t& p_k, // posterior (output)
        const gamma_distr_t& q_i, const gamma_distr_t& q_j, const gamma_distr_t& q_k,  // cavities (input)
        const gamma_distr_t& prior, const gamma_distr_t& l_j, const gamma_distr_t& l_k, // likelihood and prior
        const double coal_rate, const std::size_t order
    ) {
      /* Estimate moments of surrogate via Gaussian quadrature, conditioned on
         t_i > t_j > t_k. */

      double t_i, t_j, t_k, lp, alpha, beta, scale;
      const std::size_t dim = std::pow(order, 3);
      std::vector<double> weight, age_i, log_i, age_j, log_j, age_k, log_k; 
      weight.reserve(dim); age_i.reserve(dim); log_i.reserve(dim);
      age_j.reserve(dim); log_j.reserve(dim); age_k.reserve(dim); log_k.reserve(dim);

      alpha = q_i.alpha + q_k.alpha + q_j.alpha + l_k.alpha + l_j.alpha + prior.alpha + 2;
      gauss_laguerre_t quad_i (order, alpha);
      gauss_jacobi_t quad_j (order, l_j.alpha + prior.alpha, q_j.alpha + q_k.alpha + 1);
      gauss_jacobi_t quad_k (order, 0.0, q_k.alpha);
      for (auto& j : quad_j.grid) {
        for (auto& k : quad_k.grid) {
          beta = q_i.beta + q_j.beta*j.abscissa + 
            (l_j.beta + prior.beta*coal_rate)*(1 - j.abscissa) + 
            q_k.beta*j.abscissa*k.abscissa + l_k.beta*(1 - j.abscissa*k.abscissa);
          assert (beta < 0.0);
          scale = log(-beta) * (-alpha - 1.0); 
          for (auto& i : quad_i.grid) {
            t_i = -i.abscissa / beta;
            t_j = j.abscissa * t_i;
            t_k = k.abscissa * t_j;
            lp = log(1.0 - j.abscissa * k.abscissa)*l_k.alpha + 
              i.logweight + j.logweight + k.logweight + scale;
            weight.push_back(lp);
            age_i.push_back(t_i); log_i.push_back(log(t_i));
            age_j.push_back(t_j); log_j.push_back(log(t_j));
            age_k.push_back(t_k); log_k.push_back(log(t_k));
          }
        }
      }

      /* normalizing constant */
      double maxlp, norm, logconst;
      maxlp = *std::max_element(weight.begin(), weight.end());
      std::for_each(weight.begin(), weight.end(), [&maxlp](double& x){ x = exp(x - maxlp); });
      norm = std::accumulate(weight.begin(), weight.end(), 0.0);
      logconst = log(norm) + maxlp;

      /* project onto gamma */
      double mn_i, ln_i, mn_j, ln_j, mn_k, ln_k;
      mn_i = std::inner_product(weight.begin(), weight.end(), age_i.begin(), 0.0) / norm;
      ln_i = std::inner_product(weight.begin(), weight.end(), log_i.begin(), 0.0) / norm;
      mn_j = std::inner_product(weight.begin(), weight.end(), age_j.begin(), 0.0) / norm;
      ln_j = std::inner_product(weight.begin(), weight.end(), log_j.begin(), 0.0) / norm;
      mn_k = std::inner_product(weight.begin(), weight.end(), age_k.begin(), 0.0) / norm;
      ln_k = std::inner_product(weight.begin(), weight.end(), log_k.begin(), 0.0) / norm;
      p_i = gamma_distr_t::from_sufficient(mn_i, ln_i);
      p_j = gamma_distr_t::from_sufficient(mn_j, ln_j);
      p_k = gamma_distr_t::from_sufficient(mn_k, ln_k);

      return logconst;
    }

    static
    double quadrature(
        gamma_distr_t& p_i, gamma_distr_t& p_j, gamma_distr_t& p_k, // posterior (output)
        const gamma_distr_t& q_i, const gamma_distr_t& q_j, const gamma_distr_t& q_k,  // cavities (input)
        const gamma_distr_t& prior, const gamma_distr_t& l_j, const gamma_distr_t& l_k, // likelihood and prior
        const timescale_t& timescale, const std::size_t order
    ) {
      /* Estimate moments of surrogate via quadrature. The mutation counts (l_j.alpha, l_k.alpha)
         are silently converted to integers. */

      //std::cout << q_i << " " << q_j << " " << q_k << " " << prior << " " << l_j << " " << l_k << std::endl;//DEBUG

      const bool demography = timescale.dim > 1; // should make this toggle-able
      double logconst;

      if (q_j && q_k) { /* ages of both children not fixed */
        gamma_distr_t p_i_jk, p_j_jk, p_k_jk, p_i_kj, p_j_kj, p_k_kj;
        double lp_jk, lp_kj, w_jk, w_kj, maxlp, mn_i, ln_i, mn_j, ln_j, mn_k, ln_k;
        /* conditioning on t_i > t_j > t_k */
        lp_jk = demography ?
          quadrature_conditional(p_i_jk, p_j_jk, p_k_jk, q_i, q_j, q_k, prior, l_j, l_k, timescale, order) :
          quadrature_conditional(p_i_jk, p_j_jk, p_k_jk, q_i, q_j, q_k, prior, l_j, l_k, 1.0 / timescale.pairwise_tmrca, order);
        /* conditioning on t_i > t_k > t_j */
        lp_kj = demography ?
          quadrature_conditional(p_i_kj, p_k_kj, p_j_kj, q_i, q_k, q_j, prior, l_k, l_j, timescale, order) :
          quadrature_conditional(p_i_kj, p_k_kj, p_j_kj, q_i, q_k, q_j, prior, l_k, l_j, 1.0 / timescale.pairwise_tmrca, order);
        /* sum over conditionals */
        maxlp = std::max(lp_jk, lp_kj);
        logconst = log(exp(lp_jk - maxlp) + exp(lp_kj - maxlp)) + maxlp;
        w_jk = exp(lp_jk - logconst);
        w_kj = exp(lp_kj - logconst);
        /* project to gamma */
        mn_i = p_i_jk.mean() * w_jk + p_i_kj.mean() * w_kj;
        mn_j = p_j_jk.mean() * w_jk + p_j_kj.mean() * w_kj;
        mn_k = p_k_jk.mean() * w_jk + p_k_kj.mean() * w_kj;
        ln_i = p_i_jk.logmean() * w_jk + p_i_kj.logmean() * w_kj;
        ln_j = p_j_jk.logmean() * w_jk + p_j_kj.logmean() * w_kj;
        ln_k = p_k_jk.logmean() * w_jk + p_k_kj.logmean() * w_kj;
        p_i = gamma_distr_t::from_sufficient(mn_i, ln_i);
        p_j = gamma_distr_t::from_sufficient(mn_j, ln_j);
        p_k = gamma_distr_t::from_sufficient(mn_k, ln_k);
      } else if (q_j) {
        /* conditioning on t_i > t_j > t_k = 0 */
        logconst = demography ?
          quadrature_conditional(p_i, p_j, q_i, q_j, prior, l_j, l_k, timescale, order) :
          quadrature_conditional(p_i, p_j, q_i, q_j, prior, l_j, l_k, 1.0 / timescale.pairwise_tmrca, order);
      } else if (q_k) {
        /* conditioning on t_i > t_k > t_j = 0 */
        logconst = demography ? 
          quadrature_conditional(p_i, p_k, q_i, q_k, prior, l_k, l_j, timescale, order) :
          quadrature_conditional(p_i, p_k, q_i, q_k, prior, l_k, l_j, 1.0 / timescale.pairwise_tmrca, order);
      } else {
        /* conditioning on t_i > t_k = t_j = 0 */
        logconst = demography ?
          quadrature_conditional(p_i, q_i, prior, l_k, l_j, timescale, order) :
          quadrature_conditional(p_i, q_i, prior, l_k, l_j, 1.0 / timescale.pairwise_tmrca, order);
      }

      return logconst;
    }

    double update (const trio_t& t, const std::vector<gamma_distr_t>& likelihood, const std::vector<gamma_distr_t>& prior, const timescale_t& timescale, const std::size_t order)
    {
      /* update posterior via cavity method, using quadrature.
         complexity is O(order^2) and solution is near exact */
      double logconst;
      gamma_distr_t parent_cavity = posterior[t.parent] - parent_message[t.parent];
      gamma_distr_t lchild_cavity = posterior[t.lchild] - lchild_message[t.parent];
      gamma_distr_t rchild_cavity = posterior[t.rchild] - rchild_message[t.parent];
      logconst = quadrature(
        posterior[t.parent], posterior[t.lchild], posterior[t.rchild],
        parent_cavity, lchild_cavity, rchild_cavity,
        prior[t.parent], likelihood[t.lchild], likelihood[t.rchild],
        timescale, order
      );
      parent_message[t.parent] = posterior[t.parent] - parent_cavity;
      lchild_message[t.parent] = posterior[t.lchild] - lchild_cavity;
      rchild_message[t.parent] = posterior[t.rchild] - rchild_cavity;
      return logconst;
    }


    /* --- IMPORTANCE SAMPLING SCHEME --- */

    static
    double importance_sample(
        gamma_distr_t& p_i, gamma_distr_t& p_j, gamma_distr_t& p_k, // posterior (output)
        const gamma_distr_t& q_i, const gamma_distr_t& q_j, const gamma_distr_t& q_k,  // cavities (input)
        const gamma_distr_t& prior, const gamma_distr_t& l_j, const gamma_distr_t& l_k, // likelihood and prior
        const timescale_t& timescale, const std::size_t n, RNG& rng
    ) {
      /* Estimate moments of surrogate via importance sampling */

      // TODO make this streaming until desired ESS rather than fixed storage
      double t_x, c_i;
      auto rng_i = prior.random(), rng_j = q_j.random(), rng_k = q_k.random();
      std::vector<double> t_i(n), t_j(n), t_k(n), w(n), lp(n);
      std::fill(w.begin(), w.end(), -std::numeric_limits<double>::infinity());
      std::fill(lp.begin(), lp.end(), -std::numeric_limits<double>::infinity());
      for (std::size_t itt = 0; itt < n; ++itt) {
        /* child ages sampled from variational cavity */
        t_j[itt] = rng_j(rng); 
        t_k[itt] = rng_k(rng);

        /* branch length from eldest child to parent sampled from coalescent */
        c_i = rng_i(rng);
        t_x = std::max(t_j[itt], t_k[itt]);
        t_i[itt] = timescale.to_generations(timescale.to_coalescents(t_x) + c_i);

        if (t_i[itt] - t_x > IS_MINGAP) {
          /* likelihood weights */
          w[itt] = q_i.loglik(t_i[itt]); // variational prior for parent
          w[itt] += l_j.loglik(t_i[itt] - t_j[itt]); // mutation likelihood
          w[itt] += l_k.loglik(t_i[itt] - t_k[itt]); // mutation likelihood

          /* normalizing constant */
          lp[itt] = w[itt] + q_j.loglik(t_j[itt]) + 
            q_k.loglik(t_k[itt]) + prior.loglik(c_i);
        }
      }

      /* importance sampling */
      double maxlp, maxw, norm, logconst;
      maxlp = *std::max_element(lp.begin(), lp.end());
      logconst = std::accumulate(lp.begin(), lp.end(), 0.0, 
          [&maxlp](double x, double y) { return x + exp(y - maxlp);}
      );
      logconst = log(logconst) + maxlp; // FIXME this isn't correct
      maxw = *std::max_element(w.begin(), w.end());
      std::for_each(w.begin(), w.end(), [&maxw](double& x){x = exp(x - maxw);});
      norm = std::accumulate(w.begin(), w.end(), 0.0);
      std::for_each(w.begin(), w.end(), [&norm](double& x){x /= norm;});

      /* approximately minimize KL using simulated moments */
      double ln_i = 0.0, mn_i = 0.0, ln_j = 0.0, mn_j = 0.0, ln_k = 0.0, mn_k = 0.0;
      double ess = 0.0;
      for (std::size_t itt = 0; itt < n; ++itt) {
        ln_i += w[itt] * log(t_i[itt]); mn_i += w[itt] * t_i[itt];
        ln_j += w[itt] * log(t_j[itt]); mn_j += w[itt] * t_j[itt];
        ln_k += w[itt] * log(t_k[itt]); mn_k += w[itt] * t_k[itt];
        ess += std::pow(w[itt], 2);
      }
      ess = 1.0 / ess; // TODO: do something with ess, like asymptotic std err
                       
      /* gamma projection */
      p_i = gamma_distr_t::from_sufficient(mn_i, ln_i);
      p_j = gamma_distr_t::from_sufficient(mn_j, ln_j);
      p_k = gamma_distr_t::from_sufficient(mn_k, ln_k);

      return logconst;
    }


    double update (const trio_t& t, const std::vector<gamma_distr_t>& likelihood, const std::vector<gamma_distr_t>& prior, const timescale_t& timescale, const std::size_t num_samples, RNG& random_generator) 
    {
      /* update posterior via cavity method using importance sampling.
         complexity is O(num_samples) and solution is approximate */
      double logconst;
      gamma_distr_t parent_cavity = posterior[t.parent] - parent_message[t.parent];
      gamma_distr_t lchild_cavity = posterior[t.lchild] - lchild_message[t.parent];
      gamma_distr_t rchild_cavity = posterior[t.rchild] - rchild_message[t.parent];
      logconst = importance_sample(
        posterior[t.parent], posterior[t.lchild], posterior[t.rchild],
        parent_cavity, lchild_cavity, rchild_cavity,
        prior[t.parent], likelihood[t.lchild], likelihood[t.rchild],
        timescale, num_samples, random_generator
      );
      parent_message[t.parent] = posterior[t.parent] - parent_cavity;
      lchild_message[t.parent] = posterior[t.lchild] - lchild_cavity;
      rchild_message[t.parent] = posterior[t.rchild] - rchild_cavity;
      return logconst;
    }

  /* ------------ API ------------ */

  public:
    std::vector<gamma_distr_t> posterior;

    /* use importance sampling */
    expectation_propagation_t (const std::vector<trio_t>& trios, const std::vector<gamma_distr_t>& likelihood, const std::vector<gamma_distr_t>& prior, const timescale_t& timescale, const std::size_t num_samples, RNG& random_generator)
      : posterior (2 * trios.size() + 1)
      , parent_message (2 * trios.size() + 1)
      , lchild_message (2 * trios.size() + 1)
      , rchild_message (2 * trios.size() + 1)
    {
      assert (num_samples > 1);

      /* forward pass, from leaves to root */
      for (auto t = trios.begin(); t != trios.end(); ++t) {
        update(*t, likelihood, prior, timescale, num_samples, random_generator);
      }

      /* backward pass, from root to leaves */
      for (auto t = trios.rbegin() + 1; t != trios.rend(); ++t) {
        update(*t, likelihood, prior, timescale, num_samples, random_generator);
      }
    }

    /* use quadrature */
    expectation_propagation_t (const std::vector<trio_t>& trios, const std::vector<gamma_distr_t>& likelihood, const std::vector<gamma_distr_t>& prior, const timescale_t& timescale, const std::size_t order)
      : posterior (2 * trios.size() + 1)
      , parent_message (2 * trios.size() + 1)
      , lchild_message (2 * trios.size() + 1)
      , rchild_message (2 * trios.size() + 1)
    {
      assert (order > 1);

      /* forward pass, from leaves to root */
      for (auto t = trios.begin(); t != trios.end(); ++t) {
        update(*t, likelihood, prior, timescale, order);
      }

      /* backward pass, from root to leaves */
      for (auto t = trios.rbegin() + 1; t != trios.rend(); ++t) {
        update(*t, likelihood, prior, timescale, order);
      }
    }
};


// -------------- RELATE-FACING API ------------------ //

struct EstimateBranchLengthsVariational
{
  /* 
     Date trees via expectation propagation.  Moment matching is done either
     via importance sampling or quadrature, toggled by the choice of
     constructor. In general quadrature should be preferred, but importance
     sampling may be useful as a fallback.

     We're also assuming that all samples are contemporary. This could be
     possibly be relaxed in the future.
  */

  private:
  Data* const data;
  std::mt19937 rng;
  std::size_t order = 0; // order of quadrature rule
  std::size_t importance_samples = 0; // number of importance samples
  timescale_t timescale;

  double mutational_target(const Node& node) const {
    /* calculate node span shared across trees */
    double span = 0.0;
    for (auto j = node.SNP_begin; j < node.SNP_end; ++j) {
      span += data->dist[j];
    }
    if (node.SNP_begin > 0) span += 0.5 * data->dist[node.SNP_begin - 1];
    if (node.SNP_end < data->L - 1) span += 0.5 * data->dist[node.SNP_end];
    span *= data->mu;
    return span;
  }

  public:

  static
  std::vector<double> ages_from_lengths(const Tree& tree, const std::vector<double>& branch_lengths)
  {
    /* extract parent-child trios */
    std::vector<trio_t> trios; trios.reserve(tree.nodes.size());
    for (auto& i : tree.nodes) {
      if (i.child_left != NULL and i.child_right != NULL) {
        trios.emplace_back(i.label, i.child_left->label, i.child_right->label);
      }
    }
  
    /* order trios by leaves subtended by parent */
    std::vector<Leaves> leaves; tree.FindAllLeaves(leaves);
    std::sort(trios.begin(), trios.end(), [&leaves](const trio_t& x, const trio_t& y){
        return leaves[x.parent].num_leaves < leaves[y.parent].num_leaves;
    });

    std::size_t E = trios.size(), N = 2*E + 1;

    std::vector<double> ages (N, 0.0);
    assert (branch_lengths.size() == N); //length of branch subtended by node
    for (auto& i : trios) {
      ages[i.parent] = (ages[i.lchild] + branch_lengths[i.lchild] + ages[i.rchild] + branch_lengths[i.rchild]) / 2.0;
    }

    return ages;
  }

  /* ------- CONSTRUCTORS ------- */

  /* use quadrature and fixed coalescence rate */
  EstimateBranchLengthsVariational (Data* const data, const double& coalescence_rate, const std::size_t order)
    : data (data)
    , order (order)
    , timescale(std::vector<double>({0.0}), std::vector<double>({coalescence_rate}))
  {}

  /* use quadrature and variable coalescence rate */
  EstimateBranchLengthsVariational (Data* const data, const std::vector<double>& epoch_start, const std::vector<double>& epoch_rate, const std::size_t order)
    : data (data)
    , order (order)
    , timescale(epoch_start, epoch_rate)
  {
    //for (auto& i : epoch_start) { std::cout << "epoch_start.emplace_back(" << i << ");" << std::endl; }//DEBUG
    //for (auto& i : epoch_rate) { std::cout << "epoch_rate.emplace_back(" << i << ");" << std::endl; }//DEBUG
  }

  /* use importance sampling and fixed coalescence rate */
  EstimateBranchLengthsVariational (Data* const data, const double& coalescence_rate, const std::size_t importance_samples, const std::size_t random_seed)
    : data (data)
    , rng (random_seed)
    , importance_samples (importance_samples)
    , timescale(std::vector<double>({0.0}), std::vector<double>({coalescence_rate}))
  {}

  /* use importance sampling and variable coalescence rate */
  EstimateBranchLengthsVariational (Data* const data, const std::vector<double>& epoch_start, const std::vector<double>& epoch_rate, const std::size_t importance_samples, const std::size_t random_seed)
    : data (data)
    , rng (random_seed)
    , importance_samples (importance_samples)
    , timescale(epoch_start, epoch_rate)
  {}


  /* ------- MAIN ALGORITHM ------- */

  void EP (Tree& tree) {
    /* extract parent-child trios */
    std::vector<trio_t> trios; trios.reserve(tree.nodes.size());
    for (auto& i : tree.nodes) {
      if (i.child_left != NULL and i.child_right != NULL) {
        trios.emplace_back(i.label, i.child_left->label, i.child_right->label);
      }
    }
  
    /* order trios by leaves subtended by parent */
    std::vector<Leaves> leaves; tree.FindAllLeaves(leaves);
    std::sort(trios.begin(), trios.end(), [&leaves](const trio_t& x, const trio_t& y){
        return leaves[x.parent].num_leaves < leaves[y.parent].num_leaves;
    });
    
    //for (auto& i : trios) { std::cout << "trios.emplace_back(" << i.parent << "," << i.lchild << "," << i.rchild << ");" << std::endl; }//DEBUG

    /* convert mutational counts and areas to gamma natural parameters */
    std::vector<gamma_distr_t> likelihoods;
    for (auto& i : tree.nodes) {
      likelihoods.emplace_back(i.num_events, -mutational_target(i));
    }

    //for (auto& i : likelihoods) { std::cout << "mutation_spans.emplace_back(" << i.alpha << "," << i.beta << ");" << std::endl; }//DEBUG

    /* calculate topology-conditioned coalescent prior */
    conditional_coalescent_prior_t prior (trios);
  
    /* date via expectation propagation */
    std::vector<double> node_age; node_age.reserve(tree.nodes.size());
    if (order) {
      /* use quadrature */
      expectation_propagation_t<std::mt19937> ep (trios, likelihoods, prior.edges, timescale, order);
      for (auto& i : ep.posterior) node_age.push_back(i.mean());
    } else {
      /* use importance sampling */
      assert (importance_samples);
      expectation_propagation_t<std::mt19937> ep (trios, likelihoods, prior.edges, timescale, importance_samples, rng);
      for (auto& i : ep.posterior) node_age.push_back(i.mean());
    }
  
    /* copy posterior mean into tree */
    for (auto& i : tree.nodes) {
      if (i.parent != NULL) {
        i.branch_length = node_age[i.parent->label] - node_age[i.label];
      }
    }
    // TODO: what about storing posterior uncertainty? variance? std dev?
    // this would be best done internally in EP, by creating gamma projection
    // for each branch during IS/quadrature
  }
};

#endif
