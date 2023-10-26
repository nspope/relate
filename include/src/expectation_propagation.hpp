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

#include "hypergeometric.hpp"
#include "quadrature.hpp"
#include "anc.hpp"
#include "data.hpp"

#define OPTIM_MAXITT 100
#define OPTIM_RELTOL 1e-8
#define IS_MINGAP 1e-12

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

  private:
    const std::size_t dim;
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
        u = betainc_t::lentz(b.shape(), a.shape() + 1, 1.0 - x);
        v = betainc_t::lentz(a.shape(), b.shape() + 1, x);
        mn = a.mean() * u + b.mean() * v;
        u = betainc_t::lentz(b.shape(), a.shape() + 2, 1.0 - x);
        v = betainc_t::lentz(a.shape(), b.shape() + 2, x);
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

    /* --- QUADRATURE SCHEME --- */

    static
    double quadrature_conditional(
        gamma_distr_t& p_i, // posterior (output)
        const gamma_distr_t& q_i, // cavity (input)
        const gamma_distr_t& prior, const gamma_distr_t& l_j, const gamma_distr_t& l_k, // likelihood and prior
        const timescale_t& timescale, const std::size_t order
    ) {
      /* Estimate moments of surrogate via Gaussian quadrature, conditioned on
         t_i > t_j = t_k = 0. The integral (omitting the prior) is,

         \int_0^\infty 
           t_i^{a_i + a_ij + a_ik} e^{t_i (b_i + b_ij + b_ik)}
         dt_i
       */

      /* prior on difference in ages between parent and oldest child,
         on the coalescent timescale */
      auto prior_i = [&prior,&timescale] (const double t_i) {
        return prior.loglik(timescale.to_coalescents(t_i));
      };

      /* approximate prior assuming constant population size */
      gamma_distr_t proposal(prior.alpha, prior.beta / timescale.pairwise_tmrca);
      auto prop_i = [&proposal,&timescale] (const double t_i) {
        return proposal.loglik(t_i);
      };

      /* quadrature over ages of parent and older child */
      double t_i, lp, alpha, beta, scale;
      const std::size_t dim = order;
      std::vector<double> weight, age_i, log_i;
      weight.reserve(dim); age_i.reserve(dim); log_i.reserve(dim);

      alpha = q_i.alpha + l_k.alpha + l_j.alpha + proposal.alpha;
      beta = q_i.beta + l_k.beta + l_j.beta + proposal.beta;
      assert (alpha > -1.0);
      assert (beta < 0.0);
      scale = log(-beta) * (-alpha - 1.0);
      gauss_laguerre_t quad_i (order, alpha);
      for (auto& i : quad_i.grid) {
        t_i = -i.abscissa / beta;
        lp = prior_i(t_i) - prop_i(t_i) + i.logweight + scale;
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
         t_i > t_j and t_k = 0. The integral (omitting prior) is,

         \int_0^\infty \int_0^t_i
           t_i^{a_i} t_j^{a_j} (t_i - t_j)^{a_ij} t_i^{a_ik} \times
           e^{b_i t_i + b_j t_j + (t_i - t_j) b_ij + t_i b_ik}
         dt_j dt_i

         (v := t_j/t_i)

         = \int_0^\infty \int_0^1
           t_i^{a_i + 1} (t_i v)^{a_j} (t_i - t_i v)^{a_ij} t_i^{a_ik} \times
           e^{b_i t_i + b_j t_i v + (t_i - t_i v) b_ij + t_i b_ik}
         dt_j dt_i

         = \int_0^1 v^{a_j} (1 - v)^{a_ij} \int_0^\infty
           t_i^{a_i + a_j + a_ij + a_ik + 1} \times
           e^{t_i (b_i + b_ik + v b_j + (1 - v) b_ij)}
         dt_i dt_j

       */

      /* prior on difference in ages between parent and oldest child,
         on the coalescent timescale */
      auto prior_ij = [&prior,&timescale] (const double t_i, const double t_j) {
        return prior.loglik(
            timescale.to_coalescents(t_i) - timescale.to_coalescents(t_j)
        );
      };

      /* approximate prior assuming constant population size */
      gamma_distr_t proposal(prior.alpha, prior.beta / timescale.pairwise_tmrca);
      auto prop_ij = [&proposal,&timescale] (const double t_i, const double t_j) {
        return proposal.loglik(t_i - t_j);
      };

      /* quadrature over ages of parent and older child */
      double t_i, t_j, lp, alpha, beta, scale;
      const std::size_t dim = std::pow(order, 2);
      std::vector<double> weight, age_i, log_i, age_j, log_j;
      weight.reserve(dim); age_i.reserve(dim); log_i.reserve(dim);
      age_j.reserve(dim); log_j.reserve(dim);

      alpha = q_i.alpha + q_j.alpha + l_k.alpha + l_j.alpha + proposal.alpha + 1.0;
      gauss_laguerre_t quad_i (order, alpha);
      gauss_jacobi_t quad_j (order, l_j.alpha + proposal.alpha, q_j.alpha);
      for (auto& j : quad_j.grid) {
        beta = q_i.beta + l_k.beta + (l_j.beta + proposal.beta)*(1 - j.abscissa) + q_j.beta*j.abscissa;
        assert (beta < 0.0);
        scale = log(-beta) * (-alpha - 1.0);
        for (auto& i : quad_i.grid) {
          t_i = -i.abscissa / beta;
          t_j = j.abscissa * t_i;
          /* collect integrand and samples */
          lp = prior_ij(t_i, t_j) - prop_ij(t_i, t_j) + i.logweight + j.logweight + scale;
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
         t_i > t_j > t_k. The integral (omitting prior) is,

         \int_0^\infty \int_0^t_i 
         t_i^{a_i} t_j^{a_j} (t_i - t_j)^{a_ij} \times
         e^{b_i t_i + b_j t_j + (t_i - t_j) b_ij + t_i b_ik} \times
         \int 0^t_j t_k^{a_k} (t_i - t_k)^{a_ik} e^{t_k (b_k - b_ik)} 
         dt_k dt_j dt_i

         (u := t_k/t_j)

         = \int_0^\infty \int_0^t_i 
         t_i^{a_i} t_j^{a_j} (t_i - t_j)^{a_ij} \times
         e^{b_i t_i + b_j t_j + (t_i - t_j) b_ij + t_i b_ik} \times
         t_j^{1 + a_k} t_i^{a_ik} \int 0^1 u^{a_k} (1 - t_j/t_i u)^{a_ik} e^{t_j u (b_k - b_ik)} 
         du dt_j dt_i

         (v := t_j/t_i)

         = \int_0^\infty \int_0^1
         t_i^{a_i + 1} (t_i v)^{a_j} (t_i - t_i v)^{a_ij} \times
         e^{b_i t_i + b_j t_i v + (t_i - t_i v) b_ij + t_i b_ik} \times
         (t_i v)^{1 + a_k} t_i^{a_ik} \int 0^1 u^{a_k} (1 - v u)^{a_ik} e^{t_i v u (b_k - b_ik)} 
         du dv dt_i

         = \int_0^\infty \int_0^1
         t_i^{a_i + a_j + a_k + a_ij + a_ik + 2} v^{a_j + a_k + 1} (1 - v)^{a_ij} \times
         e^{b_i t_i + b_j t_i v + (t_i - t_i v) b_ij + t_i b_ik + t_i v u (b_k - b_ik)} \times
         \int 0^1 u^{a_k} (1 - v u)^{a_ik}
         du dv dt_i

         Rearranging,

         = \int 0^1 v^{a_j + a_k + 1} (1 - v)^{a_ij} \int_0^1 u^{a_k} (1 - v u)^{a_ik} 
         \int_0^\infty t_i^{a_i + a_j + a_k + a_ij + a_ik + 2}
         e^{t_i (b_i + b_j v + (1 - v) b_ij + b_ik + v u (b_k - b_ik))}
         dt_i du dv
       */

      /* exact prior on difference in ages between parent and oldest child,
         on the coalescent timescale */
      auto prior_ij = [&prior,&timescale] (const double t_i, const double t_j) {
        return prior.loglik(timescale.to_coalescents(t_i) - timescale.to_coalescents(t_j));
      };

      /* approximate prior assuming constant population size */
      gamma_distr_t proposal(prior.alpha, prior.beta / timescale.pairwise_tmrca);
      auto prop_ij = [&proposal,&timescale] (const double t_i, const double t_j) {
        return proposal.loglik(t_i - t_j);
      };

      /* quadrature over ages of parent and older child */
      double t_i, t_j, t_k, ln_t_k, lp, alpha, beta, scale;
      const std::size_t dim = std::pow(order, 3);
      std::vector<double> weight, age_i, log_i, age_j, log_j, age_k, log_k; 
      weight.reserve(dim); age_i.reserve(dim); log_i.reserve(dim);
      age_j.reserve(dim); log_j.reserve(dim); age_k.reserve(dim); log_k.reserve(dim);

      alpha = q_i.alpha + q_k.alpha + q_j.alpha + l_k.alpha + l_j.alpha + proposal.alpha + 2;
      gauss_laguerre_t quad_i (order, alpha);
      gauss_jacobi_t quad_j (order, l_j.alpha + proposal.alpha, q_j.alpha + q_k.alpha + 1);
      gauss_jacobi_t quad_k (order, 0.0, q_k.alpha);
      for (auto& j : quad_j.grid) {
        for (auto& k : quad_k.grid) {
          beta = q_i.beta + q_j.beta*j.abscissa + 
            (l_j.beta + proposal.beta)*(1 - j.abscissa) + 
            q_k.beta*j.abscissa*k.abscissa + l_k.beta*(1 - j.abscissa*k.abscissa);
          assert (beta < 0.0);
          scale = log(-beta) * (-alpha - 1.0); 
          for (auto& i : quad_i.grid) {
            t_i = -i.abscissa / beta;
            t_j = j.abscissa * t_i;
            t_k = k.abscissa * t_j;
            /* collect integrand and samples */
            lp = prior_ij(t_i, t_j) - prop_ij(t_i, t_j) + 
              log(1.0 - j.abscissa * k.abscissa)*l_k.alpha + 
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

      double logconst;

      if (q_j && q_k) { /* ages of both children not fixed */
        gamma_distr_t p_i_jk, p_j_jk, p_k_jk, p_i_kj, p_j_kj, p_k_kj;
        double lp_jk, lp_kj, w_jk, w_kj, maxlp, mn_i, ln_i, mn_j, ln_j, mn_k, ln_k;
        /* conditioning on t_i > t_j > t_k */
        lp_jk = quadrature_conditional(p_i_jk, p_j_jk, p_k_jk, q_i, q_j, q_k, prior, l_j, l_k, timescale, order);
        /* conditioning on t_i > t_k > t_j */
        lp_kj = quadrature_conditional(p_i_kj, p_k_kj, p_j_kj, q_i, q_k, q_j, prior, l_k, l_j, timescale, order);
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
        logconst = quadrature_conditional(p_i, p_j, q_i, q_j, prior, l_j, l_k, timescale, order);
      } else if (q_k) {
        /* conditioning on t_i > t_k > t_j = 0 */
        logconst = quadrature_conditional(p_i, p_k, q_i, q_k, prior, l_k, l_j, timescale, order);
      } else {
        /* conditioning on t_i > t_k = t_j = 0 */
        logconst = quadrature_conditional(p_i, q_i, prior, l_k, l_j, timescale, order);
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

      /* approximate prior with constant effective population size */
      gamma_distr_t proposal (prior.alpha, prior.beta / timescale.pairwise_tmrca);

      auto rng_i = prior.random(), rng_j = q_j.random(), rng_k = q_k.random();
      std::vector<double> t_i(n), t_j(n), t_k(n), c_i(n), w(n), lp(n);
      for (std::size_t itt = 0; itt < n; ++itt) {
        /* child ages sampled from variational cavity */
        t_j[itt] = rng_j(rng); 
        t_k[itt] = rng_k(rng);

        /* branch length from eldest child to parent sampled from coalescent */
        t_i[itt] = std::max(t_j[itt], t_k[itt]);
        c_i[itt] = std::max(rng_i(rng), IS_MINGAP); // minimum age gap
        t_i[itt] = timescale.to_coalescents(t_i[itt]) + c_i[itt];
        t_i[itt] = timescale.to_generations(t_i[itt]);

        /* likelihood weights */
        w[itt] = q_i.loglik(t_i[itt]); // variational prior for parent
        w[itt] += l_j.loglik(t_i[itt] - t_j[itt]);
        w[itt] += l_k.loglik(t_i[itt] - t_k[itt]); 

        /* normalizing constant */
        lp[itt] = w[itt] + q_j.loglik(t_j[itt]) + q_k.loglik(t_k[itt]) + prior.loglik(c_i[itt]);
      }

      /* importance sampling */
      double maxlp, maxw, norm, logconst;
      maxlp = *std::max_element(lp.begin(), lp.end());
      logconst = std::accumulate(lp.begin(), lp.end(), 0.0, 
          [&maxlp](double x, double y) { return x + exp(y - maxlp);}
      );
      logconst = log(logconst) + maxlp; // TODO is this correct?
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


  public:
    std::vector<gamma_distr_t> posterior;

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

// -------------- RELATE API ------------------ //

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
  {}

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
    // TODO: this assumes that node.label are in [0 ... num_samples*2 - 1] ...
    // is this always true?
  
    /* order trios by leaves subtended by parent */
    std::vector<Leaves> leaves; tree.FindAllLeaves(leaves);
    std::sort(trios.begin(), trios.end(), [&leaves](const trio_t& x, const trio_t& y){
        return leaves[x.parent].num_leaves < leaves[y.parent].num_leaves;
    });
    
    for (auto& i : trios) { std::cout << "trios.emplace_back(" << i.parent << "," << i.lchild << "," << i.rchild << ");" << std::endl; }//DEBUG

    /* convert mutational counts and areas to gamma natural parameters */
    std::vector<gamma_distr_t> likelihoods;
    for (auto& i : tree.nodes) {
      likelihoods.emplace_back(i.num_events, -mutational_target(i));
    }

    for (auto& i : likelihoods) { std::cout << "mutation_spans.emplace_back(" << i.alpha << "," << i.beta << ");" << std::endl; }//DEBUG

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
