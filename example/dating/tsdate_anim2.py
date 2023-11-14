import tsdate
import tskit
import numpy as np
import os

# load topologies
mu = 1.29e-8
ts = tskit.load("sim_varne2/2000/chr1.trees")

# setup demographic model
ne_est = 0.25 * ts.diversity() / mu

# setup prior
prior_class = tsdate.prior.MixturePrior(ts, prior_distribution="gamma", progress=True)
prior = prior_class.make_parameter_grid(ne_est)
prior_mixture = prior_class.to_gamma_mixture(5, population_size=ne_est)
global_prior = tsdate.approx.average_gammas(prior.grid_data[:,0], prior.grid_data[:,1])
prior.grid_data[:] = [1, 0]

if not os.path.exists("anim2"): os.mkdir("anim2")

for i in range(1,50,2):
    ts_mix = tsdate.date(ts, mutation_rate=mu, priors=prior, global_prior=False, global_distribution=prior_mixture, method="variational_gamma", max_iterations=i, progress=True, max_shape=100)
    nm = format(i, "03")
    ts_mix.dump(f"anim2/chr1.tsdate.iter{nm}.trees")

