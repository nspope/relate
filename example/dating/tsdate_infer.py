import tsdate
import tskit
import numpy as np

# read in zigazg
#exec(open("zigzag.py").read())

# load topologies
mu = 1.29e-8
ts = tskit.load("sim_varne2/2000/chr1.infer.trees")

# setup demographic model
#demog = tsdate.demography.PopulationSizeHistory(ne, br)
ne_est = 0.25 * ts.diversity() / mu

# setup prior
prior_class = tsdate.prior.MixturePrior(ts, prior_distribution="gamma", progress=True)
prior = prior_class.make_parameter_grid(ne_est) #prior_class.make_parameter_grid(demog)
prior_mixture = prior_class.to_gamma_mixture(5, population_size=ne_est)
global_prior = tsdate.approx.average_gammas(prior.grid_data[:,0], prior.grid_data[:,1])
prior.grid_data[:] = [1, 0]

#ts_org = tsdate.date(ts, mutation_rate=mu, population_size=ne_est, progress=True)
#ts_org.dump(f"sim_varne2/chr1.infer.tsorig.trees")
ts_mix = tsdate.date(ts, mutation_rate=mu, priors=prior, global_prior=False, global_distribution=prior_mixture, method="variational_gamma", max_iterations=50, progress=True, max_shape=100)
ts_mix.dump(f"sim_varne2/chr1.infer.tsdate.trees")

