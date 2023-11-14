import tsdate
import tskit
import numpy as np
import scipy.integrate
from math import exp, log10
import stdpopsim

import matplotlib.pyplot as plt
import msprime

homsap = stdpopsim.get_species("HomSap")
zigzag = homsap.get_demographic_model("Zigzag_1S14")
dbg = zigzag.model.debug()
epoch_start = np.logspace(1, log10(40000), 1000)
epoch_rate, _ = dbg.coalescence_rate_trajectory(lineages={"generic":2}, steps=epoch_start)

#plt.step(np.log10(np.array(epoch_start)), np.log10(2./np.array(epoch_rate)))
plt.step(np.array(epoch_start), np.log10(2./np.array(epoch_rate)))
plt.title("\"Zigzag\" model")
plt.xlabel("Generations ago")
plt.ylabel("Eff. population size (log10)")
plt.xscale('log')
plt.savefig("/home/natep/public_html/tsdate-paper/zigzag2.png")
plt.clf()


ts = tskit.load("sim_varne2/2000/chr1.trees")

def plot_loghist(x, bins):
  hist, bins = np.histogram(x[x > 0], bins=bins)
  logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
  plt.hist(x, bins=logbins)
  plt.xscale('log')
  plt.xlabel("Node age")
  plt.ylabel("Density")
  plt.savefig("/home/natep/public_html/tsdate-paper/zigzag_hist.png")

plot_loghist(ts.nodes_time, 30)
