set -e
#grep "TIMING" sim_constne/*/benchmarks/benchmarks >sim_constne/timings
#grep "DATING" sim_constne/*/benchmarks/benchmarks >sim_constne/datings
awk '$1 ~ "DATING" {print $10" "$14}' sim_constne/*/benchmarks/benchmarks >sim_constne/datings_100

echo '
import numpy as np
import matplotlib.pyplot as plt
hi = np.loadtxt("sim_constne/datings_100")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5, 2))
fig.suptitle("EP vs MCMC\n(200 trees, 1000 diploids, 100 MCMC samples)")
ax1.hexbin(hi[:,0], hi[:,1], gridsize=50, mincnt=1, bins="log")
ax1.axline((0,0), slope=1, c="black", linestyle="dashed")
ax1.set_xlabel(f"E[age] (MCMC)")
ax1.set_ylabel(f"E[age] (EP)")
ax2.hexbin(np.log10(hi[:,0]), np.log10(hi[:,1]), gridsize=50, mincnt=1, bins="log")
ax2.axline((0,0), slope=1, c="black", linestyle="dashed")
ax2.set_xlabel(f"E[log10 age] (MCMC)")
ax2.set_ylabel(f"E[log10 age] (EP)")
fig.savefig("/home/natep/public_html/tsdate-paper/EP_vs_MCMC.png")
fig.clf()
' | python
