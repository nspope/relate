import tskit
import matplotlib.pyplot as plt
import numpy as np

def RMSE(x, y):
    ok_x = np.logical_and(x > 0, np.isfinite(x))
    ok_y = np.logical_and(y > 0, np.isfinite(y))
    ok = np.logical_and(ok_x, ok_y)
    return np.mean(np.abs(x[ok] - y[ok]))

def extract_mutations(ts, offset):
    muts = {}
    bad = set()
    for m in ts.mutations():
        p = ts.sites_position[m.site] + offset
        if p in muts:
            bad.add(p)
        else:
            edge = m.edge
            parent = ts.edges_parent[edge]
            child = ts.edges_child[edge]
            muts[p] = (ts.nodes_time[parent] + ts.nodes_time[child])/2.0
    for p in bad:
        del muts[p]
    return muts

def plot_mutation_ages(ts1, ts2, title, xlabel, ylabel, outfile, offset=[0,0], do_log=True, xlim=None, ylim=None):
    muts1 = extract_mutations(ts1, offset[0])
    muts2 = extract_mutations(ts2, offset[1])
    pos1 = set(muts1.keys())
    pos2 = set(muts2.keys())
    pos = pos1.intersection(pos2)
    x = []
    y = []
    for p in pos:
        x.append(muts1[p])
        y.append(muts2[p])
    x = np.array(x)
    y = np.array(y)
    if do_log:
        x = np.log10(x)
        y = np.log10(y)
    rmse = round(RMSE(x,y), 3)
    fig, (ax1) = plt.subplots(1, 1, figsize=(6, 5))
    fig.suptitle(title)
    ax1.hexbin(x, y, gridsize=50, mincnt=1, bins='log')
    ax1.axline((0,0), slope=1, c="black", linestyle="dashed")
    plt.text(.01, .99, f'RMSE: {rmse}', ha='left', va='top', transform=ax1.transAxes)
    ax1.set_xlabel(f"{xlabel}")
    ax1.set_ylabel(f"{ylabel}")
    if xlim is not None:
        ax1.set_xlim(left=xlim[0], right=xlim[1])
    if ylim is not None:
        ax1.set_ylim(bottom=ylim[0], top=ylim[1])
    fig.savefig(outfile)
    fig.clf()

ts1 = tskit.load("sim_varne2/2000/chr1.trees")
for i in range(1,39,2):
    nm = format(i, "03")
    ts2 = tskit.load(f"anim2/chr1.tsdate.iter{nm}.trees")
    plot_mutation_ages(ts1, ts2, f"Expectation propagation\n(iteration {i})", "True mutation age (log10)", "Inferred mutation age (log10)", outfile=f"anim2/{nm}.png", xlim=[-1, 6], ylim=[-1,6])


