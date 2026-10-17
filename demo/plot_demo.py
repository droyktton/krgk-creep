"""Figures of the demo (demo/run_demo.sh). Usage: python3 plot_demo.py RUNS_DIR FIG_DIR

single_realization.png (one disorder realization, L = 32, M = 64):
  the metastable states of one period of the steady-state cycle, one creep event (alpha, saddle,
  gamma), the barriers along the cycle, and the dominant configuration of each drive compared
  with the ground state (F = 0) and the critical configuration (F = F_c).
statistics.png (L = 64, f = 0.6, three realizations, steady state):
  the structure factor of the metastable states and the statistics of the events.

At T -> 0+ the time spent in a metastable state grows as exp(U_B/T), so the steady-state cycle
is dominated by the state with the largest barrier: the dominant configuration.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(d):
    par = dict(l.split(None, 1) for l in open(f"{d}/params.txt"))
    L, M = int(par["L"]), int(par["M"])
    ev = np.loadtxt(f"{d}/events.dat", ndmin=2)
    # events.bin: per event alpha, saddle, beta, gamma (int32, L each)
    conf = np.fromfile(f"{d}/events.bin", dtype=np.int32)
    conf = conf[:len(conf) // (4 * L) * 4 * L].reshape(-1, 4, L)
    n = min(len(ev), len(conf))
    return par, L, M, ev[:n], conf[:n]


def cycle(gamma, M):
    """Period P (in events) of the steady state, gamma[n + P] = gamma[n] + M, from the end."""
    n = len(gamma)
    for P in range(1, n // 2 + 1):
        if np.all(gamma[n - P:] == gamma[n - 2 * P:n - P] + M):
            return P
    return None


def align(x, ref, M):
    """Shift x by a multiple of M (the period of the disorder) to lie closest to ref."""
    return x - M * np.round(np.mean(x - ref) / M)


def single(runs, figs):
    fs = sorted(float(d[len("single_f"):]) for d in os.listdir(runs) if d.startswith("single_f"))
    data = {}
    for f in fs:
        d = f"{runs}/single_f{f:g}"
        par, L, M, ev, conf = load(d)
        P = cycle(conf[:, 3], M)
        if P is None:
            print(f"f = {f}: no periodic cycle yet, run longer"); continue
        last = slice(len(ev) - P, len(ev))
        k = len(ev) - P + int(np.argmax(ev[last, 1]))  # event with the largest barrier
        data[f] = dict(L=L, M=M, ev=ev, conf=conf, P=P, last=last, k=k, Fc=float(par["Fc"]),
                       gs=np.loadtxt(f"{d}/ground_state.dat"), crit=np.loadtxt(f"{d}/critical.dat"))
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(data)))
    i = np.arange(next(iter(data.values()))["L"])

    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    # (a) metastable states along one period at an intermediate drive
    f0 = min(data, key=lambda f: abs(f - 0.4))
    D = data[f0]
    g = D["conf"][D["last"], 3]
    for j, x in enumerate(g):
        ax[0, 0].plot(i, x, color=plt.cm.Greys(0.3 + 0.6 * j / len(g)), lw=0.7)
    a, s, b, c = D["conf"][D["k"]]
    ax[0, 0].plot(i, a, "C0", lw=2.5, label=r"dominant state $\alpha$ (largest $U_B$)")
    ax[0, 0].plot(i, s, "C3--", lw=1.5, label="saddle (sets $U_B$)")
    ax[0, 0].plot(i, c, "C2", lw=2, label=r"next metastable state $\gamma$")
    ax[0, 0].set(xlabel="$i$", ylabel="$u_i$",
                 title=f"(a) metastable states over one period, f = {f0:g}\n({D['P']} events, advance M = {D['M']})")
    ax[0, 0].legend(fontsize=8, loc="upper left")

    # (b) barriers along the cycle
    for (f, D), col in zip(data.items(), cols):
        U = D["ev"][D["last"], 1]
        ax[0, 1].semilogy(np.arange(len(U)) / len(U), U, ".-", color=col, ms=3, lw=0.6, label=f"f = {f:g}")
    ax[0, 1].set(xlabel="event / period", ylabel="$U_B$",
                 title="(b) barriers along one period of the cycle\n(time per state $\\propto e^{U_B/T}$)")
    ax[0, 1].legend(fontsize=8)

    # (c) dominant configurations against the ground state and the critical configuration
    D = data[fs[0]]
    gs = D["gs"]
    ax[1, 0].plot(i, gs, "k", lw=3, label="ground state (F = 0)")
    ax[1, 0].plot(i, align(D["crit"], gs, D["M"]), "C3", lw=2.5, label="critical configuration ($F_c$)")
    for (f, D), col in zip(data.items(), cols):
        ax[1, 0].plot(i, align(D["conf"][D["k"], 0], gs, D["M"]), color=col, lw=1.2, label=f"dominant, f = {f:g}")
    ax[1, 0].set(xlabel="$i$", ylabel="$u_i$ (mod M, aligned)",
                 title="(c) dominant configuration vs ground and critical states")
    ax[1, 0].legend(fontsize=7, ncol=2)

    # (d) barrier of each state of the cycle against its distance to the ground state
    dist = lambda x, y, M: np.sqrt(np.mean((align(x, y, M) - y) ** 2))
    for (f, D), col in zip(data.items(), cols):
        r = np.arange(len(D["ev"]))[D["last"]]
        dg = [dist(D["conf"][k, 0], D["gs"], D["M"]) for k in r]
        ax[1, 1].semilogy(dg, D["ev"][r, 1], "o", color=col, ms=3, alpha=0.6, label=f"f = {f:g}")
        ax[1, 1].semilogy(dist(D["conf"][D["k"], 0], D["gs"], D["M"]), D["ev"][D["k"], 1], "*", color=col, ms=14, mec="k")
    ax[1, 1].set(xlabel="rms distance to the ground state", ylabel="$U_B$",
                 title="(d) barrier of each state of the cycle (star: dominant)")
    ax[1, 1].legend(fontsize=8)
    fig.suptitle(f"One disorder realization, L = {D['L']}, M = {D['M']}, $F_c$ = {D['Fc']:.4g}")
    fig.tight_layout()
    fig.savefig(f"{figs}/single_realization.png", dpi=120)
    for f, D in data.items():
        x = D["conf"][D["k"], 0]
        print(f"f = {f:g}: period {D['P']} events, dominant U_B = {D['ev'][D['k'], 1]:.3f}, rms distance "
              f"to ground state {dist(x, D['gs'], D['M']):.2f}, to critical {dist(x, D['crit'], D['M']):.2f}")


def logbins(x, n=20, integer=False):
    x = x[x > 0]
    e = np.logspace(np.log10(x.min()), np.log10(x.max() + 1), n)
    if integer:  # edges at half-integers, so that each bin holds whole values
        e = np.unique(np.round(e)) - 0.5
    h, e = np.histogram(x, e, density=True)
    c = np.sqrt(e[1:] * e[:-1])
    return c[h > 0], h[h > 0]


def statistics(runs, figs, skip=1000):
    evs, S = [], []
    for d in sorted(x for x in os.listdir(runs) if x.startswith("stats_s")):
        par, L, M, ev, conf = load(f"{runs}/{d}")
        ss = ev[:, 5] - ev[0, 5] > skip  # steady state: discard the first advance `skip`
        evs.append(ev[ss])
        g = conf[ss, 3].astype(float)
        uq = np.fft.rfft(g - g.mean(1, keepdims=True), axis=1)
        S.append(np.abs(uq[:, 1:]) ** 2 / L)
        f = float(par["F"]) / float(par["Fc"])
    ev, S = np.vstack(evs), np.vstack(S).mean(0)
    q = 2 * np.pi * np.arange(1, L // 2 + 1) / L
    Ub, Lopt, Lrel, area = ev[:, 1], ev[:, 2], ev[:, 3], ev[:, 4]

    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    ax[0, 0].loglog(q, S, "ko", ms=4, label="metastable states")
    for z, ls in ((1.25, "--"), (2 / 3, ":")):
        ax[0, 0].loglog(q, S[1] * (q / q[1]) ** (-1 - 2 * z), "C3" + ls, label=f"$q^{{-(1+2\\zeta)}}$, $\\zeta$ = {z:.3g}")
    ax[0, 0].set(xlabel="$q$", ylabel="$S(q)$", title="(a) structure factor")
    ax[0, 0].legend(fontsize=8)
    for a, x, lab, integer in ((ax[0, 1], Lrel, "$L_{rel}$ (sites moved)", True), (ax[0, 2], area, "area", True),
                               (ax[1, 0], Ub, "barrier $U_B$", False)):
        c, h = logbins(x, integer=integer)
        a.loglog(c, h, "o-", ms=4)
        a.set(xlabel=lab, ylabel="probability density")
    ax[0, 1].set_title("(b) event length")
    ax[0, 2].set_title("(c) event area")
    ax[1, 0].set_title("(d) barrier")
    ax[1, 1].loglog(Lrel, area, ".", ms=2, alpha=0.3)
    big = Lrel >= 8
    b, a0 = np.polyfit(np.log(Lrel[big]), np.log(area[big]), 1)
    ll = np.array([2, L])
    ax[1, 1].loglog(ll, np.exp(a0) * ll ** b, "C3--", label=f"fit $L_{{rel}} \\geq 8$: area $\\sim L^{{{b:.2f}}}$")
    ax[1, 1].set(xlabel="$L_{rel}$", ylabel="area", title="(e) area vs length")
    ax[1, 1].legend(fontsize=8)
    ax[1, 2].loglog(Lopt * np.exp(0.08 * np.random.default_rng(0).standard_normal(len(Lopt))), Ub, ".", ms=2, alpha=0.2)
    lo = np.unique(Lopt)
    lo = lo[[np.sum(Lopt == l) >= 20 for l in lo]]
    ax[1, 2].loglog(lo, [np.median(Ub[Lopt == l]) for l in lo], "C3o-", label="median")
    ax[1, 2].legend(fontsize=8)
    ax[1, 2].set(xlabel="$L_{opt}$ (sites of the thermally activated excitation)", ylabel="$U_B$",
                 title="(f) barrier vs size of the excitation")
    fig.suptitle(f"Steady-state creep, L = {L}, f = {f:.2g}, {len(evs)} realizations, {len(ev)} events")
    fig.tight_layout()
    fig.savefig(f"{figs}/statistics.png", dpi=120)
    print(f"statistics: {len(ev)} events; mean L_rel {Lrel.mean():.1f}, mean area {area.mean():.1f}, "
          f"mean U_B {Ub.mean():.3f}")


if __name__ == "__main__":
    runs, figs = sys.argv[1], sys.argv[2]
    os.makedirs(figs, exist_ok=True)
    single(runs, figs)
    statistics(runs, figs)
