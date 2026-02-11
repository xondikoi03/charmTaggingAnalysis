import matplotlib.pyplot as plt 
import dask
import awkward as ak
import hist.dask as dhist
import hist as h
import mplhep as mh
import dask_awkward as dak
import numpy as np

# Defining the style for plots (CMS Style)
from mplhep.plot import AnchoredText
plt.style.use(mh.style.CMS)
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12
})

# Histogram Function:
def make_hist(
        data, nBins, lo, hi, xlabel, ylabel, label, fname=None, logy=False, ax=None, show=True, return_hist=False, show_stats=True, histtype='step'
):
    histogram = dhist.Hist(h.axis.Regular(nBins, lo, hi))
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(7, 5))
    else:
        fig = ax.figure
    histogram.fill(data)
    hist_comp = histogram.compute()

    mh.histplot(hist_comp, ax=ax, label=label, histtype=histtype, lw=1)
    mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2016, loc=0)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()

    # log scale option:
    if logy:
        ax.set_yscale("log")

    # Statbox option:
    if show_stats:
        n_entries = len(data.compute())
        mean_val = np.mean(data.compute())
        std_val = np.std(data.compute())
        stat_text = (f"Entries: {n_entries}\n"
                     f"Mean:    {mean_val:.2f}\n"
                     f"Std Dev: {std_val:.2f}")
        at = AnchoredText(stat_text, prop=dict(size=10), frameon=True, loc="center right")
        ax.add_artist(at)

    # Save figure option:
    if fname:
        plt.savefig(fname)

    # Show figure option:
    if show:
        plt.show()

    if return_hist:
        return fig, ax, hist_comp
    return fig, ax


# Making Difference Histogram Function:
def make_diff_hist(
    a, b, *, bins=50, lo=0, hi=200, ylo=-1, yhi=2000, name="x",
    xlabel=None, ylabel="Entries", title=None, diff_label="A - B", color="k", histtype="errorbar", yerr=False, zero_line=True, normalize=False, weights_a=None,weights_b=None, flatten=False, cms_label=True, ax=None,
    figsize=(7, 5), grid=True, show=True, fname=None,
    return_hists=False,
):
    # Detect hist input
    try:
        from hist.basehist import BaseHist
        a_is_hist = isinstance(a, BaseHist)
        b_is_hist = isinstance(b, BaseHist)
    except Exception:
        a_is_hist = False
        b_is_hist = False

    def _prep(x):
        if flatten:
            try:
                return dak.flatten(x)
            except Exception:
                return x
        return x

    # Build/obtain hA and hB
    if a_is_hist and b_is_hist:
        hA = a
        hB = b
    elif (not a_is_hist) and (not b_is_hist):
        axis = h.axis.Regular(
            bins, lo, hi,
            name=name,
            label=(xlabel if xlabel is not None else name),
        )
        hA_d = h.dask.Hist(axis)
        hB_d = h.dask.Hist(axis)

        xA = _prep(a)
        xB = _prep(b)

        if weights_a is None:
            hA_d.fill(**{name: xA})
        else:
            hA_d.fill(**{name: xA}, weight=_prep(weights_a))

        if weights_b is None:
            hB_d.fill(**{name: xB})
        else:
            hB_d.fill(**{name: xB}, weight=_prep(weights_b))

        hA = hA_d.compute()
        hB = hB_d.compute()
    else:
        raise TypeError("make_diff_hist: both inputs must be arrays OR both must be hist objects.")

    if normalize:
        sA = float(hA.values().sum())
        sB = float(hB.values().sum())
        if sA > 0:
            hA = hA / sA
        if sB > 0:
            hB = hB / sB

    hDiff = hA - hB

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    if cms_label:
        mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2016, loc=0)

    mh.histplot(hDiff, ax=ax, label=diff_label, histtype=histtype, color=color, yerr=yerr)

    if zero_line:
        ax.axhline(0.0, color="0.5", lw=1, linestyle="--")

    if title:
        ax.set_title(title)

    ax.set_xlabel(xlabel if xlabel is not None else name)
    ax.set_ylabel(ylabel)
    ax.set_ylim(ylo, yhi)
    ax.legend()
    if grid:
        ax.grid(True, alpha=0.3)
    if fname:
        plt.savefig(fname)
    if show:
        plt.show()

    if return_hists:
        return fig, ax, hDiff, hA, hB
    return fig, ax, hDiff


# Overlay Histogram Function:
def make_overlay_hist(
    data_a, data_b, *, bins=50, x_range=(0, 200), name="x",
    histtype="step", xlabel=None, ylabel=None, title=None,
    normalize=True, labels=("A", "B"), colors=("C0", "C1"),
    weights_a=None, weights_b=None, flatten=False, cms_label=True, ax=None, figsize=(7, 5), grid=True,
    fname=None, show=True, yerr=False, linewidth=1.5,     
):
    def _prep(x):
        if flatten:
            try:
                return dak.flatten(x)
            except Exception:
                return x
        return x

    axis = h.axis.Regular(
        bins, x_range[0], x_range[1],
        name=name,
        label=(xlabel if xlabel is not None else name),
    )

    hA_d = h.dask.Hist(axis)
    hB_d = h.dask.Hist(axis)

    xA = _prep(data_a)
    xB = _prep(data_b)

    if weights_a is None:
        hA_d.fill(**{name: xA})
    else:
        hA_d.fill(**{name: xA}, weight=_prep(weights_a))

    if weights_b is None:
        hB_d.fill(**{name: xB})
    else:
        hB_d.fill(**{name: xB}, weight=_prep(weights_b))

    hA = hA_d.compute()
    hB = hB_d.compute()

    if normalize:
        sA = float(hA.values().sum())
        sB = float(hB.values().sum())
        if sA > 0:
            hA = hA / sA
        if sB > 0:
            hB = hB / sB

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    mh.histplot(
        hA,
        ax=ax,
        label=labels[0],
        color=colors[0],
        histtype=histtype,
        yerr=yerr,
        linewidth=linewidth if histtype == "step" else None,
    )
    mh.histplot(
        hB,
        ax=ax,
        label=labels[1],
        color=colors[1],
        histtype=histtype,
        yerr=yerr,
        linewidth=linewidth if histtype == "step" else None,
    )

    if cms_label:
        mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2016, loc=0)

    if title:
        ax.set_title(title)

    ax.set_xlabel(xlabel if xlabel is not None else name)
    ax.set_ylabel(
        ylabel if ylabel is not None else ("Normalized entries" if normalize else "Entries")
    )
    ax.legend()
    if grid:
        ax.grid(True, alpha=0.3)
    if fname:
        plt.savefig(fname)
    if show:
        plt.show()

    return fig, ax, hA, hB

# OS-SS histogram function:

def os_ss_hist(x_os, x_ss, *, bins=20, lo=20, hi=200, name="pt", label=None, normalize=False):
    axis = h.axis.Regular(bins, lo, hi, name=name, label=(label if label else name))
    hOS_d = h.dask.Hist(axis)
    hSS_d = h.dask.Hist(axis)
    hOS_d.fill(**{name: x_os})
    hSS_d.fill(**{name: x_ss})
    hDiff = hOS_d.compute() - hSS_d.compute()
    if normalize:
        s = float(hDiff.values().sum())
        if s != 0:
            hDiff = hDiff / s
    return hDiff