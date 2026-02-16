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
###############################
# Histogram Function:
###############################
def make_hist(
        data, nBins, lo, hi, xlabel, ylabel, label, fname=None, logy=False, ax=None, show=True, return_hist=False, show_stats=True, histtype='step', title=None
):
    """
    Build and plot a 1D histogram from a (possibly dask-awkward) array.

    Parameters:
        data: Array-like (awkward/dask-awkward) values to histogram.
        nBins, lo, hi: Histogram binning configuration.
        xlabel, ylabel: Axis labels for the plot.
        label: Legend label for the histogram.
        fname: Optional path to save the figure.
        logy: Use log scale on the y-axis.
        ax: Optional matplotlib Axes to plot into.
        show: Display the figure.
        return_hist: Return computed hist object in addition to fig/ax.
        show_stats: Add a simple stat box with entries/mean/std.
        histtype: Passed to mplhep histplot (e.g. "step", "fill").
        title: Optional plot title.

    Returns:
        (fig, ax) or (fig, ax, hist) when return_hist is True.
    """
    histogram = dhist.Hist(h.axis.Regular(nBins, lo, hi))
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(8, 6))
    else:
        fig = ax.figure
    histogram.fill(data)
    hist_comp = histogram.compute()

    mh.histplot(hist_comp, ax=ax, label=label, histtype=histtype, lw=2, color="C1")
    mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2016, loc=0)
    if title:
        ax.set_title(title)
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

########################################
# Making Difference Histogram Function:
########################################
def make_diff_hist(
    a, b, *, bins=50, lo=0, hi=200, ylo=-1, yhi=2000, name="x",
    xlabel=None, ylabel="Entries", title=None, diff_label="A - B", color="k", histtype="errorbar", yerr=False, zero_line=True, normalize=False, weights_a=None,weights_b=None, flatten=False, cms_label=True, ax=None,
    figsize=(8, 6), grid=True, show=True, fname=None,
    return_hists=False, return_errors=False,
):
    """
    Plot the difference of two histograms (A - B).

    Accepts either two arrays or two hist objects. If arrays are provided,
    histograms are built with the supplied binning and optional weights.

    Parameters:
        a, b: Input arrays or hist objects (both must be the same type). If not a hist, it converts it hist objects.
        bins, lo, hi: Histogram binning when inputs are arrays.
        ylo, yhi: y-axis limits.
        name: Axis name for histogram filling.
        xlabel, ylabel: Axis labels for the plot.
        title: Optional plot title.
        diff_label: Legend label for the difference histogram.
        color: Line/marker color for the difference plot.
        histtype: mplhep histplot style (e.g. "errorbar", "step").
        yerr: Plot y error bars when supported.
        zero_line: Draw a horizontal line at zero.
        normalize: Normalize both inputs before subtraction.
        weights_a, weights_b: Optional weights for array inputs.
        flatten: Flatten awkward arrays before filling.
        cms_label: Draw CMS label on the plot.
        ax: Optional matplotlib Axes to plot into.
        figsize: Figure size when creating a new figure.
        grid: Show grid lines.
        show: Display the figure.
        fname: Optional path to save the figure.
        return_hists: Return input hists along with the difference.
        return_errors: Return per-bin uncertainties for the difference.

    Returns:
        (fig, ax, hDiff), (fig, ax, hDiff, err), or
        (fig, ax, hDiff, err, hA, hB) when return_hists/return_errors are True.
    """

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
    err = None
    if yerr or return_errors:
        try:
            err = np.sqrt(hA.variances() + hB.variances())
        except Exception:
            err = None

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    if cms_label:
        mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2016, loc=0)

    plot_yerr = err if err is not None else yerr
    mh.histplot(hDiff, ax=ax, label=diff_label, histtype=histtype, color=color, yerr=plot_yerr)

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

    if return_hists and return_errors:
        return fig, ax, hDiff, err, hA, hB
    if return_hists:
        return fig, ax, hDiff, hA, hB
    if return_errors:
        return fig, ax, hDiff, err
    return fig, ax, hDiff

###############################
# Overlay Histogram Function:
###############################
def make_overlay_hist(
    data_a, data_b, *, bins=50, x_range=(0, 200), name="x",
    histtype="step", xlabel=None, ylabel=None, title=None,
    normalize=True, labels=("A", "B"), colors=("C0", "C1"),
    weights_a=None, weights_b=None, flatten=False, cms_label=True, ax=None, figsize=(8, 6), grid=True,
    fname=None, show=True, yerr=False, linewidth=1.5,     
):
    """
    Plot two histograms overlaid for comparison.

    Parameters:
        data_a, data_b: Input arrays to histogram and overlay.
        bins, x_range: Histogram binning when inputs are arrays.
        name: Axis name for histogram filling.
        histtype: mplhep histplot style (e.g. "step", "fill").
        xlabel, ylabel: Axis labels for the plot.
        title: Optional plot title.
        normalize: Normalize each histogram to unit area.
        labels: Legend labels for the two histograms.
        colors: Plot colors for the two histograms.
        weights_a, weights_b: Optional weights for inputs.
        flatten: Flatten awkward arrays before filling.
        cms_label: Draw CMS label on the plot.
        ax: Optional matplotlib Axes to plot into.
        figsize: Figure size when creating a new figure.
        grid: Show grid lines.
        fname: Optional path to save the figure.
        show: Display the figure.
        yerr: Plot y error bars when supported.
        linewidth: Line width when using step histtype.

    Returns:
        (fig, ax, hA, hB)
    """
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

###############################
# OS-SS histogram function:
###############################

def os_ss_hist(x_os, x_ss, *, bins=20, lo=20, hi=200, name="pt", label=None, normalize=False):
    """
    Build an OS-SS (opposite-sign minus same-sign) histogram.

    Parameters:
        x_os, x_ss: Input arrays for opposite-sign and same-sign samples.
        bins, lo, hi: Histogram binning configuration.
        name: Axis name for histogram filling.
        label: Axis label override (defaults to name).
        normalize: Normalize the OS-SS histogram to unit area.

    Returns:
        Computed hist object for OS-SS.
    """
    """
    
    """
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

