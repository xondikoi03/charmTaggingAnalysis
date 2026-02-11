import os
import hist as h
import dask
import awkward as ak
import hist.dask as dhist
import dask_awkward as dak
import uproot
import numpy as np
import mplhep as mh
from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

# =================================================================
# 1. SETUP
# =================================================================
OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"plots will be saved to: {os.path.abspath(OUTPUT_DIR)}/")

NanoAODSchema.warn_missing_crossrefs = False 
plt.style.use(mh.style.CMS)
plt.rcParams.update({
    "font.size": 12, "axes.labelsize": 12, "axes.titlesize": 12, 
    "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12
})

def make_hist(data, nBins, lo, hi, xLabel, yLabel, label, filename, logy=False, color=None):
    """
    Creates, plots, and saves a histogram to the output directory.
    """
    if isinstance(data, dak.Array) and data.ndim > 1:
        data_flat = dak.flatten(data)
    else:
        data_flat = data
        
    histogram = dhist.Hist(h.axis.Regular(nBins, lo, hi))
    histogram.fill(data_flat)
    
    # Compute
    computed_hist = histogram.compute()
    computed_data = data_flat.compute()
    
    n_entries = len(computed_data)
    mean_val = np.mean(computed_data)
    std_val = np.std(computed_data)
    
    fig, ax = plt.subplots(1 ,1, figsize=(10, 8))
    mh.histplot(computed_hist, ax=ax, label=label, yerr=False, histtype="step", linewidth=2, color=color) 
    mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2016, loc=0)
    
    stat_text = (f"Entries: {n_entries}\n"
                 f"Mean:    {mean_val:.2f}\n"
                 f"Std Dev: {std_val:.2f}")
    at = AnchoredText(stat_text, prop=dict(size=10), frameon=True, loc='upper right')
    # at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)
    
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    ax.legend(loc='upper left')
    
    if logy: ax.set_yscale("log")
    
    save_path = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(save_path)
    print(f"Saved: {save_path}")
    plt.close()

# =================================================================
# 2. LOAD DATA
# =================================================================
print("Loading Events...")
events = NanoEventsFactory.from_root(
    {"./datasets/nano106X_on_mini106X_2017_mc_NANOAOD_W1Jets_to_LNu_250K.root": "Events"},
    schemaclass=NanoAODSchema,
    metadata={"dataset":"WJets"},
    mode="dask"
).events()

# =================================================================
# 3. PRE-SELECTION & CONTROL PLOTS
# =================================================================
print("Running Pre-selection...")

nMu_raw = dak.num(events.Muon.pt)
make_hist(nMu_raw, 10, 0, 10, "Number of Muons", "Events", "Raw Events", "01_nMu_Raw.png", logy=True)

events_2mu = events[nMu_raw >= 2]
Muon_2mu = events_2mu.Muon
MET = events_2mu.MET

make_hist(Muon_2mu.pt[:, 0], 50, 0, 200, r"Leading Muon $p_T$ [GeV]", "Events", r"$>= 2 Mu$ Events", "02_LeadingMuPt.png")

# =================================================================
# 4. OBJECT DEFINITION
# =================================================================
print("Defining Physics Objects...")

Mt = np.sqrt(2 * Muon_2mu.pt * MET.pt * (1 - np.cos(Muon_2mu.phi - MET.phi)))
W_Mask = (
    (Muon_2mu.pt > 30) & 
    (abs(Muon_2mu.eta) < 2.4) &
    (Muon_2mu.tightId == True) &
    (Muon_2mu.pfRelIso04_all < 0.15) &
    (Mt > 50)
)
W_Muons = Muon_2mu[W_Mask]

Soft_Mask = (
    (Muon_2mu.pt < 20) &
    (abs(Muon_2mu.eta) < 2.4) &
    (Muon_2mu.pfRelIso04_all > 0.2) &
    (Muon_2mu.tightId == True) &
    (Muon_2mu.jetIdx >= 0) 
)
Soft_Muons = Muon_2mu[Soft_Mask]
Matched_Jets = events_2mu.Jet[Soft_Muons.jetIdx]

Jet_PassCuts = (Matched_Jets.pt > 30) & (abs(Matched_Jets.eta) < 2.4)
Good_C_Jets = Matched_Jets[Jet_PassCuts]
Good_Soft_Muons = Soft_Muons[Jet_PassCuts]

# =================================================================
# 5. FINAL SELECTION & TAGGING PLOTS
# =================================================================
nW = dak.num(W_Muons.pt)
nC = dak.num(Good_C_Jets.pt)

make_hist(nW, 5, 0, 5, "W Candidates / Event", "Events", "W Tag", "03_nW_Candidates.png", logy=True)
make_hist(nC, 5, 0, 5, "C Candidates / Event", "Events", "C Tag", "04_nC_Candidates.png", logy=True)

final_mask = (nW == 1) & (nC >= 1)

Final_W = W_Muons[final_mask]
Final_C_Jets = Good_C_Jets[final_mask]
Final_Soft_Muons = Good_Soft_Muons[final_mask]

make_hist(dak.num(Final_C_Jets.pt), 5, 0, 5, "C-Jets (Final)", "Events", "Final Sel", "05_nC_Final.png", logy=True)

# =================================================================
# 6. SS / OS ANALYSIS & SUBTRACTION
# =================================================================
print("Running SS/OS Analysis...")

# Charge Comparison
w_charge = Final_W.charge[:, 0] 
soft_charge = Final_Soft_Muons.charge

is_OS = (w_charge != soft_charge)
is_SS = (w_charge == soft_charge)

OS_Jets = Final_C_Jets[is_OS]
SS_Jets = Final_C_Jets[is_SS]

# Create Histograms 
axis = h.axis.Regular(20, 30, 150, name="pt", label=r"C-Jet $p_T$ [GeV]")

# Use storage=h.storage.Weight() to ensure variances are calculated
h_os = dhist.Hist(axis, storage=h.storage.Weight()).fill(dak.flatten(OS_Jets.pt)).compute()
h_ss = dhist.Hist(axis, storage=h.storage.Weight()).fill(dak.flatten(SS_Jets.pt)).compute()

# Extract numbers safely (handling the .value attribute if it exists)
n_os = h_os.sum().value if hasattr(h_os.sum(), "value") else h_os.sum()
n_ss = h_ss.sum().value if hasattr(h_ss.sum(), "value") else h_ss.sum()

# --- Plot Overlay (OS vs SS) ---
fig, ax = plt.subplots(figsize=(10, 8))
mh.histplot(h_os, ax=ax, label=f"OS (N={n_os:.0f})", histtype="step", linewidth=2, color="tab:blue")
mh.histplot(h_ss, ax=ax, label=f"SS (N={n_ss:.0f})", histtype="step", linewidth=2, color="tab:red")
mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2017, loc=0)
ax.set_ylabel("Jets")
ax.legend()
plt.savefig(os.path.join(OUTPUT_DIR, "06_OS_vs_SS.png"))
print(f"Saved: {os.path.join(OUTPUT_DIR, '06_OS_vs_SS.png')}")
plt.close()

# --- Plot Subtracted (OS - SS) ---
# FIX: Instead of subtracting objects (h_os - h_ss) which caused the error,
# we subtract the NumPy arrays directly.

# 1. Extract Values and Variances
val_os = h_os.values()
val_ss = h_ss.values()

var_os = h_os.variances()
var_ss = h_ss.variances()

# 2. Compute Difference and Error
# Diff = OS - SS
diff_val = val_os - val_ss

# Error = Sqrt(Var_OS + Var_SS) (Standard error propagation for subtraction)
# If variances are somehow None (unweighted), fallback to sqrt(N)
if var_os is None: var_os = val_os
if var_ss is None: var_ss = val_ss
diff_err = np.sqrt(var_os + var_ss)

n_diff = np.sum(diff_val)

# 3. Plot using mplhep.histplot with arrays
fig, ax = plt.subplots(figsize=(10, 8))

# Pass bins explicitly since we are passing arrays
mh.histplot(
    diff_val, 
    bins=axis.edges, 
    yerr=diff_err, 
    ax=ax, 
    label=f"OS - SS (N={n_diff:.0f})", 
    histtype="errorbar", 
    color="black", 
    marker='o'
)

# Add a zero line
ax.axhline(0, color='gray', linestyle='--', linewidth=1)

mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2017, loc=0)
ax.set_ylabel(r"Events (OS - SS)")
ax.set_xlabel(r"C-Jet $p_T$ [GeV]")
ax.legend()

plt.savefig(os.path.join(OUTPUT_DIR, "07_OS_minus_SS.png"))
print(f"Saved: {os.path.join(OUTPUT_DIR, '07_OS_minus_SS.png')}")
plt.close()

print("Analysis Complete.")