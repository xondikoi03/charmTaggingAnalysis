import os
import numpy as np
import awkward as ak
import dask_awkward as dak
import dask
import hist
import hist.dask as dhist
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import matplotlib.pyplot as plt
import mplhep as mh
from matplotlib.offsetbox import AnchoredText

# =================================================================
# 1. SETUP & UTILITIES
# =================================================================
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Plots will be saved to: {os.path.abspath(OUTPUT_DIR)}/")

# CMS Style
plt.style.use(mh.style.CMS)
plt.rcParams.update({"font.size": 12})

def delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    return (dphi + np.pi) % (2 * np.pi) - np.pi

def save_plot(fig, filename):
    path = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(path, bbox_inches="tight")
    print(f"Saved: {path}")
    plt.close(fig)

# --- CUSTOM PLOTTING FUNCTION WITH STATS ---
def plot_with_stats(h_computed, data_computed, label, xlabel, filename, logy=False, color="black"):
    """
    Plots a histogram with a CMS-style statistics box (Entries, Mean, Std Dev).
    Uses pre-computed data to avoid triggering dask computes during plotting.
    """
    # Calculate Stats
    n_entries = len(data_computed)
    mean_val = np.mean(data_computed)
    std_val = np.std(data_computed)

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # Plot Histogram
    mh.histplot(h_computed, ax=ax, label=label, yerr=False, histtype="step", linewidth=2, color=color)
    mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2016, loc=0)

    # Create Stats Box
    stat_text = (f"Entries: {n_entries}\n"
                 f"Mean:    {mean_val:.2f}\n"
                 f"Std Dev: {std_val:.2f}")
    
    at = AnchoredText(stat_text, prop=dict(size=10), frameon=True, loc='upper right')
    ax.add_artist(at)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Events")
    ax.legend(loc='upper left')

    if logy: ax.set_yscale("log")
    
    save_plot(fig, filename)

# =================================================================
# 2. LOAD DATA
# =================================================================
fname = "./datasets/nano106X_on_mini106X_2017_mc_NANOAOD_W1Jets_to_LNu_250K.root"
print(f"Loading: {fname}")

events = NanoEventsFactory.from_root(
    {fname: "Events"},
    schemaclass=NanoAODSchema,
    metadata={"dataset": "WJets"},
    mode="dask"
).events()

# Collections to store things we want to compute
hists_to_compute = {}   # The histograms
stats_to_compute = {}   # The raw arrays (for Mean/Std calc)
counters = {}           # The Cutflow numbers

# --- CUTFLOW 1: INITIAL ---
counters["1. Sample size"] = dak.num(events, axis=0)

# =================================================================
# 3. CONTROL PLOTS DEFINITION
# =================================================================
print("Defining Physics Objects...")

muons = events.Muon
jets = events.Jet
met = events.MET

# Transverse Mass
dphi_w = delta_phi(muons.phi, met.phi)
mt_w = np.sqrt(2 * muons.pt * met.pt * (1 - np.cos(dphi_w)))

nMu_raw = dak.num(muons, axis=1)
has_muon = nMu_raw > 0

# -- 1. Raw Jet Count --
hists_to_compute["nJet"] = dhist.Hist(hist.axis.Regular(16, 0, 16, name="nJet")).fill(dak.num(jets, axis=1))
stats_to_compute["nJet"] = dak.num(jets, axis=1)

# -- 2. Jet pT --
hists_to_compute["Jet_pt"] = dhist.Hist(hist.axis.Regular(50, 0, 150, name="Cjet_pt")).fill(dak.flatten(jets.pt))
stats_to_compute["Jet_pt"] = dak.flatten(jets.pt)

# -- 3. Raw Muon Count --
hists_to_compute["0nMu"] = dhist.Hist(hist.axis.Regular(8, 0, 8, name="0nmu")).fill(nMu_raw)
stats_to_compute["0nMu"] = nMu_raw # Keep raw data for stats

# -- 4. Muon Count --
hists_to_compute["nMu"] = dhist.Hist(hist.axis.Regular(8, 0, 8, name="nmu")).fill(nMu_raw[has_muon])
stats_to_compute["nMu"] = nMu_raw[has_muon] # Keep raw data for stats

# -- 5. Muon pT --
hists_to_compute["Mu_pt"] = dhist.Hist(hist.axis.Regular(50, 0, 150, name="Wmu_pt")).fill(dak.flatten(muons.pt[has_muon]))
stats_to_compute["Mu_pt"] = dak.flatten(muons.pt[has_muon])

# -- 6. Transverse Mass --
hists_to_compute["mt"] = dhist.Hist(hist.axis.Regular(50, 0, 150, name="mt")).fill(dak.flatten(mt_w[has_muon]))
stats_to_compute["mt"] = dak.flatten(mt_w[has_muon])



# =================================================================
# 4. SELECTION
# =================================================================
print("Defining selection cuts")

# W Selection
w_muon_mask = (
    (muons.pt > 30) &             
    (abs(muons.eta) < 2.4) &
    (muons.tightId) &             
    (muons.pfRelIso04_all < 0.15) & 
    (mt_w > 50)                   
)
W_Muons = muons[w_muon_mask]

# --- CUTFLOW 2: W BOSON ---
has_w = dak.num(W_Muons, axis=1) >= 1
counters["2. W_Mu cuts"] = dak.sum(has_w)

# Charm Selection
soft_muon_mask = (
    (muons.pt < 20) &             
    #(muons.pt > 3) &              
    (abs(muons.eta) < 2.4) &
    # (muons.looseId) &
    (muons.tightId) &
    (muons.pfRelIso04_all > 0.2) &
    (muons.jetIdx >= 0)           
)

Soft_Muons = muons[soft_muon_mask]
Matched_Jets = events.Jet[Soft_Muons.jetIdx]

c_jet_mask = (Matched_Jets.pt > 30) & (abs(Matched_Jets.eta) < 2.4)
Good_C_Jets = Matched_Jets[c_jet_mask]
Good_Soft_Muons = Soft_Muons[c_jet_mask]

# --- CUTFLOW 3: C TAG ---
has_c = dak.num(Good_C_Jets, axis=1) >= 1
counters["3. Soft C_Mu cuts"] = dak.sum(has_c)

# =================================================================
# 5. ANALYSIS REGION
# =================================================================
nW = dak.num(W_Muons, axis=1)
has_w = nW >= 1
# -- 7. n Mu after W-cut --
hists_to_compute["W-cut_nMu"] = dhist.Hist(hist.axis.Regular(8, 0, 8, name="nW")).fill(dak.num(W_Muons, axis=1)[has_w])
stats_to_compute["W-cut_nMu"] = dak.num(W_Muons, axis=1)[has_w] # Keep raw data for stats


nC = dak.num(Good_C_Jets, axis=1)
has_c = nC >= 1
# -- 8. n C after C-cut --
hists_to_compute["C-cut_nC"] = dhist.Hist(hist.axis.Regular(8, 0, 8, name="nC")).fill(nC[has_c])
stats_to_compute["C-cut_nC"] = nC[has_c] # Keep raw data for stats

# -- 9. Soft Muon Pt (Probe) --
hists_to_compute["softMu_pt"] = dhist.Hist(hist.axis.Regular(20, 0, 25, name="pt")).fill(dak.flatten(Good_Soft_Muons.pt[has_c]))
stats_to_compute["softMu_pt"] = dak.flatten(Good_Soft_Muons.pt[has_c])

# -- 10. JEt Muon Pt (Probe) --
hists_to_compute["c_pt"] = dhist.Hist(hist.axis.Regular(50, 30, 200, name="pt")).fill(dak.flatten(Good_C_Jets.pt[has_c]))
stats_to_compute["c_pt"] = dak.flatten(Good_C_Jets.pt[has_c])

event_mask = (nW == 1) & (nC == 1)

# --- CUTFLOW 4: SIGNAL REGION ---
counters["4. Signal Region"] = dak.sum(event_mask)

Final_W = W_Muons[event_mask]
Final_C = Good_C_Jets[event_mask]
Final_Soft = Good_Soft_Muons[event_mask]

# OS / SS Logic
w_charge = Final_W.charge[:, 0]
soft_charge = Final_Soft.charge
w_charge_broadcast = dak.broadcast_arrays(w_charge, soft_charge)[0]

is_OS = w_charge_broadcast != soft_charge
is_SS = w_charge_broadcast == soft_charge

# Final Physics Histograms
axis_final = hist.axis.Regular(15, 30, 150, name="pt", label=r"Charm Jet $p_T$ [GeV]")
h_os = dhist.Hist(axis_final, storage=hist.storage.Weight()).fill(dak.flatten(Final_C[is_OS].pt))
h_ss = dhist.Hist(axis_final, storage=hist.storage.Weight()).fill(dak.flatten(Final_C[is_SS].pt))

# =================================================================
# 6. MASSIVE COMPUTE
# =================================================================
print("Running Compute (Physics + Stats + Cutflow)...")

# We zip everything into one list to compute efficiently
comp_keys_h = list(hists_to_compute.keys())
comp_vals_h = list(hists_to_compute.values())

comp_keys_s = list(stats_to_compute.keys())
comp_vals_s = list(stats_to_compute.values())

comp_keys_c = list(counters.keys())
comp_vals_c = list(counters.values())

# The Big Compute
results = dask.compute(
    h_os, h_ss,         # 0, 1
    *comp_vals_h,       # Histograms
    *comp_vals_s,       # Raw Arrays (for stats)
    *comp_vals_c        # Counters
)

# Unpack
res_os, res_ss = results[0], results[1]

# Unpack Lists using slicing
offset = 2
res_hists = dict(zip(comp_keys_h, results[offset : offset + len(comp_keys_h)]))
offset += len(comp_keys_h)
res_stats = dict(zip(comp_keys_s, results[offset : offset + len(comp_keys_s)]))
offset += len(comp_keys_s)
res_counts = dict(zip(comp_keys_c, results[offset : offset + len(comp_keys_c)]))

print("Computation Complete.")

# =================================================================
# 7. PLOTTING
# =================================================================

# --- A. Control Plots with Stats Box ---
print("Generating Control Plots...")

plot_with_stats(
    res_hists["nJet"], res_stats["nJet"],
    "Raw Events", r"Number of Reco Jets", "01_Raw_nJet.png", logy=True
)
plot_with_stats(
    res_hists["Jet_pt"], res_stats["Jet_pt"],
    "Raw Events", r"Jet $p_T$ [GeV]", "02_Raw_JetPt.png", logy=True
)

plot_with_stats(
    res_hists["0nMu"], res_stats["0nMu"], 
    "Raw Events", "Number of Reco Muons", "03_Raw_nMu.png", logy=True
)
plot_with_stats(
    res_hists["nMu"], res_stats["nMu"], 
    "W Candidates", "Number of Reco Muons", "04_Raw_nMu.png", logy=True
)
plot_with_stats(
    res_hists["Mu_pt"], res_stats["Mu_pt"],
    "W Candidates", r"W Jet $p_T$ [GeV]", "05_Raw_MuPt.png"
)
plot_with_stats(
    res_hists["mt"], res_stats["mt"],
    "W Candidates", r"Transverse Mass $M_T$ [GeV]", "06_Raw_mT.png"
)

plot_with_stats(
    res_hists["W-cut_nMu"], res_stats["W-cut_nMu"],
    "Muons passing W Cuts", r"Number of Muons", "07_W-cut_nMu.png"
)

plot_with_stats(
    res_hists["C-cut_nC"], res_stats["C-cut_nC"],
    "Muons passing C Cuts", r"Number of C Candidates", "08_C-cut_nC.png"
)
plot_with_stats(
    res_hists["softMu_pt"], res_stats["softMu_pt"],
    "Muons passing C Cuts", r"Muons passing C Cuts $p_T$ [GeV]", "09_C-cut_SoftMu_pt.png"
)
plot_with_stats(
    res_hists["c_pt"], res_stats["c_pt"],
    "C Jets", r"C Jet $p_T$ [GeV]", "10_C-cut_CJet_pt.png"
)

# --- B. Smart Cutflow ---
def plot_cutflow(counts_dict):
    fig, ax = plt.subplots(figsize=(10, 8))
    labels = list(counts_dict.keys())
    values = list(counts_dict.values())
    
    bars = ax.bar(labels, values, color='teal', alpha=0.7, edgecolor='black')
    initial = values[0]
    
    for i, (bar, val) in enumerate(zip(bars, values)):
        ax.text(bar.get_x() + bar.get_width()/2, val, f"{int(val):,}", 
                ha='center', va='bottom', fontweight='bold')
        if i > 0:
            eff_step = (val / values[i-1]) * 100
            eff_tot = (val / initial) * 100
            text_y = val / 2 if val > 0 else 0
            ax.text(bar.get_x() + bar.get_width()/2, text_y, 
                   f"Step: {eff_step:.1f}%\nTot: {eff_tot:.2f}%", 
                   ha='center', va='center', color='white', fontweight='bold', fontsize=10)

    ax.set_yscale('log')
    ax.set_ylim(bottom=1)
    ax.set_ylabel("Events")
    ax.set_title("Event Selection Cutflow")
    mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2016, loc=0)
    save_plot(fig, "00_Analysis_Cutflow.png")
    
    print("\n=== CUTFLOW SUMMARY ===")
    for k, v in counts_dict.items():
        print(f"{k:20s}: {int(v):,}")

plot_cutflow(res_counts)

# --- C. Signal Extraction (OS-SS) ---
# Overlay
fig, ax = plt.subplots(figsize=(10, 8))
n_os = res_os.sum(flow=True).value 
n_ss = res_ss.sum(flow=True).value

mh.histplot(res_os, ax=ax, label=f"OS (N={n_os:.0f})(Signal+Bkg)", color="tab:blue")
mh.histplot(res_ss, ax=ax, label=f"SS (N={n_ss:.0f})(Bkg Only)", color="tab:red")
mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2016, loc=0)
ax.legend()
save_plot(fig, "11_Analysis_OS_vs_SS.png")

# Subtraction
y_diff = res_os.values() - res_ss.values()
err_diff = np.sqrt(res_os.variances() + res_ss.variances())
n_sig = np.sum(y_diff)
n_err = np.sqrt(np.sum(res_os.variances() + res_ss.variances()))

fig, ax = plt.subplots(figsize=(10, 8))
mh.histplot(
    y_diff, bins=axis_final.edges, yerr=err_diff,
    ax=ax, label=r"Signal ($W+c$)", color="black", 
    histtype="errorbar", marker='o'
)
ax.axhline(0, color="gray", linestyle="--")
mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2016, loc=0)
ax.set_xlabel(r"Charm Jet $p_T$ [GeV]")
ax.set_ylabel(r"Signal Events (OS - SS)")

# Result Box
at = AnchoredText(f"Signal Yield:\n{n_sig:.1f} $\\pm$ {n_err:.1f}", loc="upper right", frameon=False)
ax.add_artist(at)
ax.legend(loc="upper left")

save_plot(fig, "12_Analysis_Signal_Subtracted.png")

print("Analysis Complete.")