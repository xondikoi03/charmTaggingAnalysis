# --------------------------------- #
# charmTagging - W + c-jet analysis #
# Original author: Raj Handique     #
# Date: November 2025               #
# --------------------------------- #

# Neccessary Libraries and Packages
import ROOT as rt
import uproot as uprt
import mplhep as mh
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np

# Defining the style for plots (CMS Style)
plt.style.use(mh.style.CMS)
plt.rcParams.update({
    "font.size": 17,
    "axes.labelsize": 17,
    "axes.titlesize": 17,
    "xtick.labelsize": 17,
    "ytick.labelsize": 17,
    "legend.fontsize": 17
})

# Defining the function for making histograms
# Histogram function
# Functions for plotting
def make_hist(data, bins, range, xlabel, ylabel, label, fname=None, logy=False):
    histo, edges = np.histogram(data, bins=bins, range=range)
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    mh.histplot(histo, edges, histtype="step", label=label, ax=ax)
    mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2016, loc=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()
    if logy:
        ax.set_yscale("log")
    if fname:
        plt.savefig(fname)
    plt.show()

# Stacked histogram function
def make_stack_hist(data, labels, bins, range, xlabel, ylabel, fname=None, logy=False):
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    plt.hist(
        data,
        bins=bins,
        range=range,
        histtype="stepfilled",
        label=labels,
        density=True,
        stacked=True
    )

    mh.cms.label("Open Data", data=True, lumi=None, com=13, year=2016, loc=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()
    if logy:
        ax.set_yscale("log")
    if fname:
        plt.savefig(fname)
    plt.show()

# Loading the ROOT File 
btvNanoAODFile_1K = uprt.open("~/Desktop/charmTaggingAnalysis/datasets/nano106X_on_mini106X_2017_mc_NANOAOD_W1Jets_to_LNu.root")
Events = btvNanoAODFile_1K["Events"]

# RecoJets collections
RecoJet_keys = Events.keys(filter_name="Jet_*")
for RecoJets in Events.iterate(RecoJet_keys, step_size=Events.num_entries, library="ak"):
    break

# Muons and electrons Collections
RecoMuon_keys = Events.keys(filter_name="Muon_*")
for RecoMuons in Events.iterate(RecoMuon_keys, step_size=Events.num_entries, library="ak"):
    break

RecoElectron_keys = Events.keys(filter_name="Electron_*")
for RecoElectrons in Events.iterate(RecoElectron_keys, step_size=Events.num_entries, library="ak"):
    break

# Photon Collections
RecoPhoton_keys = Events.keys(filter_name="Photon_*")
for RecoPhotons in Events.iterate(RecoPhoton_keys, step_size=Events.num_entries, library="ak"):
    break

# MET collection
MET_keys = Events.keys(filter_name="MET_*")
for METs in Events.iterate(MET_keys, step_size=Events.num_entries, library="ak"):
    break

# ---------------------------- #
# W+C Jet Analysis Starts Here #
# ---------------------------- #

# Event Selection : W+c 
# Muons: HLT IsoMu24 OR HLT IsoTkMu24
# Rochester Corrections for Muons
# Muon Pt > 35 GeV, |eta| < 2.4
# I_comb = Sum(Pt(chargedhadrons)) + Et(neutralhadrons) + Et(photons) - Sum(pt(pileup))
# I_comb / pt(muon) < 0.15
# M_{T} > 55 GeV

# Jet Selection Cuts
passCuts = (RecoJets["Jet_pt"] > 30) & (abs(RecoJets["Jet_eta"]) < 2.4)
RecoJet_EvtSel = RecoJets[passCuts]

# This is the isolation cut (I_comb)
Muon_Pt_Cut = (RecoMuons["Muon_pt"] > 35)
Muon_pfIso_Cut = (abs(RecoMuons["Muon_pfRelIso04_all"] < 0.15))
Muon_Eta_Cut = (abs(RecoMuons["Muon_eta"]) < 2.4)

# Missing Trasnverse Mass Cut
# M_t = sqrt(2 * pt_muon * MET_pt * (1 - cos(mu_phi - MET_phi)))
Muon_Mt = np.sqrt(2 * RecoMuons["Muon_pt"] * METs["MET_pt"] * (1 - np.cos(RecoMuons["Muon_phi"] - METs["MET_phi"])))
Muon_Mt_Cut = (Muon_Mt > 55)  # Apply Mt cut to all tight muons

# Tight Muons 
tightMuons_Selection = Muon_Pt_Cut & Muon_pfIso_Cut & Muon_Eta_Cut & Muon_Mt_Cut  
tightMuons = RecoMuons[tightMuons_Selection]

# Pt of Leading Muon (Tight Selection)                             
Sort_tightMuons = tightMuons[ak.argsort(tightMuons["Muon_pt"], ascending=False)]
Leading_Muon_pt = ak.firsts(Sort_tightMuons["Muon_pt"])
# Leading_Muon_Mt = ak.firsts(Muon_Mt)
# Muon_Mt_Cut = (Leading_Muon_Mt > 55) # Took the Mt of the leading muon only

# Number of tight muons in the event
nTightMuons = (ak.num(Sort_tightMuons["Muon_pt"]))

# Pt of Subleading Muon (Tight Selection) 
padded_muons = ak.pad_none(Sort_tightMuons, 2, axis=1)
Subleading_Muon_pt = ak.fill_none(padded_muons[:, 1].Muon_pt, -1) #Accesses the second element of the padded array and fills None with -1

Subleading_Muon_pt_Cut = (Subleading_Muon_pt > 20) # Just a placeholder cut

# Subleading_Muon_pt = ak.pad_none(Sort_tightMuons["Muon_pt"], 2, axis=1)[:, 1] # This probably won't work
# Subleading_Muon_pt = ak.fill_none(Subleading_Muon_pt, -1)

# Selection cuts to reduce Drell-Yan contributions (OS muons)
nTightMuons_mask = (nTightMuons == 2) # Events that have two tight muons

# Sum of charges of the two tight muons should be zero (opposite sign)
tightMuon_OS_Cut = (ak.sum(Sort_tightMuons["Muon_charge"], axis=1) == 0)
OS_mask = nTightMuons_mask & tightMuon_OS_Cut & Subleading_Muon_pt_Cut
Muons_to_discard = ~OS_mask

# Final Selection Cuts
Final_Muons_Cuts = Muons_to_discard # Muon_Mt_Cut doesn't seem to work for 10K events
Selected_Muons = Sort_tightMuons[Final_Muons_Cuts]

# Some plots to visualize the selected muons
Muon_pt = ak.to_numpy(ak.flatten(Selected_Muons["Muon_pt"]))

make_hist(
    nTightMuons,
    bins=50,
    range=(0, 5),
    xlabel="Number of Muons",
    ylabel="Events",
    label="Selected Tight Muons",
    fname="Plots/nMuon_WcAnalysis.png",
    logy=False
)

make_hist(
    Leading_Muon_pt,
    bins=50,
    range=(0, 200),
    xlabel="Leading Muon $p_{T}$ [GeV]",
    ylabel="Events",
    label="Selected Tight Muons",
    fname="Plots/Leading_Muon_pt_WcAnalysis.png",
    logy=False
)

make_hist(
    Subleading_Muon_pt,
    bins=50,
    range=(0, 200),
    xlabel="Subleading Muon $p_{T}$ [GeV]",
    ylabel="Events",
    label="Selected Tight Muons",
    fname="Plots/Subleading_Muon_pt_WcAnalysis.png",
    logy=False
)

make_hist(
    Muon_pt,
    bins=50,
    range=(0, 200),
    xlabel="Muon $p_{T}$ [GeV]",
    ylabel="Events",
    label="Selected Tight Muons",
    fname="Plots/Muon_pt_WcAnalysis.png",
    logy=False
)

# Number of selected muons after Wc event selection
nSelected_Muons = ak.num(Selected_Muons["Muon_pt"])

make_hist(
    nSelected_Muons,
    bins=50,
    range=(0, 5),
    xlabel="Number of Selected Muons",
    ylabel="Events",
    label="Muons after Wc Selection",
    fname="Plots/nSelected_Muons_WcAnalysis.png",
    logy=False
)