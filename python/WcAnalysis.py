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

# Loading the dataset
NanoAODFile = uprt.open("datasets/W1JetsToLNu_TuneCP5_13TeV.root")
Events = NanoAODFile["Events"]

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

# ---------------------------- #
# W+C Jet Analysis Starts Here #
# ---------------------------- #

passCuts = (RecoJets["Jet_pt"] > 30) & (abs(RecoJets["Jet_eta"]) < 2.4)
RecoJet_Selection = RecoJets[passCuts]

muonPassCuts = (RecoMuons["Muon_pt"] > 35) & (abs(RecoMuons["Muon_eta"]) < 2.4)
MuonsForAnalysis = RecoMuons[muonPassCuts]
# I_comb = (RecoMuons["Muon_pfIsoChargedHadronPt"] + RecoMuons["Muon_pfIsoNeutralHadronEt"] + RecoPhotons["Muon_pfIsoPhotonEt"])


# Muons: HLT IsoMu24 OR HLT IsoTkMu24
# Rochester Corrections for Muons
# Muon Pt > 35 GeV, |eta| < 2.4
# I_comb = pt(chargedhadrons) + Et(neutralhadrons) + Et(photons) - Pt(pileup) 
# I_comb / pt(muon) < 0.15
# MET pt > 55 GeV

