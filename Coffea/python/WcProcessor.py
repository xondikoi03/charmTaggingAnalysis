# --------------------------------- #
# Coffea Processor for W+c Analysis 
# Original Author: Raj Handique
# Date Created: 21.11.2025
# --------------------------------- #

# Importing required libraries
import dask
import hist as h
from coffea import processor
from coffea.nanoevents import NanoAODSchema, NanoEventsFactory
import dask_awkward as dak
import hist.dask as dhist
import awkward as ak

# WcProcessor:
class WCProcessor(processor.ProcessorABC):
    def __init__(self):
        pass
    
    def make_hist(self, name, nBins, lo, hi, label):
        return(
            h.Hist.new
            .Reg(nBins, lo, hi, name=name, label=label)
        )

    def WBoson_process(self, events):
        dataset = events.metadata["dataset"]
        RecoMuons = events.Muon

        # Muon Selection Cuts: 
        Selection_Cuts = (
            (RecoMuons.pt > 35) &
            (abs(RecoMuons.eta) < 2.4) &
            (RecoMuons.pfRelIso04_all) < 0.15 &
            (RecoMuons.tightID == 1)
        )

        # M_T Selection Cut:
        MET = events.MET
        Muon_Mt = np.sqrt(2 * RecoMuons.pt * MET.pt * (1 - np.cos(RecoMuons.phi - MET.phi)))
        Selection_Cuts = Selection_Cuts & (Muon_Mt > 55)
        Muons_Sel = RecoMuons[Selection_Cuts]

        # Leading Muon Pt:
        Sorted_Muons = Muons_Sel[ak.argsort(Muons_Sel.pt, ascending=False)]
        lead_Muon_Pt = ak.firsts(Sorted_Muons.pt)
        # Subleading Muon Pt:
        sublead_Muon_Pt = ak.firsts(Sorted_Muons.pt[:, 1])

        # Number of Sorted Muons & passCuts for OS:
        nMuons = ak.num(Sorted_Muons)
        OS_mask = (
            (ak.num(Sorted_Muons.charge, axis=1) == 0) |
            (nMuons >= 2) |
            (sublead_Muon_Pt > 20)
        ) 

        # Final Selection of Muons:
        fMuons = Sorted_Muons[~OS_mask]

    # def charm_process(self, events):

    def postprocess(self, accumulator):
        return accumulator

