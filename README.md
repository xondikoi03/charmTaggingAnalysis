# Charm Tagging Analysis:

## Setting up the Environment:

This specific analysis is done using `Coffea:2024.10.0` in `conda`. Here is how I set up the environment:

1. Set up the `conda` environment:

```bash
conda create -n Coffea python=3.12.12
conda activate Coffea
```

2. To install `Coffea:2024.10.0`:

```
conda install -c conda-forge coffea=2024.10.0
```

3. Is you want to install using `pip`:
```
pip install coffea==2024.10.0
```

4. Additionally to work on a jupyter notebook:
```
conda install jupyter
```

This repository contains resources used in the analysis of charm tagging of jets. This README will be updated according to the status of the Analysis.

## Status of Analysis:

- **Event Selection:** W+c Events
- **Selection Cuts**: 
    - For Jets:
        1. Jet $p_{T}$ > 30 GeV
        2. Jet $\eta$ < 2.4.
    - For W-> $\mu + \nu$ Events:
        1. Muon $p_{T}$ > 20 GeV.
        2. Muon_pfRelIso04_all < 0.15
        3. $M_{T}$ > 55 GeV
        4. Events with two tight Muons with OS $p_{T}$ > 20 GeV are to be removed.
    - For W-> $e + \nu$ Events:
        1. Electron $p_{T}$ > 20 GeV.
        2. Electron_pfRelIso03_all < 0.15
        3. $M_{T}$ > 55 GeV
        4. Events with two tight Electrons with OS $p_{T}$ > 20 GeV are to be removed.