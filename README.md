# ExamplesFCCFullSim

Starting level package for tau analysis with full sim.

The setup runs over Pythia8 samples stored in the CIEMAT pnfs system. 
- 1M events of ZTauTau Pythia8 decays.
- FullSim: CLD, modified muon reconstruction at low momentum (MinTrackCandidateEnergy=4)
- Key4Hep version: 2024-10-03
- Sample location:  /pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/ZTauTau_SMPol_25Sept_MuonFix

The details of the simple algorithm can be found in modules/tauReco

How to login to the cluster: 
```
ssh user@pcae22.ciemat.es -p 622
ssh gaeui0X (X=2,3,4,5)
``` 

Fast start for plots:
```
git clone https://github.com/mcepeda/ExamplesFCCFullSim.git
cd ExamplesFCCFullSim

# Setup of key4hep (version compatible with the fullsim files generated in Autumn)
source setupKey4Hep.sh

# Run some simple tau checks and fill histograms  
python plotTausStart.py

# Make a plot 
python simplePlot.py
```

