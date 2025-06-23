# Mini Tree with diphoton information

To load root better use the same version I used to create the tree.

```source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-10-03```

The trees have been created with the script miniTreeForAnneMarie.py. 

The source (bigger) edm4hep roottrees are at the CIEMAT cluster, I can move part of it to CERN

A sample tree to start with with only 10000 events is in 

```root -b miniTree.root```

Contents of the file:

a) some example histograms with the energy of muons, electrons, photons, genPhotons,and GenPi0s 

b) a tree with branches with 
- number of photons, genphotons, genPi0s, and taus in the event
- 4 vector components for the photons, genphotons and Pi0s, as vectors
  
```
root [3] outtree->Show(1)
======> EVENT:1
 beamE           = 45.5991
 nPhotons        = 0
 nGenPhotons     = 4
 nGenPi0s        = 0
 nGenTaus        = 2
 nRecoTausHad    = 1
 photonP         = (vector<double>*)0x2ea9fb0
 photonE         = (vector<double>*)0x2ea9f60
 photonPx        = (vector<double>*)0x31cc7a0
 photonPy        = (vector<double>*)0x2e45e20
 photonPz        = (vector<double>*)0x31ee420
 photonM         = (vector<double>*)0x2fb6710
 genPhotonP      = (vector<double>*)0x32100a0
 genPhotonE      = (vector<double>*)0x27ded90
 genPhotonPx     = (vector<double>*)0x27e02e0
 genPhotonPy     = (vector<double>*)0x3243ee0
 genPhotonPz     = (vector<double>*)0x3253920
 genPhotonM      = (vector<double>*)0x3263360
 genPi0P         = (vector<double>*)0x3272da0
 genPi0E         = (vector<double>*)0x32827e0
 genPi0Px        = (vector<double>*)0x3292220
 genPi0Py        = (vector<double>*)0x32a1c60
 genPi0Pz        = (vector<double>*)0x32b16a0
 genPi0M         = (vector<double>*)0x32c10e0
```

