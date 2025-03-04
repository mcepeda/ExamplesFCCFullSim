# ExamplesFCCFullSim

Starting level package for tau analysis with full sim.

The setup runs over Pythia8 samples stored in the CIEMAT pnfs system. 
- 1M events of ZTauTau Pythia8 decays.
- FullSim: CLD, modified muon reconstruction at low momentum (MinTrackCandidateEnergy=4)
- Key4Hep version: 2024-10-03
- Sample location:  /pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/ZTauTau_SMPol_25Sept_MuonFix

The details of the simple algorithm can be found in modules/tauReco

## First steps 

How to login to the CIEMAT cluster: 
```
ssh user@pcae22.ciemat.es -p 622 # needs a ciemat account
ssh gaeui0X  # X=2,3,4,5
``` 

Fast start: plot simple tau variables 

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

## Extra information to understand the edm4hep files

[Namespace in edm4hep ](https://edm4hep.web.cern.ch/namespaceedm4hep.html)

[key4hep main page: https://key4hep.github.io/key4hep-doc/ ](https://key4hep.github.io/key4hep-doc/ )

<details>
<summary>Accessing collection info: what can we ask of a ReconstructedParticle</summary>

Open python and ask:

```
import edm4hep
reco = edm4hep.ReconstructedParticle()
dir(reco)
```

Answer:

['__add__', '__assign__', '__bool__', '__class__', '__delattr__', '__destruct__', '__dict__', '__dir__', '__dispatch__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__invert__', '__le__', '__lt__', '__module__', '__mul__', '__ne__', '__neg__', '__new__', '__pos__', '__python_owns__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__rmul__', '__rsub__', '__rtruediv__', '__setattr__', '__sizeof__', '__smartptr__', '__str__', '__sub__', '__subclasshook__', '__truediv__', '__weakref__', 'clone', 'clusters_begin', 'clusters_end', 'clusters_size', 'getCharge', 'getClusters', 'getCovMatrix', 'getEnergy', 'getGoodnessOfPID', 'getMass', 'getMomentum', 'getObjectID', 'getPDG', 'getParticles', 'getReferencePoint', 'getStartVertex', 'getTracks', 'getType', 'id', 'isAvailable', 'isCompound', 'makeEmpty', 'particles_begin', 'particles_end', 'particles_size', 'tracks_begin', 'tracks_end', 'tracks_size', 'unlink']

```
p = edm4hep.MCParticle()
dir(p)
```
Answer: 

['BITBackscatter', 'BITCreatedInSimulation', 'BITDecayedInCalorimeter', 'BITDecayedInTracker', 'BITLeftDetector', 'BITOverlay', 'BITStopped', 'BITVertexIsNotEndpointOfParent', '__add__', '__assign__', '__bool__', '__class__', '__delattr__', '__destruct__', '__dict__', '__dir__', '__dispatch__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__invert__', '__le__', '__lt__', '__module__', '__mul__', '__ne__', '__neg__', '__new__', '__pos__', '__python_owns__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__rmul__', '__rsub__', '__rtruediv__', '__setattr__', '__sizeof__', '__smartptr__', '__str__', '__sub__', '__subclasshook__', '__truediv__', '__weakref__', 'clone', 'daughters_begin', 'daughters_end', 'daughters_size', 'getCharge', 'getColorFlow', 'getDaughters', 'getEndpoint', 'getEnergy', 'getGeneratorStatus', 'getMass', 'getMomentum', 'getMomentumAtEndpoint', 'getObjectID', 'getPDG', 'getParents', 'getSimulatorStatus', 'getSpin', 'getTime', 'getVertex', 'hasLeftDetector', 'id', 'isAvailable', 'isBackscatter', 'isCreatedInSimulation', 'isDecayedInCalorimeter', 'isDecayedInTracker', 'isOverlay', 'isStopped', 'makeEmpty', 'parents_begin', 'parents_end', 'parents_size', 'unlink', 'vertexIsNotEndpointOfParent']

</details>

<details>
<summary>What collections are there in my edm4hep file</summary>

```
podio-dump out_reco_edm4hep_edm4hep.root

input file: /eos/experiment/fcc/ee/datasets/mlpf/condor/train/Z_tautau_v6/124/out_reco_edm4hep_edm4hep.root

datamodel model definitions stored in this file: edm4hep

Frame categories in this file:
Name                      Entries
----------------------  ---------
metadata                        1
events                        100
configuration_metadata          1
################################### events: 0 ####################################
Collections:
Name                                  ValueType                                    Size  ID
------------------------------------  -----------------------------------------  ------  --------
BuildUpVertices                       edm4hep::Vertex                                 0  fd03f5d0
BuildUpVertices_RP                    edm4hep::ReconstructedParticle                  0  310a0f04
BuildUpVertices_RP_particleIDs        edm4hep::ParticleID                             0  3e282687
BuildUpVertices_V0                    edm4hep::Vertex                                 0  9bc85810
BuildUpVertices_V0_RP                 edm4hep::ReconstructedParticle                  0  7db8c9b7
BuildUpVertices_V0_RP_particleIDs     edm4hep::ParticleID                             0  ae66aab3
CalohitMCTruthLink                    edm4hep::MCRecoCaloParticleAssociation        869  1132953d
ClusterMCTruthLink                    edm4hep::MCRecoClusterParticleAssociation       7  e2727d90
DebugHits                             edm4hep::TrackerHitPlane                        0  c498ee0d
ECALBarrel                            edm4hep::CalorimeterHit                         7  4b6bf95c
ECalBarrelCollection                  edm4hep::SimCalorimeterHit                     32  407eb17b
ECalBarrelCollectionContributions     edm4hep::CaloHitContribution                  243  047dde9c
ECALEndcap                            edm4hep::CalorimeterHit                       800  868b64f8
ECalEndcapCollection                  edm4hep::SimCalorimeterHit                    972  0c09a156
ECalEndcapCollectionContributions     edm4hep::CaloHitContribution                 8789  934221f8
EfficientMCParticles                  edm4hep::MCParticle                             2  099a8717
EventHeader                           edm4hep::EventHeader                            1  d793ab91
HCALBarrel                            edm4hep::CalorimeterHit                         0  076db0c4
HCalBarrelCollection                  edm4hep::SimCalorimeterHit                     65  6daf3b41
HCalBarrelCollectionContributions     edm4hep::CaloHitContribution                  119  4450a8f2
HCALEndcap                            edm4hep::CalorimeterHit                        51  4ca8cf89
HCalEndcapCollection                  edm4hep::SimCalorimeterHit                    123  50d90aed
HCalEndcapCollectionContributions     edm4hep::CaloHitContribution                  359  82a31463
HCALOther                             edm4hep::CalorimeterHit                         1  3ba1679d
HCalRingCollection                    edm4hep::SimCalorimeterHit                     17  ad7abfca
HCalRingCollectionContributions       edm4hep::CaloHitContribution                   32  3b573866
InefficientMCParticles                edm4hep::MCParticle                             0  55be621e
InnerTrackerBarrelCollection          edm4hep::SimTrackerHit                          5  f57448f3
InnerTrackerBarrelHitsRelations       edm4hep::MCRecoTrackerHitPlaneAssociation       5  029be193
InnerTrackerEndcapCollection          edm4hep::SimTrackerHit                          1  6d724693
InnerTrackerEndcapHitsRelations       edm4hep::MCRecoTrackerHitPlaneAssociation       1  743732ae
ITrackerEndcapHits                    edm4hep::TrackerHitPlane                        1  55c7e3f3
ITrackerHits                          edm4hep::TrackerHitPlane                        5  3263f08b
LooseSelectedPandoraPFOs              edm4hep::ReconstructedParticle                  5  bc5df93f
LumiCalCollection                     edm4hep::SimCalorimeterHit                      0  45759015
LumiCalCollectionContributions        edm4hep::CaloHitContribution                    0  dab2a7ce
LumiCalHits                           edm4hep::CalorimeterHit                         0  f7a65068
MCParticles                           edm4hep::MCParticle                            20  a1cba250
MCParticlesSkimmed                    edm4hep::MCParticle                            17  16b252ac
MCPhysicsParticles                    edm4hep::MCParticle                            20  8fe81b12
MCTruthClusterLink                    edm4hep::MCRecoClusterParticleAssociation       7  10793d25
MCTruthRecoLink                       edm4hep::MCRecoParticleAssociation              7  dc1423e8
MCTruthSiTracksLink                   edm4hep::MCRecoTrackParticleAssociation         2  d5219fb7
MUON                                  edm4hep::CalorimeterHit                         7  0f355ef3
OTrackerEndcapHits                    edm4hep::TrackerHitPlane                       12  e118ca47
OTrackerHits                          edm4hep::TrackerHitPlane                        2  3a241fef
OuterTrackerBarrelCollection          edm4hep::SimTrackerHit                          2  083cc03d
OuterTrackerBarrelHitsRelations       edm4hep::MCRecoTrackerHitPlaneAssociation       2  c42bbbee
OuterTrackerEndcapCollection          edm4hep::SimTrackerHit                         13  1450dff0
OuterTrackerEndcapHitsRelations       edm4hep::MCRecoTrackerHitPlaneAssociation      12  d1211017
PandoraClusters                       edm4hep::Cluster                                5  1ada7608
PandoraClusters_particleIDs           edm4hep::ParticleID                             0  d74d176a
PandoraPFOs                           edm4hep::ReconstructedParticle                  5  fa28d9be
PandoraPFOs_particleIDs               edm4hep::ParticleID                             0  2085e24d
PandoraStartVertices                  edm4hep::Vertex                                 5  baf37a08
PFOsFromJets                          edm4hep::ReconstructedParticle                  5  1979049e
PrimaryVertices                       edm4hep::Vertex                                 1  50874189
PrimaryVertices_RP                    edm4hep::ReconstructedParticle                  1  ee0a4817
PrimaryVertices_RP_particleIDs        edm4hep::ParticleID                             0  8c256fbd
RecoMCTruthLink                       edm4hep::MCRecoParticleAssociation              7  6b81837a
RefinedVertexJets                     edm4hep::ReconstructedParticle                  2  6a2b2496
RefinedVertexJets_particleIDs         edm4hep::ParticleID                             2  6e50f897
RefinedVertexJets_rel                 edm4hep::RecoParticleVertexAssociation          0  8dac6bb6
RefinedVertexJets_vtx                 edm4hep::Vertex                                 0  514d1442
RefinedVertexJets_vtx_RP              edm4hep::ReconstructedParticle                  0  f4bb2dda
RefinedVertexJets_vtx_RP_particleIDs  edm4hep::ParticleID                             0  fd8d9f37
RefinedVertices                       edm4hep::Vertex                                 0  392f5baf
RefinedVertices_RP                    edm4hep::ReconstructedParticle                  0  61544e17
RefinedVertices_RP_particleIDs        edm4hep::ParticleID                             0  bccf2a46
RelationCaloHit                       edm4hep::MCRecoCaloAssociation                859  603a5016
RelationMuonHit                       edm4hep::MCRecoCaloAssociation                  7  df24625a
SelectedPandoraPFOs                   edm4hep::ReconstructedParticle                  5  5d4278bb
SiTracks                              edm4hep::Track                                  2  82d755bd
SiTracks_Refitted                     edm4hep::Track                                  2  f4957afa
SiTracksCT                            edm4hep::Track                                  2  d094e81e
SiTracksMCTruthLink                   edm4hep::MCRecoTrackParticleAssociation         2  9a2f9ac5
TightSelectedPandoraPFOs              edm4hep::ReconstructedParticle                  4  5fa7cf93
VertexBarrelCollection                edm4hep::SimTrackerHit                         14  85e68310
VertexEndcapCollection                edm4hep::SimTrackerHit                          0  5add7f42
VertexJets                            edm4hep::ReconstructedParticle                  2  ca369aad
VertexJets_particleIDs                edm4hep::ParticleID                             2  9764e204
VXDEndcapTrackerHitRelations          edm4hep::MCRecoTrackerHitPlaneAssociation       0  bb4cff22
VXDEndcapTrackerHits                  edm4hep::TrackerHitPlane                        0  e9876f4e
VXDTrackerHitRelations                edm4hep::MCRecoTrackerHitPlaneAssociation      13  178c9330
VXDTrackerHits                        edm4hep::TrackerHitPlane                       13  db43b881
YokeBarrelCollection                  edm4hep::SimCalorimeterHit                      0  97bfada6
YokeBarrelCollectionContributions     edm4hep::CaloHitContribution                    0  5e0fdff3
YokeEndcapCollection                  edm4hep::SimCalorimeterHit                      7  60209550
YokeEndcapCollectionContributions     edm4hep::CaloHitContribution                   13  68f65222

Parameters:
Name               Type           Elements
-----------------  -----------  ----------
alphaQCD           std::string           1
alphaQED           std::string           1
event_scale        std::string           1
EventWeights       double                1
GenCrossSection    std::string           1
GenPdfInfo         std::string           1
signal_process_id  std::string           1


```
</details>






