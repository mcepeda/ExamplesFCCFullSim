# More complex example, not yet cleaned 
# Looks for MuTau / ETau combinations 
# Then checks the tau polatization using the hadronic tau 


import sys, os, math 
from array import array
import ROOT
from ROOT import TFile, TTree, TH1F, TH2F
import numpy as np
from podio import root_io
import edm4hep
from pathlib import Path
import ctypes

from modules import tauReco 
from modules import weightsPol
from modules import optimalVariabRho

import argparse
parser = argparse.ArgumentParser(description="Configure the analysis",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-f","--sample",default="ZTauTau_SMPol_7")
parser.add_argument("-o","--outfile",default="templates")
parser.add_argument("-d","--decay",default=2,type=int)
parser.add_argument("-p","--photonCut",default=0.1,type=float)
parser.add_argument("-R","--dRMax",default=0.4,type=float)
parser.add_argument("-n","--neutronCut",default=1,type=float)
parser.add_argument("-t","--tauCut",default=2,type=float)

args = parser.parse_args()
config = vars(args)
print(config)

selectDecay=args.decay

sample="Bhabha_test"
POL="BHABHA"

#sample="ZTauTau_SMPol_25Sept_MuonFix"
#POL="Test"

fileOutName=args.outfile+"_"+str(args.decay)+"_"+str(args.photonCut)+"_"+str(args.tauCut)+POL+".root"
#sample=args.sample
dRMax=args.dRMax
minP=args.photonCut
PNeutron=args.neutronCut
minTauCut=args.tauCut

# get all the files 
path="/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/"
file="out_reco_edm4hep_edm4hep"
filenames=[]
dir_path=path+"/"+sample
names = ROOT.std.vector('string')()
nfiles=len(os.listdir(dir_path))

nfiles=10 
#1000
badfiles=[-1] #11,15]

print (dir_path)
for i in range(1,nfiles+1):
    if (i in badfiles): continue # this one is broken? why?

#    filename=dir_path+"/{}".format(i)+"/"+file+".root"
    filename=dir_path+"/"+file+"_{}.root".format(i)
    print (filename)

    my_file = Path(filename)
    if my_file.is_file():
        root_file = tauReco.open_root_file(filename)
        if not root_file or root_file.IsZombie(): 
            continue 
        #print (my_file)
        filenames.append(filename)

print ("Read %d files" %len(filenames))
reader = root_io.Reader(filenames)

# collections to use 
genparts = "MCParticles"
#pfobjects="PandoraPFOs"
pfobjects ="TightSelectedPandoraPFOs"


treeName="outtree"
variabs=["genTauP","genMesonP","genPionP","genTauE","genTauM","genMesonE","genMesonM","genPionE","genPionM","genTauTheta","genTauPhi","genMesonTheta","genMesonPhi","genPionTheta","genPionPhi","gen_cos_theta","gen_cos_psi","gen_cos_beta","gen_w","weight_P1","weight_M1","genOmega","gen_cos_theta_tau","recoMesonP","recoPionP","recoMesonTheta","recoMesonPhi","recoPionTheta","recoPionPhi","recoMesonE","recoMesonM","recoPionE","recoPionM","cos_theta","cos_psi","cos_beta","omega","cos_theta_rho","genTauID","recoTauID","ZMass","GenZMass","GenZVisMass","beamE"]
outfile=ROOT.TFile(fileOutName,"RECREATE")
new_tree = ROOT.TTree(treeName,"processed variables")

branches = {}

for var in variabs:
        branches[var] = ctypes.c_double(0.0)  # Single double variable
        new_tree.Branch(var, ctypes.addressof(branches[var]), f"{var}/D")

# GEN
hMesonEOverBeamE=TH1F("MesonEOverBeamE","",50,0,1.5)
hMesonEOverBeamE_P1=TH1F("MesonEOverBeamE_P1","",50,0,1.5)
hMesonEOverBeamE_M1=TH1F("MesonEOverBeamE_M1","",50,0,1.5)

hMesonCosTheta=TH1F("MesonCosTheta","",201,-1,1)
hMesonCosTheta_P1=TH1F("MesonCosTheta_P1","",201,-1,1)
hMesonCosTheta_M1=TH1F("MesonCosTheta_M1","",201,-1,1)

hGenVisZMass=TH1F("GenVisZMass","",200,0,100)

# RECO

hRecoMesonEOverBeamE_ALL=TH1F("RecoMesonEOverBeamE_ALL","",50,0,1.5)
hRecoMesonEOverBeamE=TH1F("RecoMesonEOverBeamE_SIGNAL","",50,0,1.5)
hRecoMesonEOverBeamE_P1=TH1F("RecoMesonEOverBeamE_SIGNAL_P1","",50,0,1.5)
hRecoMesonEOverBeamE_M1=TH1F("RecoMesonEOverBeamE_SIGNAL_M1","",50,0,1.5)

hRecoMesonCosTheta_ALL=TH1F("RecoMesonCosTheta_ALL","",201,-1,1)
hRecoMesonCosTheta=TH1F("RecoMesonCosTheta_SIGNAL","",201,-1,1)
hRecoMesonCosTheta_P1=TH1F("RecoMesonCosTheta_SIGNAL_P1","",201,-1,1)
hRecoMesonCosTheta_M1=TH1F("RecoMesonCosTheta_SIGNAL_M1","",201,-1,1)

hRecoZMass=TH1F("RecoZMass","",200,0,100)

hRecoMeson_X=TH1F("RecoMeson_X_SIGNAL","",50,-1.5,1.5)
hRecoMeson_X_P1=TH1F("RecoMeson_X_SIGNAL_P1","",50,-1.5,1.5)
hRecoMeson_X_M1=TH1F("RecoMeson_X_SIGNAL_M1","",50,-1.5,1.5)
hRecoMeson_X_BG=TH1F("RecoMeson_X_BG","",50,-1.5,1.5)
hRecoMeson_X_ALL=TH1F("RecoMeson_X_ALL","",50,-1.5,1.5)


# BG checks
hMesonEOverBeamE_BG=TH1F("MesonEOverBeamE_BG","",50,0,1.5)
hMesonEOverBeamE_BGMuon=TH1F("MesonEOverBeamE_BGMuon","",50,0,1.5)
hMesonEOverBeamE_BGEle=TH1F("MesonEOverBeamE_BGEle","",50,0,1.5)
hMesonEOverBeamE_BGPion=TH1F("MesonEOverBeamE_BGPion","",50,0,1.5)
hMesonEOverBeamE_BGRho=TH1F("MesonEOverBeamE_BGRho","",50,0,1.5)
hMesonEOverBeamE_BGA1=TH1F("MesonEOverBeamE_BGA1","",50,0,1.5)
hMesonEOverBeamE_BGOther=TH1F("MesonEOverBeamE_BGOther","",50,0,1.5)

hRecoMesonEOverBeamE_BG=TH1F("RecoMesonEOverBeamE_BG","",50,0,1.5)
hRecoMesonEOverBeamE_BGMuon=TH1F("RecoMesonEOverBeamE_BGMuon","",50,0,1.5)
hRecoMesonEOverBeamE_BGEle=TH1F("RecoMesonEOverBeamE_BGEle","",50,0,1.5)
hRecoMesonEOverBeamE_BGPion=TH1F("RecoMesonEOverBeamE_BGPion","",50,0,1.5)
hRecoMesonEOverBeamE_BGRho=TH1F("RecoMesonEOverBeamE_BGRho","",50,0,1.5)
hRecoMesonEOverBeamE_BGA1=TH1F("RecoMesonEOverBeamE_BGA1","",50,0,1.5)
hRecoMesonEOverBeamE_BGOther=TH1F("RecoMesonEOverBeamE_BGOther","",50,0,1.5)

hMesonCosTheta_BG=TH1F("MesonCosTheta_BG","",201,-1,1)
hMesonCosTheta_BGMuon=TH1F("MesonCosTheta_BGMuon","",201,-1,1)
hMesonCosTheta_BGEle=TH1F("MesonCosTheta_BGEle","",201,-1,1)
hMesonCosTheta_BGPion=TH1F("MesonCosTheta_BGPion","",201,-1,1)
hMesonCosTheta_BGRho=TH1F("MesonCosTheta_BGRho","",201,-1,1)
hMesonCosTheta_BGA1=TH1F("MesonCosTheta_BGA1","",201,-1,1)
hMesonCosTheta_BGOther=TH1F("MesonCosTheta_BGOther","",201,-1,1)

hRecoMesonCosTheta_BG=TH1F("RecoMesonCosTheta_BG","",201,-1,1)
hRecoMesonCosTheta_BGMuon=TH1F("RecoMesonCosTheta_BGMuon","",201,-1,1)
hRecoMesonCosTheta_BGEle=TH1F("RecoMesonCosTheta_BGEle","",201,-1,1)
hRecoMesonCosTheta_BGPion=TH1F("RecoMesonCosTheta_BGPion","",201,-1,1)
hRecoMesonCosTheta_BGRho=TH1F("RecoMesonCosTheta_BGRho","",201,-1,1)
hRecoMesonCosTheta_BGA1=TH1F("RecoMesonCosTheta_BGA1","",201,-1,1)
hRecoMesonCosTheta_BGOther=TH1F("RecoMesonCosTheta_BGOther","",201,-1,1)

hRecoMesonE_BGMuon=TH1F("RecoMesonE_BGMuon","",50,0,50)
hAllPFMuonsE =TH1F("hAllPFMuonsE","",50,0,50)
hAllPFElectronsE =TH1F("hAllPFElectronsE","",50,0,50)

hMesonType_BG = TH1F("MesonType_BG","Type of the BG taus (misId)",31,-15,15)
hRecoMesonType_BG = TH1F("RecoMesonType_BG","RecoType of the BG taus (misId)",31,-15,15)
hMesonType = TH1F("MesonType","Type",31,-15,15)
hRecoMesonType = TH1F("RecoMesonType","RecoType",31,-15,15)

hRecoMeson_BGMuon_PhiTheta = TH2F ("hRecoMeson_BGMuon_PhiTheta","",100,0,3.15,100,-3.15,3.15)
hRecoMeson_BGMuon_PtTheta = TH2F ("hRecoMeson_BGMuon_PtTheta","",100,0,3.15,100,0,50)

hRecoMeson_BGEle_PhiTheta = TH2F ("hRecoMeson_BGEle_PhiTheta","",100,0,3.15,100,-3.15,3.15)
hRecoMeson_BGEle_PtTheta = TH2F ("hRecoMeson_BGEle_PtTheta","",100,0,3.15,100,0,50)

hOmega=TH1F("Omega_SIGNAL","",100,-1,1)
hOmega_P1=TH1F("Omega_SIGNAL_P1","",100,-1,1)
hOmega_M1=TH1F("Omega_SIGNAL_M1","",100,-1,1)
hOmega_BG=TH1F("Omega_BG","",100,-1,1)
hOmega_BGMuon=TH1F("Omega_BGMuon","",100,-1,1)
hOmega_BGEle=TH1F("Omega_BGEle","",100,-1,1)
hOmega_BGPion=TH1F("Omega_BGPion","",100,-1,1)
hOmega_BGRho=TH1F("Omega_BGRho","",100,-1,1)
hOmega_BGA1=TH1F("Omega_BGA1","",100,-1,1)
hOmega_BGOther=TH1F("Omega_BGOther","",100,-1,1)
hOmega_ALL=TH1F("Omega_ALL","",100,-1,1)

hOmega_GEN = TH1F("Omega_GEN_SIGNAL","",100,-1,1)
hOmega_GEN_P1 = TH1F("Omega_GEN_SIGNAL_P1","",100,-1,1)
hOmega_GEN_M1 = TH1F("Omega_GEN_SIGNAL_M1","",100,-1,1)
hOmega_GEN_BG = TH1F("Omega_GEN_BG","",100,-1,1)
hOmega_GEN_BGMuon = TH1F("Omega_GEN_BGMuon","",100,-1,1)
hOmega_GEN_BGEle = TH1F("Omega_GEN_BGEle","",100,-1,1)
hOmega_GEN_BGPion = TH1F("Omega_GEN_BGPion","",100,-1,1)  
hOmega_GEN_BGRho = TH1F("Omega_GEN_BGRho","",100,-1,1)  
hOmega_GEN_BGA1 = TH1F("Omega_GEN_BGA1","",100,-1,1)
hOmega_GEN_BGOther = TH1F("Omega_GEN_BGOther","",100,-1,1)
hOmega_GEN_ALL = TH1F("Omega_GEN_ALL","",100,-1,1)

hOmegaCosTheta=TH2F("OmegaCosTheta_SIGNAL","",100,-1,1,100,-1,1)
hOmegaCosTheta_P1=TH2F("OmegaCosTheta_SIGNAL_P1","",100,-1,1,100,-1,1)
hOmegaCosTheta_M1=TH2F("OmegaCosTheta_SIGNAL_M1","",100,-1,1,100,-1,1)
hOmegaCosTheta_BG=TH2F("OmegaCosTheta_BG","",100,-1,1,100,-1,1)
hOmegaCosTheta_BGMuon=TH2F("OmegaCosTheta_BGMuon","",100,-1,1,100,-1,1)
hOmegaCosTheta_BGEle=TH2F("OmegaCosTheta_BGEle","",100,-1,1,100,-1,1)
hOmegaCosTheta_BGPion=TH2F("OmegaCosTheta_BGPion","",100,-1,1,100,-1,1)
hOmegaCosTheta_BGRho=TH2F("OmegaCosTheta_BGRho","",100,-1,1,100,-1,1)
hOmegaCosTheta_BGA1=TH2F("OmegaCosTheta_BGA1","",100,-1,1,100,-1,1)
hOmegaCosTheta_BGOther=TH2F("OmegaCosTheta_BGOther","",100,-1,1,100,-1,1)
hOmegaCosTheta_ALL=TH2F("OmegaCosTheta_ALL","",100,-1,1,100,-1,1)

hOmegaCosThetaTau_GEN = TH2F("OmegaCosThetaTau_GEN_SIGNAL","",100,-1,1,100,-1,1)
hOmegaCosThetaTau_GEN_P1 = TH2F("OmegaCosThetaTau_GEN_SIGNAL_P1","",100,-1,1,100,-1,1)
hOmegaCosThetaTau_GEN_M1 = TH2F("OmegaCosThetaTau_GEN_SIGNAL_M1","",100,-1,1,100,-1,1)
hOmegaCosThetaTau_GEN_BG = TH2F("OmegaCosThetaTau_GEN_BG","",100,-1,1,100,-1,1)
hOmegaCosThetaTau_GEN_BGMuon = TH2F("OmegaCosThetaTau_GEN_BGMuon","",100,-1,1,100,-1,1)
hOmegaCosThetaTau_GEN_BGEle = TH2F("OmegaCosThetaTau_GEN_BGEle","",100,-1,1,100,-1,1)
hOmegaCosThetaTau_GEN_BGPion = TH2F("OmegaCosThetaTau_GEN_BGPion","",100,-1,1,100,-1,1)
hOmegaCosThetaTau_GEN_BGRho = TH2F("OmegaCosThetaTau_GEN_BGRho","",100,-1,1,100,-1,1)
hOmegaCosThetaTau_GEN_BGA1 = TH2F("OmegaCosThetaTau_GEN_BGA1","",100,-1,1,100,-1,1)
hOmegaCosThetaTau_GEN_BGOther = TH2F("OmegaCosThetaTau_GEN_BGOther","",100,-1,1,100,-1,1)
hOmegaCosThetaTau_GEN_ALL = TH2F("OmegaCosThetaTau_GEN_ALL","",100,-1,1,100,-1,1)

hCosThetaTau_GEN_ALL=TH1F("CosThetaTau_GEN_ALL","",100,-1,1)
hCosThetaRho_GEN_ALL=TH1F("CosThetaRho_GEN_ALL","",100,-1,1)  
hCosThetaRho_ALL=TH1F("CosThetaRho_ALL","",100,-1,1)      

hCosTheta=TH1F("CosTheta_SIGNAL","",100,-1,1)
hCosTheta_P1=TH1F("CosTheta_SIGNAL_P1","",100,-1,1)
hCosTheta_M1=TH1F("CosTheta_SIGNAL_M1","",100,-1,1)
hCosTheta_BG=TH1F("CosTheta_BG","",100,-1,1)
hCosTheta_BGMuon=TH1F("CosTheta_BGMuon","",100,-1,1)
hCosTheta_BGEle=TH1F("CosTheta_BGEle","",100,-1,1)
hCosTheta_BGPion=TH1F("CosTheta_BGPion","",100,-1,1)
hCosTheta_BGRho=TH1F("CosTheta_BGRho","",100,-1,1)
hCosTheta_BGA1=TH1F("CosTheta_BGA1","",100,-1,1)
hCosTheta_BGOther=TH1F("CosTheta_BGOther","",100,-1,1)
hCosTheta_ALL=TH1F("CosTheta_ALL","",100,-1,1)

hCosTheta_GEN = TH1F("CosTheta_GEN_SIGNAL","",100,-1,1)
hCosTheta_GEN_P1 = TH1F("CosTheta_GEN_SIGNAL_P1","",100,-1,1)
hCosTheta_GEN_M1 = TH1F("CosTheta_GEN_SIGNAL_M1","",100,-1,1)
hCosTheta_GEN_BG = TH1F("CosTheta_GEN_BG","",100,-1,1)
hCosTheta_GEN_BGMuon = TH1F("CosTheta_GEN_BGMuon","",100,-1,1)
hCosTheta_GEN_BGEle = TH1F("CosTheta_GEN_BGEle","",100,-1,1)
hCosTheta_GEN_BGPion = TH1F("CosTheta_GEN_BGPion","",100,-1,1)
hCosTheta_GEN_BGRho = TH1F("CosTheta_GEN_BGRho","",100,-1,1)
hCosTheta_GEN_BGA1 = TH1F("CosTheta_GEN_BGA1","",100,-1,1)
hCosTheta_GEN_BGOther = TH1F("CosTheta_GEN_BGOther","",100,-1,1)
hCosTheta_GEN_ALL = TH1F("CosTheta_GEN_ALL","",100,-1,1)


hCosPsi=TH1F("CosPsi_SIGNAL","",100,-1,1)
hCosPsi_P1=TH1F("CosPsi_SIGNAL_P1","",100,-1,1)
hCosPsi_M1=TH1F("CosPsi_SIGNAL_M1","",100,-1,1)
hCosPsi_BG=TH1F("CosPsi_BG","",100,-1,1)
hCosPsi_BGMuon=TH1F("CosPsi_BGMuon","",100,-1,1)
hCosPsi_BGEle=TH1F("CosPsi_BGEle","",100,-1,1)
hCosPsi_BGPion=TH1F("CosPsi_BGPion","",100,-1,1)
hCosPsi_BGRho=TH1F("CosPsi_BGRho","",100,-1,1)
hCosPsi_BGA1=TH1F("CosPsi_BGA1","",100,-1,1)
hCosPsi_BGOther=TH1F("CosPsi_BGOther","",100,-1,1)
hCosPsi_ALL=TH1F("CosPsi_ALL","",100,-1,1)

hCosPsi_GEN = TH1F("CosPsi_GEN_SIGNAL","",100,-1,1)
hCosPsi_GEN_P1 = TH1F("CosPsi_GEN_SIGNAL_P1","",100,-1,1)
hCosPsi_GEN_M1 = TH1F("CosPsi_GEN_SIGNAL_M1","",100,-1,1)
hCosPsi_GEN_BG = TH1F("CosPsi_GEN_BG","",100,-1,1)
hCosPsi_GEN_BGMuon = TH1F("CosPsi_GEN_BGMuon","",100,-1,1)
hCosPsi_GEN_BGEle = TH1F("CosPsi_GEN_BGEle","",100,-1,1)
hCosPsi_GEN_BGPion = TH1F("CosPsi_GEN_BGPion","",100,-1,1)
hCosPsi_GEN_BGRho = TH1F("CosPsi_GEN_BGRho","",100,-1,1)
hCosPsi_GEN_BGA1 = TH1F("CosPsi_GEN_BGA1","",100,-1,1)
hCosPsi_GEN_BGOther = TH1F("CosPsi_GEN_BGOther","",100,-1,1)
hCosPsi_GEN_ALL = TH1F("CosPsi_GEN_ALL","",100,-1,1)




# normalize

# I am sure there is a function that does this but I cannot find it
#totalEvents=0
#for e in reader.get("events"):
#    totalEvents+=1

totalEvents=1#296800
lumi=70*1000 # discuss this with Michele. Prob missing a factor 2
xsecZtautau=1476.58 #pb , from https://fcc-physics-events.web.cern.ch/FCCee/delphes/winter2023/idea/ 
weight=1#xsecZtautau/totalEvents 

print ("Events? :", totalEvents," Xsec (pb) :",xsecZtautau," Weight : ",weight) 

totalEvents=0
selectedEvents=0

sumWeightsP1=0
sumWeightsM1=0
sumWeights=0

# run over all events 
for event in reader.get("events"):

    if totalEvents%10000==0:
       print (totalEvents)
    totalEvents+=1
    mc_particles = event.get( genparts )
    beamE=mc_particles[0].getEnergy()

    ## get GEN level info
    genTaus=tauReco.findAllGenTaus(mc_particles)
    nGenTaus=len(genTaus)

    GenZMass=0
    GenZVisMass=0

    if nGenTaus==2:
     ZGen=genTaus[0][3]+genTaus[1][3]
     ZVisGen=genTaus[0][0]+genTaus[1][0]
     GenZMass=ZGen.M()
     GenZVisMass=ZVisGen.M()


    #hGenVisZMass.Fill(ZVisGen.M(),weight)
    #hGenZMass.Fill(ZGen.M(),weight)

    ## get RECO level info
    pfos = event.get(pfobjects)
    unsorted_recoTaus= tauReco.findAllTaus(pfos,dRMax, minP,PNeutron)
    recoTaus= tauReco.sort_by_P(unsorted_recoTaus)
    nRecoTaus=len(recoTaus)

    countMuonsP10=0
    countElectronsP10=0
    if nRecoTaus<1:
         continue

    
    for pf in pfos:
        if (abs(pf.getPDG())==13):
            muonP4=ROOT.TLorentzVector()
            muonP4.SetXYZM(pf.getMomentum().x,pf.getMomentum().y,pf.getMomentum().z,pf.getMass())
#            hAllPFMuonsE.Fill(pf.getEnergy(),weight)
            if muonP4.P()>10:
                 countMuonsP10+=1
     #            ZP4+=muonP4
        if (abs(pf.getPDG())==11):
            electronP4=ROOT.TLorentzVector()
            electronP4.SetXYZM(pf.getMomentum().x,pf.getMomentum().y,pf.getMomentum().z,pf.getMass())
#            hAllPFElectronsE.Fill(pf.getEnergy(),weight)
            if electronP4.P()>10:
                 countElectronsP10+=1
    #             ZP4+=electronP4

    recoMesonP4=recoTaus[0][0]
    recoTauID=recoTaus[0][1]
    recoTauQ=recoTaus[0][2]
    recoTauConsts=recoTaus[0][5]
    recoPion=recoTauConsts[0]

    recoPionP4=ROOT.TLorentzVector()
    recoPionP4.SetXYZM(recoPion.getMomentum().x,recoPion.getMomentum().y,recoPion.getMomentum().z,recoPion.getMass())

   # ZP4=recoMesonP4

    #genRhoP4=genTaus[0][0]
    #genTauId=genTaus[0][1]
    #genTauQ=genTaus[0][2]
    #genTauP4=genTaus[0][3]
    #genTauConsts=genTaus[0][6]
    #genPion=genTauConsts[0] # always the first at gen level

    #genPionP4=ROOT.TLorentzVector()
    #genPionP4.SetXYZM(genPion.getMomentum().x,genPion.getMomentum().y,genPion.getMomentum().z,genPion.getMass())

        # Do some selection
    if recoMesonP4.P()< minTauCut:
       continue
    if abs(math.cos(recoMesonP4.Theta()))>0.95:
       continue
    if (countMuonsP10<1 and countElectronsP10<1):
         continue
    if recoMesonP4.Theta()>1.565 and recoMesonP4.Theta()<1.575 :  # stupid cut to remove high momentum muons at pi/2
         continue 

    ZMass=0 #ZP4.M()

#    if nRecoTaus<2:
#         continue
#    secondrecoMesonP4=recoTaus[1][0]
#    secondrecoTauID=recoTaus[1][1]
#    secondRecoTauQ=recoTaus[1][2]

#    ZReco=MesonP4+secondMesonP4
#    hRecoZMass.Fill(ZReco.M())

    # Do some selection
#    if secondrecoMesonP4.P()< 1:
#       continue
#    if abs(math.cos(secondrecoMesonP4.Theta()))>0.95:
#       continue

    if recoTauID==selectDecay:

      genIndex=-1
      closestDR=10 # find always one
      for g in range(0,nGenTaus):
          #if genTaus[g][1]!=1: continue
          genP4=genTaus[g][0]
          dR=tauReco.dRAngle(genP4,recoMesonP4)
          if dR<closestDR:
             closestDR=dR
             genIndex=g 
   
      if genIndex==-1:
        print ("why????")
        continue      

      genMesonP4=genTaus[genIndex][0]
      genRhoP4=genTaus[genIndex][0] # repeat, just make it work, then clean
      genTauID=genTaus[genIndex][1]
      genTauP4=genTaus[genIndex][3]
      genTauConst=genTaus[genIndex][6]
      genPion=genTauConst[0]        
      genPionP4=ROOT.TLorentzVector()
      genPionP4.SetXYZM(genPion.getMomentum().x,genPion.getMomentum().y,genPion.getMomentum().z,genPion.getMass())

      (gen_cos_theta,gen_cos_psi,gen_cos_beta,gen_w,weight_P1,weight_M1)=optimalVariabRho.wVariab(genTauP4,genMesonP4,genPionP4,beamE)
      (cos_theta,cos_psi,cos_beta,w)=optimalVariabRho.wVariabRECO(recoMesonP4,recoPionP4,beamE)
      gen_cos_theta_tau=math.cos(genTauP4.Theta())
      cos_theta_rho=math.cos(recoMesonP4.Theta())

      #(gen_cos_theta,gen_cos_psi,gen_cos_beta,gen_w,weight_P1,weight_M1)=optimalVariabRho.wVariab(genTauP4,genRhoP4,genPionP4,beamE)
      #(cos_theta,cos_psi,cos_beta,w)=optimalVariabRho.wVariabRECO(recoMesonP4,recoPionP4,beamE)
      #gen_cos_theta_tau=math.cos(genTauP4.Theta())


      if abs(cos_theta)==1: continue  # border cases in which the calculation of cosTheta failed
      if abs(cos_psi)==1: continue

      selectedEvents+=1

#variabs=["genTauP","genMesonP","genPionP","genTauTheta","genTauPhi","genMesonTheta","genMesonPhi","genPionTheta","genPionPhi","gen_cos_theta","gen_cos_psi","gen_cos_beta","gen_w","weight_P1","weight_M1","genOmega","gen_cos_theta_tau","recoMesonP","recoPionP","recoMesonTheta","recoMesonPhi","recoPionTheta","recoPionPhi","cos_theta","cos_psi","cos_beta","omega","cos_theta_rho","genTauType","recotauType","ZMass","GenZMass","GenZVisMass","beamE"]

      branches["genTauP"].value=genTauP4.P()
      branches["genTauTheta"].value=genTauP4.Theta()
      branches["genTauPhi"].value=genTauP4.Phi()
      branches["genTauE"].value=genTauP4.E()
      branches["genTauM"].value=genTauP4.M()
      branches["genMesonP"].value=genMesonP4.P()
      branches["genMesonTheta"].value=genMesonP4.Theta()
      branches["genMesonPhi"].value=genMesonP4.Phi()
      branches["genMesonE"].value=genMesonP4.E()
      branches["genMesonM"].value=genMesonP4.M()
      branches["genPionP"].value=genPionP4.P()
      branches["genPionTheta"].value=genPionP4.Theta()
      branches["genPionPhi"].value=genPionP4.Phi()
      branches["genPionE"].value=genPionP4.E()
      branches["genPionM"].value=genPionP4.M()
      branches["gen_cos_theta"].value=gen_cos_theta
      branches["gen_cos_psi"].value=gen_cos_psi
      branches["gen_cos_beta"].value=gen_cos_beta
      branches["genOmega"].value=gen_w 
      branches["weight_P1"].value=weight_P1
      branches["weight_M1"].value=weight_M1
      branches["genTauID"].value=genTauID
      branches["gen_cos_theta_tau"].value=gen_cos_theta_tau
      branches["recoMesonP"].value=recoMesonP4.P()
      branches["recoMesonTheta"].value=recoMesonP4.Theta()
      branches["recoMesonPhi"].value=recoMesonP4.Phi()
      branches["recoPionP"].value=recoPionP4.P()
      branches["recoPionTheta"].value=recoPionP4.Theta()
      branches["recoPionPhi"].value=recoPionP4.Phi()
      branches["recoMesonE"].value=recoMesonP4.E()
      branches["recoMesonM"].value=recoMesonP4.M()
      branches["recoPionE"].value=recoPionP4.E()
      branches["recoPionM"].value=recoPionP4.M()
      branches["cos_theta"].value=cos_theta
      branches["cos_psi"].value=cos_psi
      branches["cos_beta"].value=cos_beta
      branches["omega"].value=w
      branches["recoTauID"].value=recoTauID
      branches["cos_theta_rho"].value=cos_theta_rho
      branches["ZMass"].value=ZMass
      branches["GenZMass"].value=GenZMass
      branches["GenZVisMass"].value=GenZVisMass
      branches["beamE"].value=beamE

      new_tree.Fill()


      hRecoMesonEOverBeamE_ALL.Fill(recoMesonP4.E()/beamE,weight)
      hRecoMesonCosTheta_ALL.Fill(math.cos(recoMesonP4.Theta()),weight)

      hCosTheta_ALL.Fill(cos_theta,weight)
      hCosPsi_ALL.Fill(cos_psi,weight)

      hOmega_GEN_ALL.Fill(gen_w,weight)
      hCosTheta_GEN_ALL.Fill(gen_cos_theta,weight)
      hCosPsi_GEN_ALL.Fill(gen_cos_psi,weight)

      hCosThetaTau_GEN_ALL.Fill(math.cos(genTauP4.Theta()),weight)
      hCosThetaRho_GEN_ALL.Fill(math.cos(genRhoP4.Theta()),weight)
      hCosThetaRho_ALL.Fill(math.cos(recoMesonP4.Theta()),weight)


      hOmegaCosTheta_ALL.Fill(w,cos_theta_rho,weight)
      hOmegaCosThetaTau_GEN_ALL.Fill(gen_w,gen_cos_theta_tau,weight)
      

      x=2*recoMesonP4.E()/beamE-1
      hRecoMeson_X_ALL.Fill(x,weight)
   
      if genIndex!=-1:
      

       selectGEN=selectDecay
       if selectDecay==2:
          selectGEN=1

       #print("GEN %d %4.2f  %4.2f %4.2f %4.2f %4.2f" %(genTauID,gen_cos_theta,gen_cos_theta_tau,gen_cos_beta, gen_cos_psi,gen_w))
       #print("RECO %d %4.2f %4.2f %4.2f %4.2f %4.2f" %(recoTauID,cos_theta,cos_theta_rho, cos_beta, cos_psi,w))

       if genTauID==selectGEN:
        weight_P1 = weightsPol.newAtauRHO( genTauP4, genMesonP4 , beamE, genTauConst, genTauID, +1)
        weight_M1 = weightsPol.newAtauRHO( genTauP4, genMesonP4 , beamE, genTauConst, genTauID, -1)

        hMesonEOverBeamE.Fill(genMesonP4.E()/beamE,weight)
        hMesonEOverBeamE_P1.Fill(genMesonP4.E()/beamE,weight_P1*weight)
        hMesonEOverBeamE_M1.Fill(genMesonP4.E()/beamE,weight_M1*weight)

        hOmega.Fill(w,weight)
        hOmega_P1.Fill(w,weight*weight_P1)
        hOmega_M1.Fill(w,weight*weight_M1)

        hOmega_GEN.Fill(gen_w,weight)
        hOmega_GEN_P1.Fill(gen_w,weight*weight_P1)
        hOmega_GEN_M1.Fill(gen_w,weight*weight_M1)

        hOmegaCosTheta.Fill(w,cos_theta_rho,weight)
        hOmegaCosTheta_P1.Fill(w,cos_theta_rho,weight*weight_P1)
        hOmegaCosTheta_M1.Fill(w,cos_theta_rho,weight*weight_M1)

        hOmegaCosThetaTau_GEN.Fill(gen_w,gen_cos_theta_tau,weight)
        hOmegaCosThetaTau_GEN_P1.Fill(gen_w,gen_cos_theta_tau,weight*weight_P1)
        hOmegaCosThetaTau_GEN_M1.Fill(gen_w,gen_cos_theta_tau,weight*weight_M1)

        hCosTheta.Fill(cos_theta,weight)
        hCosTheta_M1.Fill(cos_theta,weight_M1*weight)
        hCosTheta_P1.Fill(cos_theta,weight_P1*weight)

        hCosTheta_GEN.Fill(gen_cos_theta,weight)
        hCosTheta_GEN_M1.Fill(gen_cos_theta,weight_M1*weight)
        hCosTheta_GEN_P1.Fill(gen_cos_theta,weight_P1*weight)

        hCosPsi.Fill(cos_psi,weight)
        hCosPsi_M1.Fill(cos_psi,weight_M1*weight)
        hCosPsi_P1.Fill(cos_psi,weight_P1*weight)

        hCosPsi_GEN.Fill(gen_cos_psi,weight)
        hCosPsi_GEN_M1.Fill(gen_cos_psi,weight_M1*weight)
        hCosPsi_GEN_P1.Fill(gen_cos_psi,weight_P1*weight)

        hMesonCosTheta.Fill(math.cos(genMesonP4.Theta()),weight)
        hMesonCosTheta_P1.Fill(math.cos(genMesonP4.Theta()),weight_P1*weight)
        hMesonCosTheta_M1.Fill(math.cos(genMesonP4.Theta()),weight_M1*weight)

        hRecoMesonEOverBeamE.Fill(recoMesonP4.E()/beamE,weight)
        hRecoMesonEOverBeamE_P1.Fill(recoMesonP4.E()/beamE,weight_P1*weight)
        hRecoMesonEOverBeamE_M1.Fill(recoMesonP4.E()/beamE,weight_M1*weight)

        hRecoMeson_X.Fill(x,weight)
        hRecoMeson_X_P1.Fill(x,weight*weight_P1)
        hRecoMeson_X_M1.Fill(x,weight*weight_M1)

        hRecoMesonCosTheta.Fill(math.cos(recoMesonP4.Theta()),weight)
        hRecoMesonCosTheta_P1.Fill(math.cos(recoMesonP4.Theta()),weight_P1*weight)
        hRecoMesonCosTheta_M1.Fill(math.cos(recoMesonP4.Theta()),weight_M1*weight)

        hMesonType.Fill(genTauID,weight)
        hRecoMesonType.Fill(recoTauID,weight)
      
       else:
        hRecoMesonEOverBeamE_BG.Fill(recoMesonP4.E()/beamE,weight)
        hRecoMesonCosTheta_BG.Fill(math.cos(recoMesonP4.Theta()),weight)
        hMesonEOverBeamE_BG.Fill(genMesonP4.E()/beamE,weight)
        hMesonCosTheta_BG.Fill(math.cos(genMesonP4.Theta()),weight)
        hRecoMeson_X_BG.Fill(x,weight)
        hMesonType_BG.Fill(genTauID,weight)
        hRecoMesonType_BG.Fill(recoTauID,weight)
        hOmega_BG.Fill(w,weight)
        hCosTheta_BG.Fill(cos_theta,weight)
        hCosPsi_BG.Fill(cos_psi,weight)
        hOmegaCosTheta_BG.Fill(w,cos_theta_rho,weight)
        hOmegaCosThetaTau_GEN_BG.Fill(gen_w,gen_cos_theta_tau,weight) 

        hOmega_GEN_BG.Fill(gen_w,weight)
        hCosTheta_GEN_BG.Fill(gen_cos_theta,weight)
        hCosPsi_GEN_BG.Fill(gen_cos_psi,weight)
      

        if genTauID==-13:
          hRecoMesonEOverBeamE_BGMuon.Fill(recoMesonP4.E()/beamE,weight)
          hMesonEOverBeamE_BGMuon.Fill(genMesonP4.E()/beamE,weight)
          hRecoMesonE_BGMuon.Fill(recoMesonP4.E(),weight)
          hRecoMesonCosTheta_BGMuon.Fill(math.cos(recoMesonP4.Theta()),weight)
          hMesonCosTheta_BGMuon.Fill(math.cos(genMesonP4.Theta()),weight)
          hRecoMeson_BGMuon_PhiTheta.Fill(recoMesonP4.Theta(),recoMesonP4.Phi(),weight)
          hRecoMeson_BGMuon_PtTheta.Fill(recoMesonP4.Theta(),recoMesonP4.Pt(),weight)

          hOmega_BGMuon.Fill(w,weight)
          hOmegaCosTheta_BGMuon.Fill(w,cos_theta_rho,weight)
          hCosTheta_BGMuon.Fill(cos_theta,weight)
          hCosPsi_BGMuon.Fill(cos_psi,weight)
          hOmegaCosThetaTau_GEN_BGMuon.Fill(gen_w,gen_cos_theta_tau,weight)
          hOmega_GEN_BGMuon.Fill(gen_w,weight)
          hCosTheta_GEN_BGMuon.Fill(gen_cos_theta,weight)
          hCosPsi_GEN_BGMuon.Fill(gen_cos_psi,weight)

        elif genTauID==-11:
          hRecoMesonEOverBeamE_BGEle.Fill(recoMesonP4.E()/beamE,weight)
          hMesonEOverBeamE_BGEle.Fill(genMesonP4.E()/beamE,weight)
          hRecoMesonCosTheta_BGEle.Fill(math.cos(recoMesonP4.Theta()),weight)
          hMesonCosTheta_BGEle.Fill(math.cos(genMesonP4.Theta()),weight)
          hRecoMeson_BGEle_PhiTheta.Fill(recoMesonP4.Theta(),recoMesonP4.Phi(),weight)
          hRecoMeson_BGEle_PtTheta.Fill(recoMesonP4.Theta(),recoMesonP4.Pt(),weight)

          hOmega_BGEle.Fill(w,weight)
          hOmegaCosTheta_BGEle.Fill(w,cos_theta_rho,weight)
          hOmegaCosThetaTau_GEN_BGEle.Fill(gen_w,gen_cos_theta_tau,weight)
          hCosTheta_BGEle.Fill(cos_theta,weight)
          hCosPsi_BGEle.Fill(cos_psi,weight)
          hOmega_GEN_BGEle.Fill(gen_w,weight)
          hCosTheta_GEN_BGEle.Fill(gen_cos_theta,weight)
          hCosPsi_GEN_BGEle.Fill(gen_cos_psi,weight)

        elif genTauID==0:
          hRecoMesonEOverBeamE_BGPion.Fill(recoMesonP4.E()/beamE,weight)
          hMesonEOverBeamE_BGPion.Fill(genMesonP4.E()/beamE,weight)
          hRecoMesonCosTheta_BGPion.Fill(math.cos(recoMesonP4.Theta()),weight)
          hMesonCosTheta_BGPion.Fill(math.cos(genMesonP4.Theta()),weight)

          hOmega_BGPion.Fill(w,weight)
          hOmegaCosTheta_BGPion.Fill(w,cos_theta_rho,weight)
          hOmegaCosThetaTau_GEN_BGPion.Fill(gen_w,gen_cos_theta_tau,weight)
          hCosTheta_BGPion.Fill(cos_theta,weight)
          hCosPsi_BGPion.Fill(cos_psi,weight)
          hOmega_GEN_BGPion.Fill(gen_w,weight)
          hCosTheta_GEN_BGPion.Fill(gen_cos_theta,weight)
          hCosPsi_GEN_BGPion.Fill(gen_cos_psi,weight)

        elif genTauID==1:
          hRecoMesonEOverBeamE_BGRho.Fill(recoMesonP4.E()/beamE,weight)
          hMesonEOverBeamE_BGRho.Fill(genMesonP4.E()/beamE,weight)
          hRecoMesonCosTheta_BGRho.Fill(math.cos(recoMesonP4.Theta()),weight)
          hMesonCosTheta_BGRho.Fill(math.cos(genMesonP4.Theta()),weight)

          hOmega_BGRho.Fill(w,weight)
          hOmegaCosTheta_BGRho.Fill(w,cos_theta_rho,weight)
          hOmegaCosThetaTau_GEN_BGRho.Fill(gen_w,gen_cos_theta_tau,weight) 
          hCosTheta_BGRho.Fill(cos_theta,weight)
          hCosPsi_BGRho.Fill(cos_psi,weight)
          hOmega_GEN_BGRho.Fill(gen_w,weight)
          hCosTheta_GEN_BGRho.Fill(gen_cos_theta,weight)
          hCosPsi_GEN_BGRho.Fill(gen_cos_psi,weight)

        elif genTauID==10:
          hRecoMesonEOverBeamE_BGA1.Fill(recoMesonP4.E()/beamE,weight)
          hMesonEOverBeamE_BGA1.Fill(genMesonP4.E()/beamE,weight)
          hRecoMesonCosTheta_BGA1.Fill(math.cos(recoMesonP4.Theta()),weight)
          hMesonCosTheta_BGA1.Fill(math.cos(genMesonP4.Theta()),weight)

          hOmega_BGA1.Fill(w,weight)
          hOmegaCosTheta_BGA1.Fill(w,cos_theta_rho,weight)
          hOmegaCosThetaTau_GEN_BGA1.Fill(gen_w,gen_cos_theta_tau,weight)
          hCosTheta_BGA1.Fill(cos_theta,weight)
          hCosPsi_BGA1.Fill(cos_psi,weight)
          hOmega_GEN_BGA1.Fill(gen_w,weight)
          hCosTheta_GEN_BGA1.Fill(gen_cos_theta,weight)
          hCosPsi_GEN_BGA1.Fill(gen_cos_psi,weight)

        else:
          hRecoMesonEOverBeamE_BGOther.Fill(recoMesonP4.E()/beamE,weight)
          hMesonEOverBeamE_BGOther.Fill(genMesonP4.E()/beamE,weight)
          hRecoMesonCosTheta_BGOther.Fill(math.cos(recoMesonP4.Theta()),weight)
          hMesonCosTheta_BGOther.Fill(math.cos(genMesonP4.Theta()),weight)

          hOmega_BGOther.Fill(w,weight)
          hOmegaCosTheta_BGOther.Fill(w,cos_theta_rho,weight)
          hOmegaCosThetaTau_GEN_BGOther.Fill(gen_w,gen_cos_theta_tau,weight)
          hCosTheta_BGOther.Fill(cos_theta,weight)
          hCosPsi_BGOther.Fill(cos_psi,weight)
          hOmega_GEN_BGOther.Fill(gen_w,weight)
          hCosTheta_GEN_BGOther.Fill(gen_cos_theta,weight)
          hCosPsi_GEN_BGOther.Fill(gen_cos_psi,weight)


      else:
        hRecoMesonEOverBeamE_BG.Fill(recoMesonP4.E()/beamE,weight)
        hRecoMesonCosTheta_BG.Fill(math.cos(recoMesonP4.Theta()),weight)
        hMesonEOverBeamE_BG.Fill(0,weight)
        hMesonType_BG.Fill(-3,weight)
        hRecoMesonType_BG.Fill(recoTauID,weight) 

        hOmega_BG.Fill(w,weight)
        hOmegaCosTheta_BG.Fill(w,cos_theta_rho,weight)
        hOmegaCosThetaTau_GEN_BG.Fill(gen_w,gen_cos_theta_tau,weight)
        hCosTheta_BG.Fill(cos_theta,weight)
        hCosPsi_BG.Fill(cos_psi,weight)
        hOmega_GEN_BG.Fill(gen_w,weight)
        hCosTheta_GEN_BG.Fill(gen_cos_theta,weight)
        hCosPsi_GEN_BG.Fill(gen_cos_psi,weight)



print ("Run over ",totalEvents," selected ",selectedEvents)#," ->",selectedEvents/totalEvents)
print ("Weights?", sumWeights, sumWeightsP1,sumWeightsM1)
print ("Writing file ",fileOutName)

outfile.cd() # =ROOT.TFile(fileOutName,"RECREATE")

hGenVisZMass.Write()     

hMesonEOverBeamE.Write()
hMesonEOverBeamE_P1.Write()
hMesonEOverBeamE_M1.Write()
hMesonCosTheta.Write()
hMesonCosTheta_P1.Write()
hMesonCosTheta_M1.Write()

hRecoZMass.Write()

hRecoMesonEOverBeamE.Write()
hRecoMesonEOverBeamE_P1.Write()
hRecoMesonEOverBeamE_M1.Write()
hRecoMesonCosTheta.Write()
hRecoMesonCosTheta_P1.Write()
hRecoMesonCosTheta_M1.Write()

hMesonEOverBeamE_BG.Write()
hRecoMesonEOverBeamE_BG.Write()
hMesonEOverBeamE_BGMuon.Write()
hRecoMesonEOverBeamE_BGMuon.Write()
hMesonEOverBeamE_BGEle.Write()
hRecoMesonEOverBeamE_BGEle.Write()
hMesonEOverBeamE_BGPion.Write()
hRecoMesonEOverBeamE_BGPion.Write()
hMesonEOverBeamE_BGRho.Write()
hRecoMesonEOverBeamE_BGRho.Write()
hMesonEOverBeamE_BGA1.Write()
hRecoMesonEOverBeamE_BGA1.Write()
hMesonEOverBeamE_BGOther.Write()
hRecoMesonEOverBeamE_BGOther.Write()

hMesonCosTheta_BG.Write()
hRecoMesonCosTheta_BG.Write()
hMesonCosTheta_BGMuon.Write()
hRecoMesonCosTheta_BGMuon.Write()
hMesonCosTheta_BGEle.Write()
hRecoMesonCosTheta_BGEle.Write()
hMesonCosTheta_BGPion.Write()
hRecoMesonCosTheta_BGPion.Write()
hMesonCosTheta_BGRho.Write()
hRecoMesonCosTheta_BGRho.Write()
hMesonCosTheta_BGA1.Write()
hRecoMesonCosTheta_BGA1.Write()
hMesonCosTheta_BGOther.Write()
hRecoMesonCosTheta_BGOther.Write()

hMesonType.Write()
hRecoMesonType.Write()

hMesonType_BG.Write()
hRecoMesonType_BG.Write()

hRecoMesonEOverBeamE_ALL.Write()
hRecoMesonCosTheta_ALL.Write()

hRecoMesonE_BGMuon.Write()
hAllPFMuonsE.Write()
hAllPFElectronsE.Write()

hRecoMeson_BGMuon_PhiTheta.Write()
hRecoMeson_BGMuon_PtTheta.Write()
hRecoMeson_BGEle_PhiTheta.Write()
hRecoMeson_BGEle_PtTheta.Write()

hRecoMeson_X_ALL.Write()
hRecoMeson_X.Write()
hRecoMeson_X_P1.Write()
hRecoMeson_X_M1.Write()
hRecoMeson_X_BG.Write()


hOmega.Write()         
hOmega_P1.Write()       
hOmega_M1.Write()       
hOmega_BG.Write()       
hOmega_BGMuon.Write()   
hOmega_BGEle.Write()    
hOmega_BGPion.Write()   
hOmega_BGRho.Write()    
hOmega_BGA1.Write()     
hOmega_BGOther.Write()  
hOmega_ALL.Write()


hCosTheta.Write()       
hCosTheta_P1.Write()    
hCosTheta_M1.Write()    
hCosTheta_BG.Write()    
hCosTheta_BGMuon.Write()
hCosTheta_BGEle.Write() 
hCosTheta_BGPion.Write()
hCosTheta_BGRho.Write() 
hCosTheta_BGA1.Write()  
hCosTheta_BGOther.Write()
hCosTheta_ALL.Write()

hOmegaCosTheta.Write()
hOmegaCosTheta_P1.Write()       
hOmegaCosTheta_M1.Write()       
hOmegaCosTheta_BG.Write()       
hOmegaCosTheta_BGMuon.Write()   
hOmegaCosTheta_BGEle.Write()    
hOmegaCosTheta_BGPion.Write()   
hOmegaCosTheta_BGRho.Write()    
hOmegaCosTheta_BGA1.Write()     
hOmegaCosTheta_BGOther.Write()  
hOmegaCosTheta_ALL.Write()

hOmegaCosThetaTau_GEN.Write()
hOmegaCosThetaTau_GEN_P1.Write()
hOmegaCosThetaTau_GEN_M1.Write()
hOmegaCosThetaTau_GEN_BG.Write()
hOmegaCosThetaTau_GEN_BGMuon.Write()
hOmegaCosThetaTau_GEN_BGEle.Write()
hOmegaCosThetaTau_GEN_BGPion.Write()
hOmegaCosThetaTau_GEN_BGRho.Write()
hOmegaCosThetaTau_GEN_BGA1.Write()
hOmegaCosThetaTau_GEN_BGOther.Write()
hOmegaCosThetaTau_GEN_ALL.Write()

hCosPsi.Write()       
hCosPsi_P1.Write()     
hCosPsi_M1.Write()     
hCosPsi_BG.Write()     
hCosPsi_BGMuon.Write() 
hCosPsi_BGEle.Write()  
hCosPsi_BGPion.Write() 
hCosPsi_BGRho.Write()  
hCosPsi_BGA1.Write()   
hCosPsi_BGOther.Write()
hCosPsi_ALL.Write()

hOmega_GEN.Write()
hOmega_GEN_P1.Write()
hOmega_GEN_M1.Write()
hOmega_GEN_BG.Write()
hOmega_GEN_BGMuon.Write()
hOmega_GEN_BGEle.Write()
hOmega_GEN_BGPion.Write()
hOmega_GEN_BGRho.Write()
hOmega_GEN_BGA1.Write()
hOmega_GEN_BGOther.Write()
hOmega_GEN_ALL.Write()

hCosTheta_GEN.Write()
hCosTheta_GEN_P1.Write()
hCosTheta_GEN_M1.Write()
hCosTheta_GEN_BG.Write()
hCosTheta_GEN_BGMuon.Write()
hCosTheta_GEN_BGEle.Write()
hCosTheta_GEN_BGPion.Write()
hCosTheta_GEN_BGRho.Write()
hCosTheta_GEN_BGA1.Write()
hCosTheta_GEN_BGOther.Write()
hCosTheta_GEN_ALL.Write()


hCosThetaTau_GEN_ALL.Write()
hCosThetaRho_GEN_ALL.Write()
hCosThetaRho_ALL.Write()

hCosPsi_GEN.Write()
hCosPsi_GEN_P1.Write()
hCosPsi_GEN_M1.Write()
hCosPsi_GEN_BG.Write()
hCosPsi_GEN_BGMuon.Write()
hCosPsi_GEN_BGEle.Write()
hCosPsi_GEN_BGPion.Write()
hCosPsi_GEN_BGRho.Write()
hCosPsi_GEN_BGA1.Write()
hCosPsi_GEN_BGOther.Write()
hCosPsi_GEN_ALL.Write()

new_tree.Write()

outfile.Close()


