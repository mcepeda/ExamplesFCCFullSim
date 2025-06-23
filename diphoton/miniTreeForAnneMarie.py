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
from modules import myutils

import argparse
parser = argparse.ArgumentParser(description="Configure the analysis",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Configuration, general
parser.add_argument("-f","--sample",default="ZTauTau_PolSM_March24_2M")
parser.add_argument("-o","--outfile",default="miniTree")

args = parser.parse_args()
config = vars(args)
print(config)

fileOutName=args.outfile+".root"
sample=args.sample

# get all the files
#path="/nfs/cms/cepeda/FCC/fullsim/" 
path="/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/"
file="out_reco_edm4hep_edm4hep"
filenames=[]
dir_path=path+"/"+sample
names = ROOT.std.vector('string')()
nfiles=len(os.listdir(dir_path))

nfiles=10# 2000 
badfiles=[-1] #11,15]

print (dir_path)
for i in range(1,nfiles+1):
    if (i in badfiles): continue # this one is broken? why?
#    filename=dir_path+"/{}".format(i)+"/"+file+".root"
    filename=dir_path+"/"+file+"_{}.root".format(i)
    print (filename)

    my_file = Path(filename)
    if my_file.is_file():
        root_file = myutils.open_root_file(filename) 
        # not critical, just paranoia since I had job failures 
        if not root_file or root_file.IsZombie():
            continue
        #print (my_file)
        filenames.append(filename)

print ("Read %d files" %len(filenames))
reader = root_io.Reader(filenames)

# collections to use 
genparts = "MCParticles"
pfobjects="PandoraPFOs"
#pfobjects ="TightSelectedPandoraPFOs"

# Configuration of the tree
treeName="outtree"
variabsVec=["photon","genPhoton","genPi0"]
vecComponents=["P","E","Px","Py","Pz","M"]
variabs=["beamE","nPhotons","nGenPhotons","nGenPi0s","nGenTaus","nRecoTausHad"]
outfile=ROOT.TFile(fileOutName,"RECREATE")
new_tree = ROOT.TTree(treeName,"processed variables")

branches = {}

for var in variabs:
        branches[var] = ctypes.c_double(0.0)  # Single double variable
        new_tree.Branch(var, ctypes.addressof(branches[var]), f"{var}/D")

for var in variabsVec:
    for comp in vecComponents:
        branches[var+comp] = ROOT.std.vector('double')()
        new_tree.Branch(var+comp, branches[var+comp])


# Accounting
hEvents = TH1F("hEvents","hEvents",2,0,2)
hPFMuonsE =TH1F("hPFMuonsE","",50,0,50)
hPFElectronsE =TH1F("hPFElectronsE","",50,0,50)
hPFPhotonsE =TH1F("hPFPhotonsE","",50,0,50)
hGenPhotonsE =TH1F("hGenPhotonsE","",50,0,50)
hGenPi0sE =TH1F("hGenPi0sE","",50,0,50)


totalEvents=0
selectedEvents=0

# run over all events 
for event in reader.get("events"):

    if totalEvents%10000==0:
       print (totalEvents)

    # clean vector branches
    for var in variabsVec:
      for comp in vecComponents:
         branches[var+comp].clear()

    # get event info
    totalEvents+=1
    mc_particles = event.get( genparts )
    beamE=mc_particles[0].getEnergy()
    pfos = event.get(pfobjects)

    ## get GEN level info
 
    # tau generator info (not needed, kept just in case)
    genTaus=tauReco.findAllGenTaus(mc_particles)
    nGenTaus=len(genTaus)
    #gentau1_visP4=genTaus[0][0]
    #gentau2_visP4=genTaus[1][0]
    #gentau1_ID=genTaus[0][1]
    #gentau2_ID=genTaus[1][1]
    #gentau1_Q=genTaus[0][2]
    #gentau2_Q=genTaus[1][2]
    #gentau1_P4=genTaus[0][3]
    #gentau2_P4=genTaus[1][3]

    nGenPhotons=0
    nGenPi0s=0

    for mc in mc_particles:
        if (abs(mc.getPDG())==22):
          if mc.getEnergy()>0.1 and mc.getGeneratorStatus()==1:
            photonP4=ROOT.TLorentzVector()
            photonP4.SetXYZM(mc.getMomentum().x,mc.getMomentum().y,mc.getMomentum().z,mc.getMass())
            hGenPhotonsE.Fill(mc.getEnergy())
            branches["genPhotonE"].push_back(mc.getEnergy())
            branches["genPhotonP"].push_back(photonP4.P())
            branches["genPhotonPx"].push_back(mc.getMomentum().x)
            branches["genPhotonPy"].push_back(mc.getMomentum().y)
            branches["genPhotonPz"].push_back(mc.getMomentum().z)
            branches["genPhotonM"].push_back(mc.getMass())

            nGenPhotons+=1
        if (abs(mc.getPDG())==111):
            pi0P4=ROOT.TLorentzVector()
            pi0P4.SetXYZM(mc.getMomentum().x,mc.getMomentum().y,mc.getMomentum().z,mc.getMass())
            hGenPi0sE.Fill(mc.getEnergy())
            branches["genPi0E"].push_back(mc.getEnergy())
            branches["genPi0P"].push_back(pi0P4.P())
            branches["genPi0Px"].push_back(mc.getMomentum().x)
            branches["genPi0Py"].push_back(mc.getMomentum().y)
            branches["genPi0Pz"].push_back(mc.getMomentum().z)
            branches["genPi0M"].push_back(mc.getMass())
            nGenPi0s+=1

    ## get RECO level info

    # build taus with DR=0.4, minP>0.1, neutron rejection (pandora bug) at 2 GeV)
    # not needed for diphoton study!
    unsorted_recoTaus= tauReco.findAllTaus(pfos,0.4, 0.1,2) 
    recoTaus= myutils.sort_by_P(unsorted_recoTaus)
    nRecoTausHad=len(recoTaus)

    nPhotons=0
    
    for pf in pfos:
        if (abs(pf.getPDG())==13):
            muonP4=ROOT.TLorentzVector()
            muonP4.SetXYZM(pf.getMomentum().x,pf.getMomentum().y,pf.getMomentum().z,pf.getMass())
            hPFMuonsE.Fill(pf.getEnergy())
        if (abs(pf.getPDG())==11):
            electronP4=ROOT.TLorentzVector()
            electronP4.SetXYZM(pf.getMomentum().x,pf.getMomentum().y,pf.getMomentum().z,pf.getMass())
            hPFElectronsE.Fill(pf.getEnergy())
        if (abs(pf.getPDG())==22):
          if pf.getEnergy()>0.1:
            photonP4=ROOT.TLorentzVector()
            photonP4.SetXYZM(pf.getMomentum().x,pf.getMomentum().y,pf.getMomentum().z,pf.getMass())
            hPFPhotonsE.Fill(pf.getEnergy())
            branches["photonE"].push_back(pf.getEnergy())
            branches["photonP"].push_back(photonP4.P())
            branches["photonPx"].push_back(pf.getMomentum().x)
            branches["photonPy"].push_back(pf.getMomentum().y)
            branches["photonPz"].push_back(pf.getMomentum().z)
            branches["photonM"].push_back(pf.getMomentum().z)
            nPhotons+=1

    # Some selection here?
    # if... 

    branches["nPhotons"].value=nPhotons
    branches["nGenPhotons"].value=nGenPhotons
    branches["nGenPi0s"].value=nGenPi0s
    branches["nGenTaus"].value=nGenTaus
    branches["nRecoTausHad"].value=nRecoTausHad
    branches["beamE"].value=beamE


    selectedEvents+=1
    new_tree.Fill()

hEvents.Fill(0,totalEvents)
hEvents.Fill(1,selectedEvents)

print ("Run over ",totalEvents," selected ",selectedEvents)#," ->",selectedEvents/totalEvents)
print ("Writing file ",fileOutName)

outfile.cd() # =ROOT.TFile(fileOutName,"RECREATE")

hEvents.Write()
hPFMuonsE.Write()
hPFElectronsE.Write()
hPFPhotonsE.Write()
hGenPhotonsE.Write()
hGenPi0sE.Write()

new_tree.Write()

outfile.Close()


