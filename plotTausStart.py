import sys, os, math 
from array import array
import ROOT
from ROOT import TFile, TTree, TH1F, TH2F
import numpy as np
from podio import root_io
import edm4hep
from pathlib import Path

from modules import tauReco
from modules import myutils 

import argparse
parser = argparse.ArgumentParser(description="Configure the analysis",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-f","--sample",default="ZTauTau_SMPol_25Sept_MuonFix")
parser.add_argument("-o","--outfile",default="firstTest_")
parser.add_argument("-d","--decay",default=-777,type=int)
parser.add_argument("-p","--photonCut",default=0.1,type=float)
parser.add_argument("-R","--dRMax",default=0.4,type=float)
parser.add_argument("-n","--neutronCut",default=1,type=float)
parser.add_argument("-t","--test",default=True,type=bool)


args = parser.parse_args()
config = vars(args)
print(config)

# get configuration parameters
dRMax=args.dRMax
minP=args.photonCut
selectDecay=args.decay
fileOutName=args.outfile
PNeutron=args.neutronCut
selectDecay=args.decay
sample=args.sample
test=args.test


decayString="decay"+str(selectDecay)+"_"+str(dRMax)+"_"+str(args.photonCut)+"_"+str(PNeutron)
if selectDecay==-777:
    decayString="decayAll_"+str(dRMax)+"_"+str(args.photonCut)+"_"+str(PNeutron)
fileOutName=args.outfile+decayString+".root"

print ("=====================================")

# get all the files 
path="/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/"
file="out_reco_edm4hep_edm4hep"
filenames=[]
dir_path=path+"/"+sample
names = ROOT.std.vector('string')()
nfiles=len(os.listdir(dir_path))

nfiles=1000
if test==True:
   nfiles=5

print ("Reading files from %s" %dir_path)
for i in range(1,nfiles+1):
    filename=dir_path+"/"+file+"_{}.root".format(i)
    print (filename)
    my_file = Path(filename)
    if my_file.is_file():
        root_file = myutils.open_root_file(filename)
        if not root_file or root_file.IsZombie():
            continue
        filenames.append(filename)

print ("Read %d files" %len(filenames))
reader = root_io.Reader(filenames)

# collections to use 
genparts = "MCParticles"
pfobjects="PandoraPFOs"
#pfobjects ="TightSelectedPandoraPFOs"

# Defining a few histogram
hGenTauP=TH1F("histoGenTauP","",50,0,50)
hGenVisTauP=TH1F("histoGenTauVisP","",50,0,50)
hGenTauType=TH1F("histoGenTauType","",40,-20,20)
hGenVisTauMass=TH1F("histoGenTauVisMass","",200,0,3)
hGenTauQ=TH1F("histoGenTauQ","",3,-1.5,1.5)
hGenTauTheta=TH1F("histoGenTauTheta","",100,0,3.15)
hGenTauPhi=TH1F("histoGenTauPhi","",100,-3.2,3.2)

hRecoTauP=TH1F("histoRecoTauP","",50,0,50)
hRecoTauType=TH1F("histoRecoTauType","",40,-20,20)
hRecoTauMass=TH1F("histoRecoTauMass","",200,0,3)
hRecoTauQ=TH1F("histoRecoTauQ","",3,-1.5,1.5)
hRecoTauTheta=TH1F("histoRecoTauTheta","",100,0,3.15)
hRecoTauPhi=TH1F("histoRecoTauPhi","",100,-3.2,3.2)

hNTaus=TH1F("histoNTaus","",6,0,6)
hNGenTaus=TH1F("histoNGenTaus","",6,0,6)
hNGenTausHad=TH1F("histoNGenTausHad","",6,0,6)


print ("-------------------------------------")
print ("Start processing!")

countEvents=0
# run over all events 
for event in reader.get("events"):

    if countEvents%500==0:
        print ("... %d" %countEvents)
    countEvents+=1

    # get the constituents
    mc_particles = event.get( genparts )
    pfos = event.get(pfobjects)

    # build the generator level taus 
    genTaus=tauReco.findAllGenTaus(mc_particles)
    nGenTaus=len(genTaus)
    nGenTausHad=0 # we do not know this yet 

    # build the reconstructed level taus 
    recoTaus= tauReco.findAllTaus(pfos,dRMax, minP,PNeutron)
    nTaus=len(recoTaus)

    for i in range(0,nGenTaus):
          # read the tau information
          genVisTauP4=genTaus[i][0] # to do: find a clearer dictionary for this
          genTauId=genTaus[i][1]
          genTauQ=genTaus[i][2]
          genTauP4=genTaus[i][3]
          #genTauDR=genTaus[i][4]
          #genTauNConsts=genTaus[i][5]
          #genTauConsts=genTaus[i][6]

          # count only the hadronic decays to compare to reco: 
          if genTauId<0:
             nGenTausHad+=1 

          # in case you want to only plot a decay 
          if selectDecay!=-777 and selectDecay!=genTauId:
             continue 

          # maybe you want to add some cuts 
          #if genVisTauP4.P()<5: continue 
          #if abs(math.cos(genVisTauP4.Theta())>0.9): continue

          # check what is in one event:
          #print ("Gen",genTauP4.P(),genVisTauP4.P(),genVisTauP4.Theta(),genVisTauP4.Phi(),genTauId,genTauQ,genTauDR)

          hGenTauP.Fill(genTauP4.P())
          hGenVisTauP.Fill(genVisTauP4.P())
          hGenVisTauMass.Fill(genVisTauP4.M())
          hGenTauType.Fill(genTauId)
          hGenTauQ.Fill(genTauQ)
          hGenTauTheta.Fill(genTauP4.Theta())
          hGenTauPhi.Fill(genTauP4.Phi())


    for j in range(0,nTaus):
          recoTauP4=recoTaus[j][0]
          recoTauId=recoTaus[j][1]
          recoTauQ=recoTaus[j][2]
          #recoTauDR=recoTaus[j][3]
          #recoTauNConsts=recoTaus[j][4]
          #recoTauConsts=recoTaus[j][5]

          # to make the code more economic we are checking gen and reco in parallel, but 
          # there is a difference in the DM labelling:
          # at reco level we count photons and at gen level pi0s: difference in the
          # decay mode (1 gen can be 1 or 2 reco, etc )
          recoDM=recoTauId
          if recoTauId==2:
             recoDM=1
          elif recoTauId>=3 and recoTauId<10:
             recoDM=3
          elif (recoTauId>=11 and recoTauId<15):
             recoDM=11

          if selectDecay!=-777 and selectDecay!=recoDM:
             continue

          # maybe you want to add some cuts 
          #if genVisTauP4.P()<5: continue 
          #if abs(math.cos(genVisTauP4.Theta())>0.9): continue

          # print information at the reco level
          #  print ("Reco?",recoTauP4.P(),recoTauP4.Theta(),recoTauP4.Phi(),recoTauId,recoTauQ,recoTauDR)

          hRecoTauType.Fill(recoTauId)
          hRecoTauP.Fill(recoTauP4.P())
          hRecoTauMass.Fill(recoTauP4.M())
          hRecoTauTheta.Fill(recoTauP4.Theta())
          hRecoTauPhi.Fill(recoTauP4.Phi())
          hRecoTauQ.Fill(recoTauQ)

    hNTaus.Fill(nTaus)
    hNGenTaus.Fill(nGenTaus)
    hNGenTausHad.Fill(nGenTausHad)


print ("-------------------------------------")
print ("Processed %d events" %countEvents)
print ("Plots saved in %s" %fileOutName)
print ("=====================================")


# save plots for later
outfile=ROOT.TFile(fileOutName,"RECREATE")

hGenTauP.Write()
hGenVisTauP.Write()
hGenTauType.Write()
hGenVisTauMass.Write()
hGenTauQ.Write()
hGenTauTheta.Write()
hGenTauPhi.Write()

hRecoTauP.Write()
hRecoTauType.Write()
hRecoTauMass.Write()
hRecoTauQ.Write()
hRecoTauTheta.Write()
hRecoTauPhi.Write()

hNTaus.Write()
hNGenTaus.Write()
hNGenTausHad.Write()

