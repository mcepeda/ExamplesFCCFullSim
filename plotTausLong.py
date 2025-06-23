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
parser.add_argument("-f","--sample",default="ZTauTau_SMPol_March24_long")#"ZTauTau_SMPol_25Sept_MuonFix")
parser.add_argument("-o","--outfile",default="effis_")
parser.add_argument("-d","--decay",default=-777,type=int) # GEN 
parser.add_argument("-p","--photonCut",default=0.5,type=float)
parser.add_argument("-R","--dRMax",default=0.4,type=float)
parser.add_argument("-n","--neutronCut",default=1,type=float)
parser.add_argument("-t","--test",default=False,type=bool)


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

nfiles=200
if test==True:
   nfiles=20

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

print ("-------------------------------------")



# collections to use 
genparts = "MCParticles"
pfobjects="PandoraPFOs"
#pfobjects ="TightSelectedPandoraPFOs"

# Defining many histogram
hGenTauPt=TH1F("histoGenTauPt","",250,0,50)
hGenVisTauPt=TH1F("histoGenTauVisPt","",250,0,50)
hGenTauP=TH1F("histoGenTauP","",250,0,50)
hGenVisTauP=TH1F("histoGenTauVisP","",250,0,50)
hGenTauType=TH1F("histoGenTauType","",21,-1,20)
hGenVisTauMass=TH1F("histoGenTauVisMass","",500,0,10)
hGenTauQ=TH1F("histoGenTauQ","",3,-1.5,1.5)
hGenTauEta=TH1F("histoGenTauEta","",100,-5,5)
hGenTauTheta=TH1F("histoGenTauTheta","",100,0,3.15)
hGenTauDR=TH1F("histoGenTauDR","Angle of Tau Constituents",100,0,1)

hMatchedGenTauPt=TH1F("histoMatchedGenTauPt","",250,0,50)
hMatchedGenVisTauPt=TH1F("histoMatchedGenTauVisPt","",250,0,50)
hMatchedGenTauP=TH1F("histoMatchedGenTauP","",250,0,50)
hMatchedGenVisTauP=TH1F("histoMatchedGenTauVisP","",250,0,50)

hMatchedGenTauType=TH1F("histoMatchedGenTauType","",21,-1,20)
hMatchedGenVisTauMass=TH1F("histoMatchedGenTauVisMass","",500,0,10)
hMatchedGenTauQ=TH1F("histoMatchedGenTauQ","",3,-1.5,1.5)
hMatchedGenTauEta=TH1F("histoMatchedGenTauEta","",100,-5,5)
hMatchedGenTauTheta=TH1F("histoMatchedGenTauTheta","",100,0,3.15)

hMatchedGenTauDR=TH1F("histoMatchedGenTauDR","Angle of Tau Constituents",100,0,1)

hRecoTauPt=TH1F("histoRecoTauPt","",250,0,50)
hRecoTauP=TH1F("histoRecoTauP","",250,0,50)

hRecoTauMass=TH1F("histoRecoTauMass","",500,0,10)
hRecoTauType=TH1F("histoRecoTauType","",21,-1,20)
hRecoTauQ=TH1F("histoRecoTauQ","",3,-1.5,1.5)
hRecoTauEta=TH1F("histoRecoTauEta","",100,-5,5)
hRecoTauTheta=TH1F("histoRecoTauTheta","",100,0,3.15)
hRecoTauDR=TH1F("histoRecoTauDR","Angle of Tau Constituents",100,0,1)

hRecoConstP=TH1F("hRecoConstP","",500,0,50)
hRecoConstPhotonP=TH1F("hRecoConstPhotonP","",500,0,50)
hRecoConstPionP=TH1F("hRecoConstPionP","",500,0,50)

h2DRecoConstPType      =TH2F("hRecoConstPType","",500,0,50,15,0,15)
h2DRecoConstPhotonPType=TH2F("hRecoConstPhotonPType","",500,0,50,15,0,15)
h2DRecoConstPionPType  =TH2F("hRecoConstPionPType","",500,0,50,15,0,15)

hRecoConstPOverTauP=TH1F("hRecoConstPOverTauP","",100,0,2)
hRecoConstPhotonPOverTauP=TH1F("hRecoConstPhotonPOverTauP","",100,0,2)
hRecoConstPionPOverTauP=TH1F("hRecoConstPionPOverTauP","",100,0,2)

hRecoConstPion1P=TH1F("hRecoConstPion1P","",500,0,50)
hRecoConstPion2P=TH1F("hRecoConstPion2P","",500,0,50)
hRecoConstPion3P=TH1F("hRecoConstPion3P","",500,0,50)

hMatchedGenConstPion1P=TH1F("hMatchedGenConstPion1P","",500,0,50)
hMatchedGenConstPion2P=TH1F("hMatchedGenConstPion2P","",500,0,50)
hMatchedGenConstPion3P=TH1F("hMatchedGenConstPion3P","",500,0,50)

hGenConstPion1P=TH1F("hGenConstPion1P","",500,0,50)
hGenConstPion2P=TH1F("hGenConstPion2P","",500,0,50)
hGenConstPion3P=TH1F("hGenConstPion3P","",500,0,50)

hGenConstP=TH1F("hGenConstP","",500,0,50)
hGenConstPi0P=TH1F("hGenConstPi0P","",500,0,50)
hGenConstPionP=TH1F("hGenConstPionP","",500,0,50)
hGenConstPhotonP=TH1F("hGenConstPhotonP","",500,0,50)

h2DGenConstPType      =TH2F("hGenConstPType","",500,0,50,15,0,15)
h2DGenConstPi0PType=TH2F("hGenConstPi0PType","",500,0,50,15,0,15)
h2DGenConstPionPType  =TH2F("hGenConstPionPType","",500,0,50,15,0,15)
h2DGenConstPhotonPType=TH2F("hGenConstPhotonPType","",500,0,50,15,0,15)

hGenConstPOverTauP=TH1F("hGenConstPOverTauP","",100,0,2)
hGenConstPi0POverTauP=TH1F("hGenConstPi0POverTauP","",100,0,2)
hGenConstPionPOverTauP=TH1F("hGenConstPionPOverTauP","",100,0,2)
hGenConstPhotonPOverTauP=TH1F("hGenConstPhotonPOverTauP","",100,0,2)

hMatchedGenConstP=TH1F("hMatchedGenConstP","",500,0,50)
hMatchedGenConstPi0P=TH1F("hMatchedGenConstPi0P","",500,0,50)
hMatchedGenConstPionP=TH1F("hMatchedGenConstPionP","",500,0,50)
hMatchedGenConstPhotonP=TH1F("hMatchedGenConstPhotonP","",500,0,50)
hMatchedGenConstPOverTauP=TH1F("hMatchedGenConstPOverTauP","",100,0,2)
hMatchedGenConstPi0POverTauP=TH1F("hMatchedGenConstPi0POverTauP","",100,0,2)
hMatchedGenConstPionPOverTauP=TH1F("hMatchedGenConstPionPOverTauP","",100,0,2)
hMatchedGenConstPhotonPOverTauP=TH1F("hMatchedGenConstPhotonPOverTauP","",100,0,2)

h2DGenTauTypeMass=TH2F("histo2DGenTauTypeMass","",21,-1,20,500,0,10)
h2DRecoTauTypeMass=TH2F("histo2DRecoTauTypeMass","",21,-1,20,500,0,10)
hMatched2DGenTauTypeMass=TH2F("histoMatched2DGenTauTypeMass","",21,-1,20,500,0,10)

h2DGenTauDRType=TH2F("histo2DGenTauDRType","",100,0,1,21,-1,20)
h2DRecoTauDRType=TH2F("histo2DRecoTauDRType","",100,0,1,21,-1,20)
hMatched2DGenTauDRType=TH2F("histoMatched2DGenTauDRType","",100,0,1,21,-1,20)
h2DGenTauDRNConst=TH2F("histo2DGenTauDRNConst","",100,0,1,10,0,10)
h2DRecoTauDRNConst=TH2F("histo2DRecoTauDRNConst","",100,0,1,10,0,10)
hMatched2DGenTauDRNConst=TH2F("histoMatched2DGenTauDRNConst","",100,0,1,10,0,10)

h2DTauPt=TH2F("histo2DTauPt","",250,0,50,250,0,50)
h2DTauP=TH2F("histo2DTauP","",250,0,50,250,0,50)
h2DTauDR=TH2F("histo2DTauDR","",100,0,1,100,0,1)
h2DTauMass=TH2F("histo2DTauMass","",500,0,10,500,0,10)
h2DTauType=TH2F("histo2DTauType","",21,-1,20,21,-1,20)
h2DTauQ=TH2F("histo2DTauQ","",4,-2,2,4,-2,2)

hResTauPt=TH1F("histoResTauPt","",500,-1,1)
hResTauP=TH1F("histoResTauP","",500,-1,1)
hResTauMass=TH1F("histoResTauMass","",500,-1,1)

hNTaus=TH1F("histoNTaus","",6,0,6)
hNGenTaus=TH1F("histoNGenTaus","",6,0,6)
hNTausType=TH1F("histoNTausType","",6,0,6)
hNGenTausType=TH1F("histoNGenTausType","",6,0,6)

countEvents=0
# run over all events 
for event in reader.get("events"):

    if countEvents%1000==0:
        print (".... %d" %countEvents)
    countEvents+=1

    mc_particles = event.get( genparts )
    pfos = event.get(pfobjects)

    genTaus=tauReco.findAllGenTaus(mc_particles)
    nGenTaus=len(genTaus)

    recoTaus= tauReco.findAllTaus(pfos,dRMax, minP,PNeutron)
    nTaus=len(recoTaus)

    foundGen=False

    nGenTausType=0
    nTausType=0
    nGenTausHad=0

    for i in range(0,nGenTaus):
          genVisTauP4=genTaus[i][0] # to do: find a clearer dictionary for this
          genTauId=genTaus[i][1]
          genTauQ=genTaus[i][2]
          genTauP4=genTaus[i][3]
          genTauDR=genTaus[i][4]
          genTauNConsts=genTaus[i][5]
          genTauConsts=genTaus[i][6]

          # remove leptonic decays 
          if genTauId<0:
             continue

          nGenTausHad+=1

          # pick only a decay mode in particular if you want 
          if selectDecay!=-777 and selectDecay!=genTauId:
             continue 

          nGenTausType+=1
          foundGen=True

          if genVisTauP4.P()<5: continue 
          if abs(math.cos(genVisTauP4.Theta()))>0.9: continue

          #print ("Gen",genTauP4.P(),genVisTauP4.P(),genVisTauP4.Theta(),genVisTauP4.Phi(),genTauId,genTauQ,genTauDR,genTauNConsts)

          hGenTauPt.Fill(genTauP4.Pt())
          hGenVisTauPt.Fill(genVisTauP4.Pt())
          hGenTauP.Fill(genTauP4.P())
          hGenVisTauP.Fill(genVisTauP4.P())
          hGenVisTauMass.Fill(genVisTauP4.M())
          hGenTauType.Fill(genTauId)
          hGenTauQ.Fill(genTauQ)
          h2DGenTauTypeMass.Fill(genTauId,genVisTauP4.M())
          hGenTauEta.Fill(genTauP4.Eta())
          hGenTauTheta.Fill(genTauP4.Theta())

          hGenTauDR.Fill(genTauDR)
          h2DGenTauDRType.Fill(genTauDR,genTauId)
          h2DGenTauDRNConst.Fill(genTauDR,genTauNConsts)

          countPionsRun=0

          #print ("all GEN")
          # Look inside the generator level tau: check the constituents (decay products)

          for c in range(0,genTauNConsts):
             const=genTauConsts[c]
             constP4=ROOT.TLorentzVector()
             constP4.SetXYZM(const.getMomentum().x,const.getMomentum().y,const.getMomentum().z,const.getMass())

             hGenConstP.Fill(constP4.P())
             h2DGenConstPType.Fill(constP4.P(),genTauId)
             hGenConstPOverTauP.Fill(constP4.P()/genVisTauP4.P())

             # PDG ID for Pi0s == 111
             if const.getPDG()==111:
                hGenConstPi0P.Fill(constP4.P())
                h2DGenConstPi0PType.Fill(constP4.P(),genTauId)
                hGenConstPi0POverTauP.Fill(constP4.P()/genVisTauP4.P())
                # Pi0s decay to Photons 
                daughtersPi0=const.getDaughters()
                for dauPhoton in daughtersPi0:
                    dauPhotonP4=ROOT.TLorentzVector()
                    dauPhotonP4.SetXYZM(dauPhoton.getMomentum().x,dauPhoton.getMomentum().y,dauPhoton.getMomentum().z,dauPhoton.getMass())
                    hGenConstPhotonP.Fill(dauPhotonP4.P())
                    h2DGenConstPhotonPType.Fill(dauPhotonP4.P(),genTauId)
                    hGenConstPhotonPOverTauP.Fill(dauPhotonP4.P()/genVisTauP4.P())             

             elif abs(const.getPDG())==211:  # PDG ID for charged pions
                hGenConstPionP.Fill(constP4.P())
                h2DGenConstPionPType.Fill(constP4.P(),genTauId)
                hGenConstPionPOverTauP.Fill(constP4.P()/genVisTauP4.P())
                if countPionsRun==0:
                     hGenConstPion1P.Fill(constP4.P())
                     countPionsRun+=1
                elif countPionsRun==1:
                     hGenConstPion2P.Fill(constP4.P())
                     countPionsRun+=1
                elif countPionsRun==2:
                     hGenConstPion3P.Fill(constP4.P())
                     countPionsRun+=1


          findMatch=-1
          dRMatch=1

          # For each generator level tau, find the reconstructed tau that is closest:
          for j in range(0,nTaus):
            recoTauP4=recoTaus[j][0] # to do: find a clearer dictionary for this
            recoTauId=recoTaus[j][1]
            recoTauQ=recoTaus[j][2]
            recoTauDR=recoTaus[j][3]
            recoTauNConsts=recoTaus[j][4]
            recoTauConsts=recoTaus[j][5]

            # we want to study migrations: keep all the decays but count how many are good 
            # careful, at reco level we count photons and at gen level pi0s: difference in the
            # decay mode (1 gen can be 1,2 reco)

            recoDM=recoTauId
            if recoTauId==2:
              recoDM=1
            elif (recoTauId>=11 and recoTauId<15):
              recoDM=11
            elif recoTauId>=3 and recoTauId<10:
              recoDM=2

            if selectDecay!=-777 and selectDecay==recoDM:
                nTausType+=1

            # but remove at least the leptonic ones / failed ID
            if recoTauId<0:
              continue

            angleMatch=myutils.dRAngle(recoTauP4,genVisTauP4)

            # find closest
            if angleMatch<dRMatch:
               dRMatch=angleMatch
               findMatch=j
          
          # If you have not found it, continue: this is a efficiency loss 
          if findMatch==-1:
             continue 

          # now, get the kinematics of the matched reco tau
          recoTauP4=recoTaus[findMatch][0] 
          recoTauId=recoTaus[findMatch][1]
          recoTauQ=recoTaus[findMatch][2]
          recoTauDR=recoTaus[findMatch][3]
          recoTauNConsts=recoTaus[findMatch][4]
          recoTauConsts=recoTaus[findMatch][5]

#          print ("Reco?",recoTauP4.P(),recoTauId,recoTauQ,recoTauDR,recoTauNConsts)


          # Now that we have a matched (gen,reco) pair, more checks for efficiency and resolution

          countPionsRun=0
          #print ("Matched GEN!")
          # GEN: Look inside the tau, constituents: 
          for c in range(0,genTauNConsts):
             const=genTauConsts[c]
             constP4=ROOT.TLorentzVector()
             constP4.SetXYZM(const.getMomentum().x,const.getMomentum().y,const.getMomentum().z,const.getMass())

             hMatchedGenConstP.Fill(constP4.P())
             hMatchedGenConstPOverTauP.Fill(constP4.P()/genVisTauP4.P())

             if const.getPDG()==111:
                hMatchedGenConstPi0P.Fill(constP4.P())
                hMatchedGenConstPi0POverTauP.Fill(constP4.P()/genVisTauP4.P())
                daughtersPi0=const.getDaughters()
                for dauPhoton in daughtersPi0:
                    dauPhotonP4=ROOT.TLorentzVector()
                    dauPhotonP4.SetXYZM(dauPhoton.getMomentum().x,dauPhoton.getMomentum().y,dauPhoton.getMomentum().z,dauPhoton.getMass())
                    hMatchedGenConstPhotonP.Fill(dauPhotonP4.P())
                    hMatchedGenConstPhotonPOverTauP.Fill(dauPhotonP4.P()/genVisTauP4.P())             

             elif abs(const.getPDG())==211:
                hMatchedGenConstPionP.Fill(constP4.P())
                hMatchedGenConstPionPOverTauP.Fill(constP4.P()/genVisTauP4.P())
                if countPionsRun==0:
                     hMatchedGenConstPion1P.Fill(constP4.P())
                     countPionsRun+=1
                elif countPionsRun==1:
                     hMatchedGenConstPion2P.Fill(constP4.P())
                     countPionsRun+=1
                elif countPionsRun==2:
                     hMatchedGenConstPion3P.Fill(constP4.P())
                     countPionsRun+=1

          countPionsRun=0
          # RECO:  Look inside the tau, constituents: 
          for c in range(0,recoTauNConsts):
             const=recoTauConsts[c]
             constP4=ROOT.TLorentzVector()
             constP4.SetXYZM(const.getMomentum().x,const.getMomentum().y,const.getMomentum().z,const.getMass())
             hRecoConstP.Fill(constP4.P())
             h2DRecoConstPType.Fill(constP4.P(),recoTauId)
             hRecoConstPOverTauP.Fill(constP4.P()/recoTauP4.P())

             if const.getCharge()==0:
                hRecoConstPhotonP.Fill(constP4.P())
                h2DRecoConstPhotonPType.Fill(constP4.P(),recoTauId)
                hRecoConstPhotonPOverTauP.Fill(constP4.P()/recoTauP4.P())

             else:
                hRecoConstPionP.Fill(constP4.P())
                h2DRecoConstPionPType.Fill(constP4.P(),recoTauId)
                hRecoConstPionPOverTauP.Fill(constP4.P()/recoTauP4.P())
                if countPionsRun==0:
                     hRecoConstPion1P.Fill(constP4.P())
                     countPionsRun+=1
                elif countPionsRun==1:
                     hRecoConstPion2P.Fill(constP4.P())
                     countPionsRun+=1
                elif countPionsRun==2:
                     hRecoConstPion3P.Fill(constP4.P())
                     countPionsRun+=1


          # Many plots 

          hRecoTauType.Fill(recoTauId)
          hRecoTauPt.Fill(recoTauP4.Pt())
          hRecoTauP.Fill(recoTauP4.P())

          hRecoTauMass.Fill(recoTauP4.M())
          hRecoTauQ.Fill(recoTauQ)
          h2DRecoTauTypeMass.Fill(recoTauId,recoTauP4.M())
          hRecoTauEta.Fill(recoTauP4.Eta())
          hRecoTauTheta.Fill(recoTauP4.Theta()) 

          hRecoTauDR.Fill(recoTauDR)
          h2DRecoTauDRType.Fill(recoTauDR,recoTauId)
          h2DRecoTauDRNConst.Fill(recoTauDR,recoTauNConsts)

          hMatchedGenTauPt.Fill(genTauP4.Pt())
          hMatchedGenVisTauPt.Fill(genVisTauP4.Pt())
          hMatchedGenTauP.Fill(genTauP4.P())
          hMatchedGenVisTauP.Fill(genVisTauP4.P())

          hMatchedGenVisTauMass.Fill(genVisTauP4.M())
          hMatchedGenTauType.Fill(genTauId)
          hMatchedGenTauQ.Fill(genTauQ)
          hMatched2DGenTauTypeMass.Fill(genTauId,genVisTauP4.M())
          hMatchedGenTauEta.Fill(genTauP4.Eta())
          hMatchedGenTauTheta.Fill(genTauP4.Theta())

          hMatchedGenTauDR.Fill(genTauDR)
          hMatched2DGenTauDRType.Fill(genTauDR,genTauId)
          hMatched2DGenTauDRNConst.Fill(genTauDR,genTauNConsts)

          h2DTauPt.Fill(recoTauP4.Pt(),genVisTauP4.Pt())
          h2DTauMass.Fill(recoTauP4.M(),genVisTauP4.M())
          h2DTauType.Fill(recoTauId,genTauId)
          h2DTauP.Fill(recoTauP4.P(),genVisTauP4.P())

          h2DTauDR.Fill(recoTauDR,genTauDR)

          # Resolution plots:
          if genVisTauP4.P()!=0:
             hResTauPt.Fill( (recoTauP4.Pt()-genVisTauP4.Pt())/genVisTauP4.Pt())
             hResTauMass.Fill( (recoTauP4.M()-genVisTauP4.M())/genVisTauP4.M())
             hResTauP.Fill( (recoTauP4.P()-genVisTauP4.P())/genVisTauP4.P())

    #print ("Taus???",nGenTaus,nTaus) 
    hNTaus.Fill(nTaus)
    hNTausType.Fill(nTausType)
    hNGenTausType.Fill(nGenTausType)
    hNGenTaus.Fill(nGenTausHad)


# Do efficiencies (divide matched gen by all gen)

hEffiGenTauPt=hMatchedGenTauPt.Clone()
hEffiGenTauPt.SetName("hEffiGenTauPt")
hEffiGenTauPt.Divide(hGenTauPt)

hEffiGenVisTauPt=hMatchedGenVisTauPt.Clone()
hEffiGenVisTauPt.SetName("hEffiGenVisTauPt")
hEffiGenVisTauPt.Divide(hGenVisTauPt)

hEffiGenTauP=hMatchedGenTauP.Clone()
hEffiGenTauP.SetName("hEffiGenTauP")
hEffiGenTauP.Divide(hGenTauP)

hEffiGenVisTauP=hMatchedGenVisTauP.Clone()
hEffiGenVisTauP.SetName("hEffiGenVisTauP")
hEffiGenVisTauP.Divide(hGenVisTauP)

hEffiGenVisTauMass=hMatchedGenVisTauMass.Clone()
hEffiGenVisTauMass.SetName("hEffiGenVisTauMass")
hEffiGenVisTauMass.Divide(hGenVisTauMass)

hEffiGenTauEta=hMatchedGenTauEta.Clone()
hEffiGenTauEta.SetName("hEffiGenTauEta")
hEffiGenTauEta.Divide(hGenTauEta)

hEffiGenTauTheta=hMatchedGenTauTheta.Clone()
hEffiGenTauTheta.SetName("hEffiGenTauTheta")
hEffiGenTauTheta.Divide(hGenTauTheta)

hEffiGenTauType=hMatchedGenTauType.Clone()
hEffiGenTauType.SetName("hEffiGenTauType")
hEffiGenTauType.Divide(hGenTauType)

print ("-------------------------------------")
print ("Processed %d events" %countEvents)
print ("Plots saved in %s" %fileOutName)
print ("=====================================")


# save plots for later
outfile=ROOT.TFile(fileOutName,"RECREATE")

hGenTauPt.Write()
hGenVisTauPt.Write()
hGenTauP.Write()
hGenVisTauP.Write()
hGenTauType.Write()
hGenVisTauMass.Write()
hGenTauQ.Write()
hGenTauEta.Write()
hGenTauTheta.Write()

hGenTauDR.Write()

hMatchedGenTauPt.Write()
hMatchedGenVisTauPt.Write()
hMatchedGenTauP.Write()
hMatchedGenVisTauP.Write()

hMatchedGenTauType.Write()
hMatchedGenVisTauMass.Write()
hMatchedGenTauQ.Write()
hMatchedGenTauEta.Write()
hMatchedGenTauTheta.Write()

hMatchedGenTauDR.Write()


hEffiGenTauPt.Write()
hEffiGenVisTauPt.Write()
hEffiGenTauP.Write()
hEffiGenVisTauP.Write()

hEffiGenTauType.Write()
hEffiGenVisTauMass.Write()
hEffiGenTauEta.Write()
hEffiGenTauTheta.Write()


hRecoTauPt.Write()
hRecoTauP.Write()

hRecoTauType.Write()
hRecoTauMass.Write()
hRecoTauQ.Write()
hRecoTauEta.Write()
hRecoTauTheta.Write()

hRecoTauDR.Write()

h2DTauPt.Write()
h2DTauP.Write()

h2DTauMass.Write()
h2DTauType.Write()

h2DTauDR.Write()

hResTauPt.Write()
hResTauP.Write()

hResTauMass.Write()

h2DGenTauTypeMass.Write()
h2DRecoTauTypeMass.Write()
hMatched2DGenTauTypeMass.Write()

h2DGenTauDRType.Write()
h2DRecoTauDRType.Write()
hMatched2DGenTauDRType.Write()
h2DGenTauDRNConst.Write()
h2DRecoTauDRNConst.Write()
hMatched2DGenTauDRNConst.Write()

hNTaus.Write()
hNGenTaus.Write()

hNTausType.Write()
hNGenTausType.Write()

hRecoConstP.Write()            
hRecoConstPhotonP.Write()      
hRecoConstPionP.Write()        

h2DRecoConstPType.Write()      
h2DRecoConstPhotonPType.Write()
h2DRecoConstPionPType.Write()  

hRecoConstPOverTauP.Write()       
hRecoConstPhotonPOverTauP.Write() 
hRecoConstPionPOverTauP.Write()   

hGenConstP.Write()
hGenConstPi0P.Write()
hGenConstPionP.Write()
hGenConstPhotonP.Write()

h2DGenConstPType.Write()
h2DGenConstPi0PType.Write()
h2DGenConstPionPType.Write()
h2DGenConstPhotonPType.Write()

hGenConstPOverTauP.Write()
hGenConstPi0POverTauP.Write()
hGenConstPionPOverTauP.Write()
hGenConstPhotonPOverTauP.Write()

hMatchedGenConstP.Write()
hMatchedGenConstPi0P.Write()
hMatchedGenConstPionP.Write()
hMatchedGenConstPhotonP.Write()

hMatchedGenConstPOverTauP.Write()
hMatchedGenConstPi0POverTauP.Write()
hMatchedGenConstPionPOverTauP.Write()
hMatchedGenConstPhotonPOverTauP.Write()

hMatchedGenConstPion1P.Write()
hGenConstPion1P.Write()
hRecoConstPion1P.Write()

hMatchedGenConstPion2P.Write()
hGenConstPion2P.Write()
hRecoConstPion2P.Write()

hMatchedGenConstPion3P.Write()
hGenConstPion3P.Write()
hRecoConstPion3P.Write()


outfile.Close() 
