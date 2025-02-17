import sys
import math
import ROOT
from array import array
from podio import root_io
import edm4hep

from modules import myutils

# Check a generator level tau candidate, find the decay, 
# and compute visible (meson) variables 
def visTauGen(candTau):
        countPionsTauGen=0
        countPi0TauGen=0
        countMuonDecay=0
        countElectronDecay=0
        countOther=0

        genTauP4=ROOT.TLorentzVector()
        genTauP4.SetXYZM(candTau.getMomentum().x,candTau.getMomentum().y,candTau.getMomentum().z,candTau.getMass())

        visTauP4=ROOT.TLorentzVector()
        visTauP4.SetXYZM(0,0,0,0)
        chargeTau=0
        daughters=candTau.getDaughters()
        tauID=-1

        maxAngleConsts=0
        nConsts=0
        const={}

        # loop over daughter particles of the tau 
        for dau in daughters:
             dauP4=ROOT.TLorentzVector()
             dauP4.SetXYZM(dau.getMomentum().x,dau.getMomentum().y,dau.getMomentum().z,dau.getMass())
             dauPDG=abs(dau.getPDG())

             #print ('...dau',dauP4.P(),dauP4.Theta(),dauP4.Phi(),dau.getMass())

             # we want to compare the reco P4 to the 'visible' gen P4: skip neutrinos
             if (dauPDG==12 or dauPDG==14 or dauPDG==16):
                continue 

             # lepton decays 
             if dauPDG==13 :
                countMuonDecay+=1
                #continue # either filter here or at the analysis level
             if dauPDG==11 :
                countElectronDecay+=1               
                #continue # either filter here or at the analysis level

             # in this Pythia sample the tau decay directly goes to pi0/pi, without the rho/a1
             # to be checked in KKMC and Whizard...
             # if there was a rho, we would need an additional step  
             if dauPDG==211 or dauPDG==321 or dauPDG==323:   # kaons and pions paired together 
                countPionsTauGen+=1
             elif dauPDG==111 :
                countPi0TauGen+=1
             elif dau.getCharge()!=0 and (dauPDG!=11 and dauPDG!=13):
                print (dau.getPDG()) 
                countOther+=1
  
             # compute the angle of the constituents (cone size) for further studies 
             dR=myutils.dRAngle(genTauP4,dauP4)
             if maxAngleConsts<dR:
                   maxAngleConsts=dR

             const[nConsts]=dau
             nConsts+=1

             chargeTau+=dau.getCharge()
             visTauP4+=dauP4

        # now encode the ID in a int
        # this could be much more elegant, simple for now
        if countMuonDecay>0:
           tauID=-13
        elif countElectronDecay>0:
           tauID=-11
        elif countOther>0: # refinement: check what these are 
           tauID=-2
        elif abs(chargeTau)==1: # at gen level if this happens something failed
           if (countPionsTauGen==1):
                     tauID=countPi0TauGen
           elif (countPionsTauGen==3):
                     tauID=countPi0TauGen+10

        # return an object with the visible pt, ID, charge, and the true Pt 
        # a future step would be to define a class for the tau
        return (visTauP4,tauID,chargeTau,genTauP4,maxAngleConsts,nConsts,const)                 

# Reversed procedure for reconstructed pfos
# Starting from a pion, find particles in a cone around it, and 
# build the tau 
def buildTauFromPion(lead,allPfs,DRCone=1,minP=0,PNeutron=10):
        countPions=1
        countPhotons=0

        chargeTau=lead.getCharge()
        tauID=-1

        leadP4=ROOT.TLorentzVector()
        leadP4.SetXYZM(lead.getMomentum().x,lead.getMomentum().y,lead.getMomentum().z,lead.getMass())
        tauP4=ROOT.TLorentzVector()
        tauP4=leadP4

        maxConeTau=0
        const={}

        nConsts=1
        const[0]=lead      
        countNeutrons=0

#        print ('...lead',leadP4.P(),leadP4.Theta(),math.cos(leadP4.Theta())) # leadP4.Phi(),lead.getMass())

        for cand in allPfs:
             if (cand==lead): continue 
 
             candP4=ROOT.TLorentzVector()
             candP4.SetXYZM(cand.getMomentum().x,cand.getMomentum().y,cand.getMomentum().z,cand.getMass())
             candPDG=abs(cand.getPDG())

#             print ('...cand',candP4.P(),candP4.Theta(),candP4.Phi(),cand.getMass())
 
             # max angle? check backgrounds as well
             dR=myutils.dRAngle(candP4,leadP4) 
             if (dR>DRCone):  # cut to be tuned
                continue 
             # how low should we go in P?
             if (candP4.P()<minP):
                continue 

             # now check ID and clean
             if candPDG==11 or candPDG==13: 
                continue
             if (candPDG==2112 and candP4.P()>PNeutron): # Pandora FIXME: pion -> neutron misID 
                countNeutrons+=1
             if candPDG==211 :   # ignoring the difference between kaons and pions for now 
                countPions+=1
             if candPDG==22 :  # careful: here counting photons and not pi0s. Account for merged/lost photons.  
                countPhotons+=1
             else: 
                continue

             if maxConeTau<dR:
                  maxConeTau=dR

             chargeTau+=cand.getCharge()
             tauP4+=candP4
             const[nConsts]=cand
             nConsts+=1

        #print (tauP4.Pt(),chargeTau,countPions,countPhotons)

        # set the ID: only valid combinations can be a tau (charge, constituents compatible with
        # tau decay). can be refined in the future. 
        if abs(chargeTau)==1:
           if (countPions==1 and countNeutrons==0):
                 if countPhotons<=4:
                     tauID=countPhotons
                 if countPhotons>4:
                      tauID=5

           elif (countPions==3): 
                 if countPhotons<=2:
                     tauID=countPhotons+10 # more or less copied from the CMS convention for tauDecay
                 if countPhotons>2:
                     tauID=3+10 # capping the number of photons

           elif (countPions==1 and countNeutrons>0): # Future FIXME: Pandora pion->neutron misID issue 
                     tauID=15 

           #print (tauP4.P(),math.cos(tauP4.Theta()),chargeTau,countPions,countPhotons,tauID)
  
           # return an object with P4, ID, Charge, AngleMax, nConsts, constIdx 
           # should be a class in the future
           return (tauP4,tauID,chargeTau,maxConeTau,nConsts,const)

        else:
           tauP4.SetXYZM(0,0,0,0) # safety, always return an object 
           return (tauP4,-1,0,0,0,0) 
                 

# loop over all gen taus 
def findAllGenTaus(mc_particles):
    genTaus={}
    nGenTaus=0
    for particle in mc_particles:

        if abs(particle.getPDG()) != 15:
          continue
        # in the pythia sample we need to check the genStatus:
        # (in some events we have several copies of the tau)
        if particle.getGeneratorStatus()!=2:
          continue
#        print ("genTau!",particle.getGeneratorStatus())

        tauP4=ROOT.TLorentzVector()
        tauP4.SetXYZM(particle.getMomentum().x,particle.getMomentum().y,particle.getMomentum().z,particle.getMass())

        genTau=visTauGen(particle)
        visTauP4=genTau[0]
        genTauId=genTau[1]

        genTaus[nGenTaus]=genTau
        nGenTaus+=1

    return genTaus

# function to find all reco taus starting from PFO collection 
def findAllTaus(pfos,dRMax,minPt,PNeutron):

    taus={}
    nTaus=0

    for pf in pfos:
        if (abs(pf.getPDG())!=211): 
           continue 

        recoTau=buildTauFromPion(pf,pfos,dRMax,minPt,PNeutron)
        candTauP4=recoTau[0]
        candTauId=recoTau[1]
        candTauCharge=recoTau[2]

        # FIXME: this is very ugly, angular separation between taus to avoid duplicates
        # in these samples most of the events have 1 tau (and 1 prong)
        # in a real scenario this could be slow and we could veto events 
        duplicate=False
        for i in range(0,nTaus):
            if (myutils.dRAngle(candTauP4,taus[i][0])<0.05): duplicate=True
        if (duplicate==True): continue

        taus[nTaus]=recoTau
        nTaus+=1
        #print ("...",pf.getObjectID().index,candTauP4.Pt(),candTauP4.Phi(),candTauP4.Theta(),candTauId,candTauCharge)

    return taus

