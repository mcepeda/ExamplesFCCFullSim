import sys
import math
import ROOT
from array import array
from podio import root_io
import edm4hep

# I'm sure this exists already 
def dRAngle(p1,p2):
   dphi=p1.Phi()-p2.Phi()
   if (dphi>math.pi) : dphi=2*math.pi-dphi
   dtheta=p1.Theta()-p2.Theta()
   dR=math.sqrt(dtheta*dtheta+dphi*dphi)
   return dR

# trick to prevent broken files (should not be a problem at CIEMAT)
def open_root_file(file_path):
    try:
        # Suppress ROOT's default error messages to the terminal
        ROOT.gErrorIgnoreLevel = ROOT.kError

        # Attempt to open the ROOT file in "READ" mode without auto-recovery
        root_file = ROOT.TFile.Open(file_path, "READ")

        # Check if the file is a zombie
        if not root_file or root_file.IsZombie():
            raise IOError(f"Error: '{file_path}' is a zombie or could not be opened.")
        
        # Check if file is recoverable (potentially corrupted)
        if root_file.TestBit(ROOT.TFile.kRecovered):
            raise IOError(f"Error: '{file_path}' is corrupted and has been recovered.")
        
        #print(f"'{file_path}' opened successfully.")
        return root_file

    except Exception as e:
        print(f"Error: {e}")
        return None

# Fuction to sort by tau P
def sort_by_P(Tau):
    tau_with_P = []

    for i in range(0,len(Tau)):
        tau_with_P.append((Tau[i], Tau[i][0].P()))
    
    # Sort the list based on the P() value in descending order
    sorted_tau_with_P = sorted(tau_with_P, key=lambda x: x[1], reverse=True)
    
    # Extract only the sorted Tau[i] objects from the tuples
    sortedTau = [tau for tau, _ in sorted_tau_with_P]
   
    return sortedTau
