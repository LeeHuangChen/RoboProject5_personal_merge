
try:
    from ompl import base as ob
    from ompl import geometric as og
    from ompl import util as ou
except:
    # if the ompl module is not in the PYTHONPATH assume it is installed in a
    # subdirectory of the parent directory called "py-bindings."
    from os.path import abspath, dirname, join
    import sys
    sys.path.insert(0, join(dirname(dirname(abspath(__file__))),'py-bindings'))
    
    from ompl import base as ob
    from ompl import geometric as og

from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import PeptideChain
from MMTK.Proteins import Protein
from MMTK.Collections import Collection
from MMTK.ForceFields import Amber94ForceField
#from MMTK.Visualization import view
from Scientific.Visualization import VMD; module = VMD

import numpy as np

global sc2
sc=None
sc2=None
numberOfResidues=None
universe=None
universe2=None
chain=None
chain2=None


class MyProjectionEvaluator(ob.ProjectionEvaluator):
    def __init__(self, space, cellSizes):
        super(MyProjectionEvaluator, self).__init__(space)
        #self.setCellSizes(cellSizes)

    def getDimension(self):
        return 2

    def project(self, state, projection):
        projection[0] = state[0].value
        projection[1] = state[1].value

def isStateValid(state):
    # Some arbitrary condition on the state (note that thanks to
    # dynamic type checking we can just call getX() and do not need
    # to convert state to an SE2State.)
    #return state[0].value < .6 or state[1].value <.6
    einit=universe.energy()
    j=0
    kb = 1.3806488e-23
    T=300
    for i in range(0,len(sc)):
    	sc2[i].phiAngle().setValue(state[j].value)
    	j=j+1

    	sc2[i].psiAngle().setValue(state[j].value)
    	j=j+1
    enext=universe2.energy()
    print "(",einit,",",enext,")"

    threshold=np.exp((einit-enext)/(kb*T))
    rng=ou.RNG()
    if(enext<=einit):
    	return True
    else:
    	rValue=rng.uniform01()
    	if(rValue<threshold):
    		return True
    	else:
    		return False

if __name__ == "__main__":
    configuration = PDBConfiguration('2YCC.pdb')
    chains = configuration.createPeptideChains()
    chain = chains[0]
    sc = chain[11:27]

    configuration2 = PDBConfiguration('2YCC.pdb')
    chains2 = configuration.createPeptideChains()
    chain2 = chains2[0]
    
    sc2 = chain2[11:27]



    numberOfResidues = len(sc)