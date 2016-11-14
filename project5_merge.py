#!/usr/bin/env python

######################################################################
# Software License Agreement (BSD License)
#
#  Copyright (c) 2010, Rice University
#  All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions
#  are met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials provided
#     with the distribution.
#   * Neither the name of the Rice University nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
#  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
#  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
#  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
#  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.
######################################################################

# Author: Mark Moll

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


sc="NULL"
numberOfResidues="NULL"
universe="NULL"

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
    	sc[i].phiAngle().setValue(state[j].value)
    	j=j+1

    	sc[i].psiAngle().setValue(state[j].value)
    	j=j+1
    enext=universe.energy()
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
    

def planWithSimpleSetup():
    # create an So2 state space
    n=numberOfResidues*2
    space = ob.CompoundStateSpace()

    for i in range(n):
    	space.addSubspace(ob.SO2StateSpace(),1)
    


    # create a simple setup object
    ss = og.SimpleSetup(space)
    ss.setStateValidityChecker(ob.StateValidityCheckerFn(isStateValid))

    start = ob.State(space)

    goal = ob.State(space)

    # Access all residues
    #goalAngles = []
    j=0
    for i in range(0,len(sc)): 
    	#print i,	
    	#goalAngles.append(sc[i].chiAngle().getValue())

    	start[j]=sc[i].phiAngle().getValue()
    	#sc[i].phiAngle().setValue(.5)
    	#start[j]=.5

    	j=j+1

    	start[j]=sc[i].psiAngle().getValue()
    	#sc[i].psiAngle().setValue(.5)
    	#start[j]=.5

    	j=j+1

    #print goalAngles
    global universe 
    universe = InfiniteUniverse(Amber94ForceField())
    protein = Protein(chain)
    universe.addObject(protein)
    energy = universe.energy()
    print energy


    

    
    #ss.setStartAndGoalStates(start, goal)
    ss.setStartState(start)
    ss.setGoal(ob.Goal(ss.getSpaceInformation()))

    planner=og.RRT(ss.getSpaceInformation())
    ss.setPlanner(planner)
    ss.setup()

    # this will automatically choose a default planner with
    # default parameters
    solved = ss.solve(10)

    if solved:
        # try to shorten the path
        ss.simplifySolution()
        # print the simplified path
        print(ss.getSolutionPath())
        
    pds=ob.PlannerDataStorage()
    plannerdata=ob.PlannerData(ss.getSpaceInformation())
    planner.getPlannerData(plannerdata)
    pds.store(plannerdata,"PlannerData.txt")
    


if __name__ == "__main__":
	configuration = PDBConfiguration('2YCC.pdb')
	chains = configuration.createPeptideChains()
	chain = chains[0]
	sc = chain[11:27]
	numberOfResidues = len(sc)
	planWithSimpleSetup()
