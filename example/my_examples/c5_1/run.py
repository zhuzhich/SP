#!/usr/bin/env python

#  Implement Bender's decomposition in Birge book
#  Chapter 5, example 1, pp184
#  Author: Zhicheng Zhu
#  Email: zzhu3@lamar.edu

# Bender's decomposition algorithm
##
##

# Python imports
from pyutilib.misc import import_file
from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.opt.base import SolverFactory
from pyomo.opt.parallel import SolverManagerFactory
from pyomo.opt.parallel.manager import solve_all_instances
import pdb

# import the master a
# initialize the master instance.
mstr_mdl = import_file("master.py").model
mstr_inst = mstr_mdl.create_instance("master.dat")

# initialize the sub-problem instances.
sb_mdl = import_file("sub.py").model
sub_insts = []
sub_insts.append(
    sb_mdl.create_instance(name="sub1", \
                           filename="sub1.dat"))
sub_insts.append(
    sb_mdl.create_instance(name="sub2", \
                           filename="sub2.dat"))
# initialize the solver manager.
solver_manager = SolverManagerFactory("serial")

# miscellaneous initialization.
mstr_inst.Min_Stage2_cost = float("-Inf")

gap = float("Inf")
max_iterations = 50

# main benders loop.
for i in range(1, max_iterations+1):
	print("\nIteration=%d"%(i))
	
	#solve the subproblem
	solve_all_instances(solver_manager, 'cplex', sub_insts)
	for instance in sub_insts:
		print("cost for scenario="+instance.name+" is"+\
				str(round(instance.oSub(),4)))
	print("")
	
	mstr_inst.CUTS.add(i)
	for s, inst in enumerate(sub_insts, 1):
		for t in mstr_inst.MAT:
			mstr_inst.dMatPro[t,s,i] = \
                inst.dual[inst.cMatPro[t]]
		for p in mstr_inst.sPRO:
			mstr_inst.dMaxPro[p,s,i] = \
                inst.dual[inst.cMaxPro[p]]
	#add master cut
	cut = sum((mstr_inst.prob[s] *\
				mstr_inst.dMatPro[t,s,i] * \
                mstr_inst.pHMatPro[t,s])
               for t in mstr_inst.MAT
               for s in mstr_inst.SCEN)
	cut += sum((mstr_inst.prob[s] *\
				mstr_inst.dMaxPro[t,s,i] * \
                mstr_inst.pHMaxPro[t,s])
               for t in mstr_inst.sPRO
               for s in mstr_inst.SCEN)
	cut +=	sum(mstr_inst.prob[s] * \
				(mstr_inst.dMatPro[t,s,i]* mstr_inst.pTMatPro[t,tx] + \
                 mstr_inst.dMaxPro[p,s,i]* mstr_inst.pTMaxPro[p,tx]) * \
				 (-mstr_inst.xMat[tx])
				 for t in mstr_inst.MAT
				 for tx in mstr_inst.MAT
				 for s in mstr_inst.SCEN)	
	mstr_inst.Cut_Defn.add(mstr_inst.Min_Stage2_cost >= cut)
	Min_Stage2_cost = 0;
	for s, inst in enumerate(sub_insts, 1):
		scen = mstr_inst.SCEN[s]
		Min_Stage2_cost += mstr_inst.prob[scen]*inst.oSub()
	print("Expected Stage2 cost= "+str(round(Min_Stage2_cost, 4)))
	print("")			 
	newgap = round(mstr_inst.Min_Stage2_cost.value - \
                   Min_Stage2_cost, 6)
	newgap = abs(newgap)
	if newgap == 0:
        # get rid -0.0, which makes this script easier
        # to test against a baseline
		newgap = 0
	print("New gap= "+str(newgap)+"\n")

	if newgap > 0.00001:
		gap = min(gap, newgap)
	else:
		print("Benders converged!")
		break
		
# re-solve the master and update the subproblem inv1 values.
	solve_all_instances(solver_manager, 'cplex', [mstr_inst])

	print("Master objective cost="+str(mstr_inst.oMaster()))
####
	for instance in sub_insts:
		for p in mstr_inst.MAT:
            # the master inventory values might be slightly
            # less than 0 (within tolerance); threshold here.
			instance.xMat[p] = max(mstr_inst.xMat[p](),mstr_inst.MinInput[p])
else:
    # this gets executed when the loop above does not break
    print("Maximum Iterations Exceeded")

print("\nConverged master solution values:")
for p in sorted(mstr_inst.MAT):
    print("x["+str(p)+"]="+str(round(mstr_inst.xMat[p](), 4)))
print("objective valule"+str(round(mstr_inst.oMaster(), 4)))
