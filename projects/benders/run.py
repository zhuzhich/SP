#!/usr/bin/env python

#  Implement Bender's decomposition multi cut in Birge book
#  Chapter 5, example 2, pp199 
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
import numpy
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
sub_insts.append(
    sb_mdl.create_instance(name="sub3", \
                           filename="sub3.dat"))
# initialize the solver manager.
solver_manager = SolverManagerFactory("serial")

# miscellaneous initialization.
for i in mstr_inst.Scen:
	mstr_inst.Min_Stage2_cost[i] = float("-Inf")

max_iterations = 50
flag_sita = [0 for i in range(3)]   
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
		mstr_inst.dSolution[s,i] = \
                inst.dual[inst.cY1]	
	#add master cut
	for s, inst in enumerate(sub_insts, 1):
		scen = mstr_inst.Scen[s]
		if flag_sita[s-1] == 1:
			continue
		cut = (mstr_inst.prob[s] *\
				mstr_inst.dSolution[s,i] * \
                mstr_inst.pH[s])
		cut +=	mstr_inst.prob[s] * \
				(mstr_inst.dSolution[s,i]* mstr_inst.pT) * \
				 (-mstr_inst.x)
		mstr_inst.Cut_Defn.add(mstr_inst.Min_Stage2_cost[s] >= cut)
		#
		scen = mstr_inst.Scen[s]
		Min_Stage2_cost = mstr_inst.prob[scen]*inst.oSub()
		print("Expected Stage2 cost ["+str(scen)+"]= " \
                               +str(round(Min_Stage2_cost, 4)))
		print("")			 
		newgap = round(mstr_inst.Min_Stage2_cost[scen].value - \
					Min_Stage2_cost, 6)
		newgap = abs(newgap)
		if newgap == 0:
        # get rid -0.0, which makes this script easier
        # to test against a baseline
			newgap = 0
		print("New gap= "+str(newgap)+"\n")

		if newgap <= 0.00001:
			flag_sita[s-1] = 1
			print("scenario["+str(s)+"]converged")

	if flag_sita.count(0) == 0:
		print("Multicut converged!!")
		break
# re-solve the master and update the subproblem inv1 values.
	solve_all_instances(solver_manager, 'cplex', [mstr_inst])
####
	for instance in sub_insts:
		instance.x = max(mstr_inst.x(),0)

else:
    # this gets executed when the loop above does not break
	print("Maximum Iterations Exceeded")

print("\nConverged master solution values:")
print("x="+str(round(mstr_inst.x(), 4)))
for p in mstr_inst.Scen:
	print("objective valule ["+str(p)+"]"+str(round(mstr_inst.Min_Stage2_cost[p].value, 4)))

