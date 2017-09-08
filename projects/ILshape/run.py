from pyutilib.misc import import_file
from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.opt.base import SolverFactory
from pyomo.opt.parallel import SolverManagerFactory
from pyomo.opt.parallel.manager import solve_all_instances
import random
import math
import pdb
from debug import *

# import the master a
# initialize the master instance.
mstr_mdl = import_file("master.py").model
mstr_inst = mstr_mdl.create_instance("master.dat")

# initialize the sub-problem instances.
sb_mdl = import_file("sub.py").model
sub_insts = []
for s in mstr_inst.Scen:
	sub_insts.append(
		sb_mdl.create_instance(name="sub"+str(s), \
                           filename="sub.dat"))
	#print ("instance name = %s" %sub_insts[-1].name)
# initialize the solver manager.
solver_manager = SolverManagerFactory("serial")
# miscellaneous initialization.
for i in mstr_inst.Scen:
	mstr_inst.sita[i] = float("-Inf")

max_iterations = 50
#indicate each subproblem cut converged or not


# main benders loop.
for ii in range(1, max_iterations+1):
	flag_sita = [0 for i in range(0,mstr_inst.NUMSCEN())] 
	print("\nIteration=%d"%(ii))
	#solve the subproblem
	solve_all_instances(solver_manager, 'cplex', sub_insts)
	#start multi-cut
	for s, inst in enumerate(sub_insts, 1):
		subObj = round(inst.oSub(),4)
		print("obj. value for scenario "+str(s)+ " = "+inst.name+" is "+\
				str(subObj)+"("+str(mstr_inst.prob*subObj)+")")
		scen = mstr_inst.Scen[s]
		#print ("scen=%d" %scen)
		#if flag_sita[s-1] == 1:
		#	continue
		cut = 0
		#constraint g
		for i in inst.scg:
			cut += inst.dual[inst.cSg[i]]
		
		#constraint p
		for i in inst.sI:
			for t in inst.sT:
				cut += inst.dual[inst.cSp[i,t]]*(-1)
		
		#constraint q,r,s,t
		for i in inst.sI:
			for t in inst.sT:
				for r in inst.sR:
					cut += inst.dual[inst.cSq[i,t,r]]
					cut += inst.dual[inst.cSr[i,t,r]]
					cut += inst.dual[inst.cSs[i,t,r]]
					cut += inst.dual[inst.cSt[i,t,r]]
		#constraint u
		for t in inst.sT_0:
			cut += inst.dual[inst.cSu[t]]	
		#constraint v
		for i in inst.sI:
			for t in inst.sT_ex:
				for r in inst.sR:
					cut += inst.dual[inst.cSv[i,t,r]]
		cut = mstr_inst.prob*cut
		print "added cut"
		print cut,
		#constraint i
		for i in inst.sI:
			print ("+%f*x[%d]" %(mstr_inst.prob*inst.dual[inst.cSi[i]],i)),
			cut += mstr_inst.prob*inst.dual[inst.cSi[i]]*mstr_inst.x[i]
		print ""
		mstr_inst.Cut_Defn.add(mstr_inst.sita[s] >= cut)			
		Lbound = mstr_inst.prob*inst.oSub()
		print("Lower bound ["+str(scen)+"]= " \
                               +str(round(Lbound, 4)))
		print("Upper bound ["+str(scen)+"]"\
			+str(round(mstr_inst.sita[scen].value, 4)))
			
		newgap = round(mstr_inst.sita[scen].value - \
					Lbound, 6)
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
	print "solving master problem:"
	for i in mstr_inst.sI:
		print("x["+str(i)+"]="+str(round(mstr_inst.x[i](), 4)))
	print ("objective value=%f"  %mstr_inst.oMaster())
	for p in mstr_inst.Scen:
		print("objective valule ["+str(p)+"]"+str(round(mstr_inst.sita[p].value, 4)))

####

	for instance in sub_insts:
		for i in instance.sI:
			instance.x[i] = max(mstr_inst.x[i](),0)

else:
    # this gets executed when the loop above does not break
	print("Maximum Iterations Exceeded")

print("\nConverged master solution values:")
for i in mstr_inst.sI:
	print("x["+str(i)+"]="+str(round(mstr_inst.x[i](), 4)))
print ("objective value=%f" %mstr_inst.oMaster())
for p in mstr_inst.Scen:
	print("objective valule ["+str(p)+"]"+str(round(mstr_inst.sita[p].value, 4)))

  	
