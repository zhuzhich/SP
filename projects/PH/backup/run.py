from pyutilib.misc import import_file
from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.opt.base import SolverFactory
from pyomo.opt.parallel import SolverManagerFactory
from pyomo.opt.parallel.manager import solve_all_instances
import random
import math
import pdb
import time

# initialize the instances.
ef_mdl = import_file("ReferenceModel.py").model
ef_insts = [ef_mdl.create_instance(name="PH_instance", filename="scen.dat")]
solver_manager = SolverManagerFactory("serial")
start_time1 = time.clock()
solve_all_instances(solver_manager, 'cplex', ef_insts)
end_time1 = time.clock()
print ("solving time = ",str(end_time1-start_time1))
"""
for i in ef_insts[0].sI:
	for t in ef_insts[0].sT:
		print ("")
		for r in ef_insts[0].sR:
			print ("(%d,%d,%d)=(%d,%d)" %(i,t,r,ef_insts[0].u[i,t,r,1](),ef_insts[0].v[i,t,r,1]())),
"""			

for i in ef_insts[0].sI:
	print("x["+str(i)+"]="+str(round(ef_insts[0].x[i](), 4)))
print ("objective value=%f" %ef_insts[0].objCost())
  	
