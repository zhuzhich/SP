from pyutilib.misc import import_file
from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.opt.base import SolverFactory
from pyomo.opt.parallel import SolverManagerFactory
from pyomo.opt.parallel.manager import solve_all_instances
import random
import math
import pdb

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
	print ("instance name = %s" %sub_insts[-1].name)
	for i in sub_insts[-1].sI:
		for r in sub_insts[-1].sR:
			#sub_insts[-1].pLT[i,r] = round((-math.log(random.uniform(0,1)))**\
			#						(1.0/sub_insts[-1].w_shape[i])*sub_insts[-1].w_scale[i])
			print "(%d,%d)=%d" %(i,r,sub_insts[-1].pLT[i,r]),
a = 1
b = 2