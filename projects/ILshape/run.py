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
"""
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
"""
def solve_callback(solver, mstr_inst):
    print "CB-Solve"
def cut_callback(solver, mstr_inst):
    print "CB-Cut"
def node_callback(solver, mstr_inst):
    print "CB-Node"
	

opt = SolverFactory('_cplex_direct')
opt.set_callback('cut-callback', cut_callback)
opt.set_callback('node-callback', node_callback)
opt.set_callback('solve-callback', solve_callback)

results = opt.solve(mstr_inst, tee=True)
#print results