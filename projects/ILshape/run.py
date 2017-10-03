#!/usr/bin/python

from __future__ import print_function

from math import fabs
import sys

import cplex
from cplex.callbacks import UserCutCallback, LazyConstraintCallback
from cplex.exceptions import CplexError

from pyutilib.misc import import_file
from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.opt.base import SolverFactory
from pyomo.opt.parallel import SolverManagerFactory
from pyomo.opt.parallel.manager import solve_all_instances

# The class BendersLazyConsCallback
# allows to add Benders' cuts as lazy constraints.
#
class BendersLazyConsCallback(LazyConstraintCallback):

    def __call__(self):

		x = self.params.x
		sita = self.params.sita
		scen = self.params.numScen
		subLP = self.subLP
        # Get the current x solution
		solX = self.get_values(x)
		print ("lazy cut!!!solX="),
		print (solX)
		subLP.separate(solX)
		solMstSita = self.get_values(sita)
		#solMstSita = float("-Inf")
		
		if self.first == True:
			self.first = False
		for s in range(scen):
			if self.first == True or abs(solMstSita[s]-subLP.subObj[s])>1e-03:
				self.add(constraint=subLP.cutLhs[s],
									 sense="G",
									 rhs=subLP.cutRhs[s])

		xx=1
		"""
		"""
# The class BendersUserCutCallback
# allows to add Benders' cuts as user cuts.
#
class BendersUserCutCallback(UserCutCallback):

    def __call__(self):
		x = self.params.x
		sita = self.params.sita
		scen = self.params.numScen
		subLP = self.subLP
        # Get the current x solution
		solX = self.get_values(x)
		print ("user cut!!!solX="),
		print (solX)
		subLP.separate(solX)
		solMstSita = self.get_values(sita)
		#solMstSita = float("-Inf")
		if not self.is_after_cut_loop():
			return		
		if self.first == True:
			self.first = False
		for s in range(scen):
			if self.first == True or abs(solMstSita[s]-subLP.subObj[s])>1e-03:
				self.add(constraint=subLP.cutLhs[s],
									 sense="G",
									 rhs=subLP.cutRhs[s])
	
		xx = 1



class subproblemLP:
    #
    def __init__(self, params):
		self.numScen = params.numScen
		self.mstX = params.x
		self.mstSita = params.sita
		sb_mdl = import_file("sub.py").model
		sub_insts = []
		for s in range(self.numScen):
			sub_insts.append(
				sb_mdl.create_instance(name="sub"+str(s), \
								   filename="sub.dat"))
		solver_manager = SolverManagerFactory("serial")
		self.model = sb_mdl
		self.instance = sub_insts
		self.solver_manager = solver_manager
		

    # This method separates Benders' cuts violated by the current x solution.
    # Violated cuts are found by solving the worker LP
    #
    def separate(self, xSol):
		sub_insts = self.instance
		for instance in sub_insts:
			for i in instance.sI:
				instance.x[i] = max(xSol[i-1],0)
		solve_all_instances(self.solver_manager, 'cplex', sub_insts)
		subObj = []
		cutRhs = []
		cutLhs = []
		for s, inst in enumerate(sub_insts, 1):
			subObj.append(round(inst.oSub(),4))
			print("obj. value for scenario "+str(s)+ " = "+inst.name+" is "+\
					str(subObj[s-1])+"("+str(1.0/self.numScen*subObj[s-1])+")")
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
			cut = 1.0/self.numScen*cut		
			print ("added cut")
			print (cut)
			cutRhs.append(cut)
			#constraint i
			thecoefs = []
			for i in inst.sI:
				print ("+%f*x[%d]" %(1.0/self.numScen*inst.dual[inst.cSi[i]],i)),
				thecoefs.append(-1.0/self.numScen*inst.dual[inst.cSi[i]])

			theind = self.mstX[:]
			theind.append(self.mstSita[s-1])
			thecoefs.append(1.0)
			cutLhs.append(cplex.SparsePair(ind=theind, val=thecoefs))
			
		self.cutLhs = cutLhs
		self.cutRhs = cutRhs
		self.subObj = subObj
		
"""
"""


# This function creates the master ILP	
def createMasterILP(cpx,params):
	cpx.objective.set_sense(cpx.objective.sense.minimize)
	numComponents = 4
	numScen = 1
	params.numScen = numScen
	setupCost = 5
	cCR=[]
	cPR=[]
	kesi=[]
	for i in range(numComponents):
		cCR.append((i+1)*2)
		cPR.append(1)
		kesi.append(0)
		if i == 0:
			kesi[i] = 1
	#variables
	varNameX=[]
	x=[]
	for i in range(numComponents):
		varNameX.append("x" + str(i))
		x.append(cpx.variables.get_num())
		cpx.variables.add(obj=[cPR[i]],
						  lb=[0.0], ub=[1.0], types=["B"],
						  names=[varNameX[i]])
	params.x = x
	
	varNameSita=[]
	sita = []
	for i in range(numScen):
		varNameSita.append("sita" + str(i))
		sita.append(cpx.variables.get_num())
		cpx.variables.add(obj=[1.0],
						  lb=[-10000.0],
						  names=[varNameSita[i]])
	params.sita = sita
	
	varNameZ = "z"    
	cpx.variables.add(obj=[setupCost],
					  lb=[0.0],ub=[1.0], types=["B"],
					  names=[varNameZ])
    
	for i in range(numComponents):
		cpx.linear_constraints.add(
			lin_expr=[cplex.SparsePair([varNameX[i]],[1.0])],
			senses=["G"], 
			range_values=[0.0],
			rhs=[kesi[i]])
	for i in range(numComponents):
		cpx.linear_constraints.add(
			lin_expr=[cplex.SparsePair([varNameX[i],varNameZ],[-1.0,1.0,])],
			senses=["G"], 
			range_values=[0.0],
			rhs=[0])
class parameters:
	def __init__(self):
		self.numScen = []
		self.x = []
		self.sita = []
def benders_main():
	# Create master ILP
	cpx = cplex.Cplex()
	params = parameters()
	createMasterILP(cpx,params)
	subLP = subproblemLP(params)

	cpx.parameters.preprocessing.presolve.set(
		cpx.parameters.preprocessing.presolve.values.off)
	cpx.parameters.threads.set(1)

	cpx.parameters.mip.strategy.search.set(
		cpx.parameters.mip.strategy.search.values.traditional)
	#lazy constraints
	
	lazyBenders = cpx.register_callback(BendersLazyConsCallback)
	lazyBenders.params = params
	lazyBenders.subLP = subLP
	lazyBenders.first = True
	#user cut
	userBenders = cpx.register_callback(BendersUserCutCallback)
	userBenders.params = params
	userBenders.subLP = subLP
	userBenders.first = True
    
	# Solve the model
	cpx.solve()

	solution = cpx.solution
	print()
	print("Solution [x1,x2,x3,x4,sita,Z]: ", solution.get_values())
	print("Objective value: ", solution.get_objective_value())
	"""	
	if solution.get_status() == solution.status.MIP_optimal:
		# Write out the optimal tour
		succ = [-1] * numNodes
		for i in range(numNodes):
			sol = solution.get_values(x[i])
			for j in range(numNodes):
				if sol[j] > 1e-03:
					succ[i] = j
		print("Optimal tour:")
		i = 0
		while succ[i] != 0:
			print("%d, " % i, end=' ')
			i = succ[i]
		print(i)
	else:
		print("Solution status is not optimal")
        """

def usage():
	print("Usage:     don't use it!!")


if __name__ == "__main__":
	benders_main()