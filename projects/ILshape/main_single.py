#Author: Zhicheng Zhu
#Email: zzhu3@lamar.edu
#
#Description:
#main file for single-cut implementation
#1. send integer master problem into CPLEX
#2. add cuts in callback functions:
#	2.1 lazy constraints: when solution is integer feasible. Both benders cut and integer L-shaped cut
#	2.2 user cuts: when solution is not integer feasible. Benders cut only
#3. Best incumbent solution until no cuts can be found in 2.
#
#Last Update: 03/21/2018
#
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
import time
import genFile
# The class BendersLazyConsCallback
# allows to add Benders' cuts as lazy constraints.
#
lazy_cut = 0
user_cut = 0
class BendersLazyConsCallback(LazyConstraintCallback):

	def __call__(self):
		x = self.params.x;		
		sita = self.params.sita
		scen = self.params.numScen
		subLP = self.subLP
		subIP = self.subIP
		# Get the current x solution
		solX = self.get_values(x)
		#print ("lazy cut!!!solX="),
		#print (solX)

		solMstSita = self.get_values(sita)
		#solMstSita = float("-Inf")
		global lazy_cut
		for i in range(len(solX)):
			solX[i] = int(round(solX[i], 0))		
		if solX not in self.params.LazyX:
			self.params.LazyX.append(solX)
			subLP.separate(solX)
			subIP.separate(solX)
			lazy_cut += 1
			"""
			print("Lazy cut!!!!"),
			print ("x="),
			print(solX)
			print ("Lshape "),
			print ("LHS="),
			print (subIP.cutLhs[0])
			print ("RHS="),
			print (subIP.cutRhs[0])
			"""
			self.add(constraint=subIP.cutLhs,
					sense="G",
					rhs=subIP.cutRhs)

			if self.first == True or abs(solMstSita - 1.0/scen*subLP.subObj)>1e-02:
				lazy_cut += 1
				"""
				print ("Benders scen=%d",s),
				print ("LHS="),
				print (subLP.cutLhs[s])
				print ("RHS="),
				print (subLP.cutRhs[s])		
				"""
				self.add(constraint=subLP.cutLhs,
									 sense="G",
									 rhs=subLP.cutRhs)
			if self.first == True:
				self.first = False
			
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
		#print ("user cut!!!solX="),
		#print (solX)
		subLP.separate(solX)
		solMstSita = self.get_values(sita)
		#solMstSita = float("-Inf")
		if not self.is_after_cut_loop():
			return		
		for i in range(len(solX)):
			solX[i] = round(solX[i],2)
		if solX not in self.params.UserX:
			self.params.UserX.append(solX)			
			if self.first == True or abs(solMstSita-1.0/scen*subLP.subObj)>1e-02:
				global user_cut
				"""
				print("Benders cut!!!!"),
				print ("x="),
				print(solX)
				print ("LHS="),
				print (subLP.cutLhs[s])
				print ("RHS="),
				print (subLP.cutRhs[s])
				"""
				user_cut += 1
				self.add(cut=subLP.cutLhs,
									 sense="G",
									 rhs=subLP.cutRhs)
			if self.first == True:
				self.first = False	
		#xx = 1

#provide integer L-shaped cut to lazy constraint.
class subproblemIP:
    #
	def __init__(self, params):
		self.numScen = params.numScen
		self.mstX = params.x
		self.mstSita = params.sita
		sb_mdl = import_file("subIP.py").model
		sub_insts = []
		global directory
		for s in range(self.numScen):
			sub_path = directory + "\\data\\ScenNode" + str(s+1) + ".dat"
			sub_insts.append(
				sb_mdl.create_instance(name="subIP"+str(s+1), \
								   filename=sub_path))
		solver_manager = SolverManagerFactory("serial")
		self.model = sb_mdl
		self.instance = sub_insts
		self.solver_manager = solver_manager
		self.L = 0	#L, global minimum.

    # This method separates integer L-shaped cuts violated by the current x solution.
    # Violated cuts are found by solving the integer sub-problems
    #
	def separate(self, xSol):
		sub_insts = self.instance
		for instance in sub_insts:
			for i in instance.sI:
				instance.x[i] = max(xSol[i-1],0)
		
		solve_all_instances(self.solver_manager, 'cplex', sub_insts)
		subObj = 0
		cutRhs = 0
		cutLhs = 0
		q_x = 0
		s_x_e1 = int(sum(xSol[i] for i in range(len(xSol))))		#abs(s), integer L-shaped
		
		for s, inst in enumerate(sub_insts, 0):
			subObj += round(inst.oSub(),4)
			q_x += 1.0/self.numScen*inst.oSub()
		
		#get right-hand-side
		tmp = -(q_x - self.L)*(s_x_e1 - 1) + self.L
		cutRhs = tmp
		#get index and coefficient of x
		thecoefs = []
		for i in range(len(xSol)):
			tmp = q_x - self.L
			if xSol[i] == 1.0:
				thecoefs.append(-tmp)	#negative is because of LHS
			else:
				thecoefs.append(tmp)				
		theind = self.mstX[:]
		
		#get index and coefficient of sita
		theind.append(self.mstSita)		
		thecoefs.append(1.0)

		cutLhs = cplex.SparsePair(ind=theind, val=thecoefs)

		self.cutLhs = cutLhs
		self.cutRhs = cutRhs
		self.subObj = subObj

		
class subproblemLP:
    #
	def __init__(self, params):
		self.numScen = params.numScen
		self.mstX = params.x
		self.mstSita = params.sita
		sb_mdl = import_file("sub.py").model
		sub_insts = []
		global directory
		for s in range(self.numScen):
			sub_path = directory + "\\data\\ScenNode" + str(s+1) + ".dat"
			sub_insts.append(
				sb_mdl.create_instance(name="sub"+str(s+1), \
								   filename=sub_path))
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
		subObj = 0
		cutRhs = 0
		cutLhs = 0
		thecoefs = [0]*len(xSol)
		for s, inst in enumerate(sub_insts, 1):
			subObj += round(inst.oSub(),4)
			#print("obj. value for scenario "+str(s)+ " = "+inst.name+" is "+\
			#		str(subObj[s-1])+"("+str(1.0/self.numScen*subObj[s-1])+")")
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
			#constraint s1
			for i in inst.sI:
				for r in inst.sR_0:
					cut += 2*inst.dual[inst.cSs1[i,r]]
			#constraint u
			for t in inst.sT_0:
				cut += inst.dual[inst.cSu[t]]	
			#constraint v
			for i in inst.sI:
				for t in inst.sT_ex:
					for r in inst.sR:
						cut += inst.dual[inst.cSv[i,t,r]]
			cut = 1.0/self.numScen*cut		
			#print ("added cut")
			#print (cut)
			cutRhs += cut
			#constraint i
			for i in inst.sI:
				#print ("+%f*x[%d]" %(1.0/self.numScen*inst.dual[inst.cSi[i]],i)),
				thecoefs[i-1] = 1.0/self.numScen*inst.dual[inst.cSi[i]]

		theind = self.mstX[:]
		theind.append(self.mstSita)
		thecoefs.append(1.0)
		cutLhs = cplex.SparsePair(ind=theind, val=thecoefs)
			
		self.cutLhs = cutLhs
		self.cutRhs = cutRhs
		self.subObj = subObj
		
"""
"""


# This function creates the master ILP	
def createMasterILP(cpx,params):
	cpx.objective.set_sense(cpx.objective.sense.minimize)
	numComponents = params.numComp
	numScen = params.numScen
	setupCost = params.setupCost
	cPR = params.PR_cost
	cCR = params.CR_cost
	kesi = params.kesi
	"""
	cCR=[]
	cPR=[]
	kesi=[]
	for i in range(numComponents):
		cCR.append((i+1)*2)
		cPR.append(1)
		kesi.append(0)
		if i == 0:
			kesi[i] = 1
	"""
	#variables
	#for x
	varNameX=[]
	x=[]
	for i in range(numComponents):
		varNameX.append("x" + str(i))
		x.append(cpx.variables.get_num())
		cpx.variables.add(obj=[cPR[i]],
						  lb=[0.0], ub=[1.0], types=["B"],
						  names=[varNameX[i]])
	params.x = x
	#for sita
	varNameSita = "sita"
	sita = cpx.variables.get_num()
	cpx.variables.add(obj=[1.0],
					  lb=[-100000.0],
					  names=[varNameSita])
	params.sita = sita
	#for z
	varNameZ = "z"    
	cpx.variables.add(obj=[setupCost],
					  lb=[0.0],ub=[1.0], types=["B"],
					  names=[varNameZ])
    
	#constraints
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
	def __init__(self, I_in, w_in, d_in,PR_cost,CR_cost,kesi):
		self.numComp = I_in
		self.numScen = w_in
		self.setupCost = d_in
		self.PR_cost = PR_cost
		self.CR_cost = CR_cost
		self.kesi = kesi		
		self.x = []
		self.sita = []
		self.LazyX = []				#note history of x that pass to lazy constraint
		self.UserX = []				#note history of x that pass to user cuts
		
def benders_main(I_in,w_in,d_in,PR_cost,CR_cost,kesi):
	# Create master ILP
	cpx = cplex.Cplex()
	params = parameters(I_in,w_in,d_in,PR_cost,CR_cost,kesi)
	createMasterILP(cpx,params)
	subLP = subproblemLP(params)
	subIP = subproblemIP(params)

	start_time = time.clock()
	cpx.parameters.preprocessing.presolve.set(
		cpx.parameters.preprocessing.presolve.values.off)
	cpx.parameters.threads.set(1)

	cpx.parameters.mip.strategy.search.set(
		cpx.parameters.mip.strategy.search.values.traditional)
		
	#lazy constraints
	lazyBenders = cpx.register_callback(BendersLazyConsCallback)
	lazyBenders.params = params
	lazyBenders.subLP = subLP
	lazyBenders.subIP = subIP
	lazyBenders.first = True
	#user cut
	userBenders = cpx.register_callback(BendersUserCutCallback)
	userBenders.params = params
	userBenders.subLP = subLP
	userBenders.first = True
    
	# Solve the model
	cpx.solve()
	end_time = time.clock()
	solution = cpx.solution
	print()
	print("Solution: ", solution.get_values())
	print("Objective value: ", solution.get_objective_value())
	print ("time = "),
	print (end_time-start_time)
	print ("lazy cut"),
	print (lazy_cut)
	print ("user cut"),
	print (user_cut)
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

#############################################################	
comp_list = [3]
time_list = [7]#[10,20,30]
scen_list = [910]#[20,50,100]
d = 5 #setup cost
counter = 0
directory = "C:\\Users\\Yisha\\Desktop\\ZZC\\ILshape" 
for I in comp_list:
	for T in time_list:
		for w in scen_list:
			PR_cost = []
			CR_cost = []
			kesi = []
			counter += 1
			print (str(counter) + "_" + str(I) + "_" + str(T) + "_" + str(w))
			genFile.create_allFiles(I,T,w,d,directory,PR_cost,CR_cost,kesi)
			benders_main(I,w,d,PR_cost,CR_cost,kesi)
"""
if __name__ == "__main__":
	benders_main()
"""