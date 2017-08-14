#  Implemention of bender's decomposition of SORP considering PM
#  Author: Zhicheng Zhu
#  Email: zzhu3@lamar.edu

# Master Problem
##
##

from pyomo.core import *
model = AbstractModel(name="Master")

#####################################
#Parameters
#####################################
model.NUMSCEN 	= Param(default=4)						#scenarios
model.Scen		= RangeSet(model.NUMSCEN)				
model.prob 		= Param(default=1.0/model.NUMSCEN.value)#equal probability
model.pI		= Param(default=4)       				#component numbers
model.sI 		= RangeSet(model.pI)
model.pT		= Param(default=5)						#time horizon
model.sT		= range(0,model.pT.value+1)				#time set [0,pT]
model.ps		= Param(default=2)						#starting time
model.pR 		= Param(default=model.pTime.value + 1)	#max individuals
model.sR 		= RangeSet(model.pR)
model.pCPR		= Param(model.sI)						#PR cost
model.pCCR		= Param(model.sI)						#CR cost
model.pKesi		= Param(model.sI)						#failure state 

#second stage variable
model.nY = Param(within=PositiveIntegers)
model.sY = RangeSet(model.nY)
model.CUTS = Set(within=PositiveIntegers, ordered=True)
#Subproblem constraint number, for dual solution.
#dual solution
model.dSolution = Param(model.Scen, model.CUTS, \
			default=0.0, mutable=True)
model.pH = Param(model.Scen)
model.pT = Param()
######################################
#Variables
######################################
model.x =  Var()
model.Min_Stage2_cost = Var(model.Scen)
######################################
#Constraints
######################################
def m_c1_rule(model):
	return model.x<=10 and model.x>=0
model.c_x_1 = Constraint(rule=m_c1_rule)	

model.Cut_Defn = ConstraintList()

######################################
#Objective
######################################
def min_obj_rule(model):
	return sum(model.Min_Stage2_cost[i] for i in model.Scen)
model.oMaster = Objective(rule=min_obj_rule,sense=minimize)
	
	
