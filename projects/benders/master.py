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
model.NUMSCEN 	= Param()								#scenarios
model.Scen		= RangeSet(model.NUMSCEN)				
model.prob 		= Param(default=1.0/model.NUMSCEN)		#equal probability
model.pI		= Param()       						#component numbers
model.sI 		= RangeSet(model.pI)
#model.pT		= Param()								#time horizon
#model.sT		= RangeSet(0,model.pT)					#time set [0,pT]
#model.ps		= Param(default=2)						#starting time
#model.pR 		= Param()								#max individuals
#model.sR 		= RangeSet(model.pR)
model.pCPR		= Param(model.sI)						#PR cost
model.pCCR		= Param(model.sI)						#CR cost
model.pKesi		= Param(model.sI)						#failure state 
model.pd		= Param()								#set-up cost
model.CUTS = Set(within=PositiveIntegers, ordered=True)
#Subproblem constraint number, for dual solution.
#dual solution
model.dSolution = Param(model.Scen, model.CUTS, \
					default=0.0, mutable=True)
model.pdH 		= Param(model.Scen)
model.pdT 		= Param()
######################################
#Variables
######################################
model.x 		= Var(model.sI,within=Binary)
model.z 		= Var(within=Binary)
model.sita 		= Var(model.Scen)
######################################
#Constraints
######################################
def m_c1_rule(model,i):
	return model.x[i] >= model.pKesi[i]
model.c_x_1 	= Constraint(model.sI, rule=m_c1_rule)	
def m_c2_rule(model,i):
	return model.z >= model.x[i]
model.c_x_2 	= Constraint(model.sI, rule=m_c2_rule)

model.Cut_Defn 	= ConstraintList()

######################################
#Objective
######################################
def min_obj_rule(model):
	return model.pd*model.z+\
		   sum(model.pCPR[i]*model.x[i] for i in model.sI)+\
		   sum((model.pCCR[i]-model.pCPR[i])*model.pKesi[i] for i in model.sI)+\
		   sum(model.sita[i] for i in model.Scen)
model.oMaster 	= Objective(rule=min_obj_rule,sense=minimize)
	
	
