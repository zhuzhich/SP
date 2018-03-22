#  Author: Zhicheng Zhu
#  Email: zzhu3@lamar.edu

# Master Problem
##
#Last Update: 03/21/2018
#

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
model.pd		= Param()								#set-up cost
#model.CUTS = Set(within=PositiveIntegers, ordered=True)

#For the large number of components
#generate parameters.
def pCPR_init(model, i):
	return 1
model.pCPR		= Param(model.sI, initialize=pCPR_init)  #PR cost	
def pCCR_init(model, i):
	return i*2					
model.pCCR		= Param(model.sI, initialize=pCCR_init)	 #CR cost
def pKesi_init(model, i):
	if i==1:
		return 1
	else:
		return 0	
model.pKesi		= Param(model.sI, initialize=pKesi_init) #failure state

######################################
#Variables
######################################
model.x 		= Var(model.sI, within=Binary)
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

#dynamic constraint list for benders cuts
model.Cut_Defn 	= ConstraintList()

######################################
#Objective
######################################
def min_obj_rule(model):
	if 0:
		return model.pd*model.z+\
			   sum(model.pCPR[i]*model.x[i] for i in model.sI)+\
			   sum((model.pCCR[i]-model.pCPR[i])*model.pKesi[i] for i in model.sI)+\
			   sum(model.sita[i] for i in model.Scen)
	else:#no constant
		return model.pd*model.z+\
			   sum(model.pCPR[i]*model.x[i] for i in model.sI)+\
			   sum(model.sita[i] for i in model.Scen)		
model.oMaster 	= Objective(rule=min_obj_rule,sense=minimize)
	
	
