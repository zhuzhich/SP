#  Implemention of bender's decomposition of SORP considering PM
#  Author: Zhicheng Zhu
#  Email: zzhu3@lamar.edu

# Master Problem
##
##

from pyomo.core import *
def pyomo_create_model(options=None, model_options=None):

	model = ConcreteModel(name="Master_Int")

	#####################################
	#Parameters
	#####################################
	model.NUMSCEN 	= Param(initialize=1.0)								#scenarios
	model.Scen		= RangeSet(model.NUMSCEN)				
	model.prob 		= Param(initialize=1.0/model.NUMSCEN)		#equal probability
	model.pI		= Param(initialize=4.0)       						#component numbers
	model.sI 		= RangeSet(model.pI)
	model.pd		= Param(initialize=5.0)								#set-up cost
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
	
	return model
	
	
