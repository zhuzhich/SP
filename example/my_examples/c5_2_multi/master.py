#  Implement Bender's decomposition multi cut in Birge book
#  Chapter 5, example 2, pp199 
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
model.NUMSCEN = Param()
model.Scen = RangeSet(model.NUMSCEN)
# first stage variable
model.nX = Param(within=PositiveIntegers)
model.sX = RangeSet(model.nX)
#second stage variable
model.nY = Param(within=PositiveIntegers)
model.sY = RangeSet(model.nY)
model.CUTS = Set(within=PositiveIntegers, ordered=True)

######################################
#Variables
######################################
model.x =  Var(model.sX)
model.Min_Stage2_cost = Var(model.Scen)
######################################
#Constraints
######################################
def m_c1_rule(model,i):
	return model.x[i]<=10 and model.x[i]>=0
model.c_x_1 = Constraint(model.sX,rule=m_c1_rule)	

model.Cut_Defn = ConstraintList()

######################################
#Objective
######################################
def min_obj_rule(model):
	return sum(model.Min_Stage2_cost[i] for i in model.Scen)

model.oMaster = Objective(rule=min_obj_rule,sense=minimize)
	
	
