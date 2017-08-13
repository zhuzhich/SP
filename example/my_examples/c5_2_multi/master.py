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
def prob_validator(model, value, s):
	return (value>=0) and (value<=1.0)
model.prob = Param(model.Scen, validate=prob_validator)
# first stage variable
#model.nX = Param(within=PositiveIntegers)
#model.sX = RangeSet(1,model.nX)
#second stage variable
model.nY = Param(within=PositiveIntegers)
model.sY = RangeSet(1,model.nY)
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
	
	
