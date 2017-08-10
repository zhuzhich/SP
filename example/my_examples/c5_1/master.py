#  Implement Bender's decomposition in Birge book
#  Chapter 5, example 1, pp184
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
model.NUMSCEN = Param(within=PositiveIntegers)
model.SCEN = RangeSet(model.NUMSCEN)
def prob_validator(model, value, s):
	return (value>=0) and (value<=1.0)
model.prob = Param(model.SCEN, validate=prob_validator)
# kinds of materials
model.NMAT = Param(within=PositiveIntegers)
model.MAT = RangeSet(model.NMAT)
#kinds of products
model.pNPRO = Param(within=PositiveIntegers)
model.sPRO = RangeSet(model.pNPRO)

# number of total materials
model.TotalMat = Param(within=PositiveIntegers)

# minimun input
model.MinInput = Param(model.MAT, within=PositiveIntegers)

# material cost
model.pMatCost = Param(model.MAT, within=PositiveIntegers)

model.CUTS = Set(within=PositiveIntegers, ordered=True)
model.dMatPro = Param(model.MAT, model.SCEN, model.CUTS, \
				default=0.0, mutable=True)
model.dMaxPro = Param(model.sPRO, model.SCEN, model.CUTS, \
				default=0.0, mutable=True)
model.pHMatPro = Param(model.MAT, model.SCEN, \
					within=NonNegativeReals)
model.pHMaxPro = Param(model.sPRO, model.SCEN, \
					within=NonNegativeReals)
model.pTMatPro = Param(model.MAT, model.MAT)
model.pTMaxPro = Param(model.sPRO, model.MAT)
######################################
#Variables
######################################
model.xMat =  Var(model.MAT, within=NonNegativeReals)
model.Min_Stage2_cost = Var()
######################################
#Constraints
######################################
def max_material_rule(model):
	return sum([model.xMat[i] for i in model.MAT])<=model.TotalMat
model.cMaxMat = Constraint(rule=max_material_rule)	

def min_material_rule(model, i):
	return model.xMat[i] >= model.MinInput[i]
model.cMinMat = Constraint(model.MAT, rule=min_material_rule)

model.Cut_Defn = ConstraintList()

######################################
#Objective
######################################
def min_cost_rule(model):
	return sum([model.pMatCost[i]*model.xMat[i] for i in model.MAT])+\
			model.Min_Stage2_cost

model.oMaster = Objective(rule=min_cost_rule,sense=minimize)
	
	
