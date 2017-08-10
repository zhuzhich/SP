#  Implement Bender's decomposition in Birge book
#  Chapter 5, example 1, pp184
#  Author: Zhicheng Zhu
#  Email: zzhu3@lamar.edu

#Subproblem
##
##

from pyomo.core import *

model = AbstractModel(name="subProblem")

#
model.rc = Suffix(direction=Suffix.IMPORT)
model.lrc = Suffix(direction=Suffix.IMPORT)
model.urc = Suffix(direction=Suffix.IMPORT)
model.dual = Suffix(direction=Suffix.IMPORT)


#####################################
#Parameters
#####################################
model.NUMSCEN = Param(within=PositiveIntegers)
model.SCEN = RangeSet(model.NUMSCEN)
def prob_validator(model, value, s):
	return (value>=0) and (value<=1.0)
model.prob = Param(model.SCEN, validate=prob_validator)

#kinds of materials
model.NMAT = Param(within=PositiveIntegers)
model.MAT = RangeSet(model.NMAT)

#kinds of products
model.pNPRO = Param(within=PositiveIntegers)
model.sPRO = RangeSet(model.pNPRO)

#cost of each materials
model.pMatUse = Param(model.MAT, model.sPRO, within=NonNegativeReals)

#
model.pMatUseX = Param(model.MAT, within=NonNegativeReals)

#material from first stage x1 x2
model.xMat = Param(model.MAT, within=NonNegativeReals, mutable=True)

#maximun number of products d1 d2
model.pMaxPro = Param(model.sPRO, within=PositiveIntegers)

#cost per product
model.pCostPerPro = Param(model.sPRO)

######################################
#Variables
######################################
model.yPro =  Var(model.sPRO, within=NonNegativeReals)

######################################
#Constraints
######################################
def mat_pro_rule(model, i):
	return sum([model.pMatUse[i,j]*model.yPro[j] for j in model.sPRO])<=\
			model.pMatUseX[i]*model.xMat[i]
model.cMatPro = Constraint(model.MAT, rule=mat_pro_rule)

def max_pro_rule(model, i):
	return model.yPro[i] <= model.pMaxPro[i]
model.cMaxPro = Constraint(model.sPRO, rule=max_pro_rule)

######################################
#Objective
######################################
def exp_stage2_rule(model):
	return sum(model.pCostPerPro[i]*model.yPro[i] for i in model.sPRO)
model.oSub = Objective(rule=exp_stage2_rule,sense=minimize)
