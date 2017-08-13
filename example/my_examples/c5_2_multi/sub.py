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
# first stage variable
#model.nX = Param(within=PositiveIntegers)
#model.sX = RangeSet(1,model.nX)
#second stage variable
model.nY = Param(within=PositiveIntegers)
model.sY = RangeSet(1,model.nY)

#kesi
model.pKesi = Param(within=PositiveIntegers)

#first stage
model.x = Param(within=NonNegativeReals, mutable=True)

######################################
#Variables
######################################
model.y =  Var(model.sY, within=NonNegativeReals)

######################################
#Constraints
######################################
def s_c1_rule(model):
	return model.y[model.sY[1]]-model.y[model.sY[2]] == model.pKesi-model.x
model.cY1 = Constraint(rule=s_c1_rule)

######################################
#Objective
######################################
def exp_stage2_rule(model):
	return sum(model.y[i] for i in model.sY)
model.oSub = Objective(rule=exp_stage2_rule,sense=minimize)
