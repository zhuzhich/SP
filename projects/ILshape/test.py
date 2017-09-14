#  Implemention of bender's decomposition of SORP considering PM
#  Author: Zhicheng Zhu
#  Email: zzhu3@lamar.edu

# Master Problem
##
##

from pyomo.core import *
model = AbstractModel(name="Master")

######################################
#Variables
######################################
model.x 		= Var(within=Integers, bounds=(0,1))
model.y			= Var(within=Integers, bounds=(0,1))

######################################
#Objective
######################################
def min_obj_rule(model):
	return (model.y-1)**2
model.oMaster 	= Objective(rule=min_obj_rule,sense=minimize)
	
	
