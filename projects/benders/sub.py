#  Implemention of bender's decomposition of SORP considering PM
#  Author: Zhicheng Zhu
#  Email: zzhu3@lamar.edu


#Subproblem
##
##

from pyomo.core import *
import math
import random

model = AbstractModel(name="subProblem")

#
model.rc = Suffix(direction=Suffix.IMPORT)
model.lrc = Suffix(direction=Suffix.IMPORT)
model.urc = Suffix(direction=Suffix.IMPORT)
model.dual = Suffix(direction=Suffix.IMPORT)


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
model.sT_End	= model.sT[:-1]							#sT-end = [0,end-1]
model.sT_0		= model.sT[1:]							#sT-0=[1,end]
model.ps		= Param(default=2)						#starting time
model.pR 		= Param(default=model.pT.value + 1)		#max individuals
model.sR 		= range(1,model.pR.value)
model.sR_End	= model.sR[:-1]
model.sR_0		= model.sR[1:]	
model.pCPR		= Param(model.sI)						#PR cost
model.pCCR		= Param(model.sI)						#CR cost
model.pKesi		= Param(model.sI)						#failure state 
model.pd		= Param(default=5)						#set-up cost
model.w_shape	= Param(model.sI)
model.w_scale	= Param(model.sI)

def pLT_init(model, i, r):
	return round((-math.log(random.uniform(0,1)))**(1.0/model.w_shape[i])*model.w_scale[i])
model.pLT		= Param(model.sI, model.sR, initialize=pLT_init)

#first stage variable
model.x			= Param(model.sI, within=NonNegativeReals, mutable=True)

######################################
#Variables
######################################
model.xs		= Var(model.sI, model.sT, model.sR, within=Binary)
model.ws		= Var(model.sI, model.sT, model.sR, within=Binary)
model.zs		= Var(model.sT_0, within=Binary)
model.u			= Var(model.sI, model.sT, model.sR, within=Binary)
model.v			= Var(model.sI, model.sT, model.sR, within=Binary)

######################################
#Constraints
######################################

def s_cb_rule(model, i, t, r):
	return model.xs[i,t,r]<=model.xs[i,t+1,r]
model.cSb = Constraint(model.sI, model.sT_End, model.sR, rule=s_cb_rule)

def s_cc_rule(model, i, t, r):
	return model.xs[i,t+1,r+1]<=model.xs[i,t,r]
model.cSc = Constraint(model.sI, model.sT_End, model.sR_End, rule=s_cc_rule)

def s_cd_rule(model, i, t):
	return sum(model.xs[i,t,r]-model.xs[i,t-1,r] for r in model.sR)<=model.zs[t]
model.cSd = Constraint(model.sI, model.sT_0, rule=s_cd_rule)

#
#skip constraint e
#
def s_cf_rule(model, i, t, r):
	print (i,t,r)
	if t<=model.pT.value-int(model.pLT[i,r+1]):
		#return model.xs[i,t,r]<=model.xs[i,t+int(model.pLT[i,r+1]),r+1]
		return True
	else:
		return True
model.cSf = Constraint(model.sI,  model.sT, model.sR_End, rule=s_cf_rule)
"""
####
def s_cg_rule(model, i, t):## not finish
	return model.xs[i,t+model.pLT[i,1],1]==1
model.cS2 = Constraint(model.sI, model.sT_End, model.sR_End, rule=s_cg_rule)
###
def s_ch_rule(model, i, r):
	return model.xs[i,0,r]=0
model.cSh = Constraint(model.sI, model.sR_0, rule=s_ch_rule)

def s_ci_rule(model, i):
	return mode.x[i]=model.xs[i,0,1]
model.cid = Constraint(model.sI, rule=s_ci_rule)

#def s_cj_rule(model, i, t):
#	return sum(model.xs[i,t,r]-model.xs[i,t-1,r] for r in model.sR)<=model.zs[t]
#model.cSj = Constraint(model.sI, model.sT_0, rule=s_cj_rule)

def s_ck_rule(model, i, t, r):
	return model.ws[i,t,r] = model.xs[i,t,r]-model.xs[i,t-1,r]
model.cSk = Constraint(model.sI, model.sT_0, model.sR, rule=s_ck_rule)

def s_cl_rule(model, i, r):
	return model.ws[i,0,r]=model.xs[i,0,r]
model.cSl = Constraint(model.sI, model.sR, rule=s_cl_rule)

def s_cm_rule(model, i, r, t):#not finish
	return smodel.ws[i,t,r]=model.ws[i,t-model.pT.value,r]
model.cSm = Constraint(model.sI, model.sR, model.sT, rule=s_cm_rule)

def s_cn_rule(model, i, r, t):#not finish
	return model.ws[i,t,r]=0
model.cSn = Constraint(model.sI, model.sR, model.sT_0, rule=s_cn_rule)

def s_co_rule(model, i, t):#not finish
	return 
model.cSo = Constraint(model.sI, model.sT_0, rule=s_co_rule)

def s_cp_rule(model, i, t):#not finish
	return sum(model.xs[i,t,r]-model.xs[i,t-1,r] for r in model.sR)<=model.zs[t]
model.cSp = Constraint(model.sI, model.sT_0, rule=s_cp_rule)

def s_cq_rule(model, i, t):
	return model.u[i,t,r]+model.v[i,t,r]<=1
model.cSq = Constraint(model.sI, model.sT, model.sR, rule=s_cq_rule)
"""

######################################
#Objective
######################################
def exp_stage2_rule(model):
	return 1
model.oSub = Objective(rule=exp_stage2_rule,sense=minimize)
