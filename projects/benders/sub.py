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
model.pT		= Param()								#time horizon
model.pT_1		= Param()								#
model.sT		= RangeSet(0,model.pT)					#time set [0,pT]
model.sT_End	= RangeSet(0,model.pT-1)				#sT_end = [0,end-1]
model.sT_0		= RangeSet(1,model.pT)					#sT_0=[1,end]
model.pT_ex		= Param()								#Extended time horizon
model.sT_ex		= RangeSet(0,model.pT_ex)				#extended time set [0,pT_ex]
model.sT_exT	= RangeSet(model.pT+1, model.pT_ex)					#set[pT+1, pT_ex]
model.ps		= Param()								#starting time
model.pR 		= Param()								#max individuals
model.sR 		= RangeSet(1,model.pR)
model.sR_End	= RangeSet(1,model.pR-1)
model.sR_0		= RangeSet(2,model.pR)
model.pCPR		= Param(model.sI)						#PR cost
model.pCCR		= Param(model.sI)						#CR cost
model.pKesi		= Param(model.sI)						#failure state 
model.pd		= Param(default=5)						#set-up cost
model.w_shape	= Param(model.sI)
model.w_scale	= Param(model.sI)

model.sTest		= RangeSet(3,1)							#Test......

def pLT_init(model, i, r):
	return int(round((-math.log(random.uniform(0,1)))**(1.0/model.w_shape[i])*model.w_scale[i]))
model.pLT		= Param(model.sI, model.sR, initialize=pLT_init)

def scf_init(model):
	return ((i,r,t) for i in model.sI for r in model.sR_End for	t in model.sT\
			if t<=model.pT.value-model.pLT[i,r+1])
model.scf		= Set(dimen=3, initialize=scf_init)

def scg_init(model):
	return ( i for i in model.sI if model.pLT[i,1]<=model.pT.value)
model.scg		= Set(initialize=scg_init)

def scm_init(model):
	return ((i,r,t) for i in model.sI for r in model.sR_End for t in model.sT_exT\
			if t<=model.pT.value+model.pLT[i,r]-1)
model.scm		= Set(dimen=3, initialize=scm_init)

def scn_init(model):
	return ((i,r,t) for i in model.sI for r in model.sR_End for t in model.sT_exT\
			if t>=model.pT.value+model.pLT[i,r])
model.scn		= Set(dimen=3, initialize=scn_init)
#first stage variable
model.x			= Param(model.sI, within=NonNegativeReals, mutable=True)

######################################
#Variables
######################################
model.xs		= Var(model.sI, model.sT, model.sR, within=Binary)
model.ws		= Var(model.sI, model.sT_ex, model.sR, within=Binary)
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

def s_cf_rule(model, i, r, t):
	return model.xs[i,t,r]<=model.xs[i,t+model.pLT[i,r+1],r+1]
model.cSf = Constraint(model.scf, rule=s_cf_rule)

def s_cg_rule(model, i):## not finish
	return model.xs[i,model.pLT[i,1],1]==1
model.cSg = Constraint(model.scg, rule=s_cg_rule)

def s_ch_rule(model, i, r):
	return model.xs[i,0,r]==0
model.cSh = Constraint(model.sI, model.sR_0, rule=s_ch_rule)

def s_ci_rule(model, i):
	return model.x[i]==model.xs[i,0,1]
model.cSi = Constraint(model.sI, rule=s_ci_rule)

#def s_cj_rule(model, i, t):
#	return sum(model.xs[i,t,r]-model.xs[i,t-1,r] for r in model.sR)<=model.zs[t]
#model.cSj = Constraint(model.sI, model.sT_0, rule=s_cj_rule)

def s_ck_rule(model, i, t, r):
	return model.ws[i,t,r] == model.xs[i,t,r]-model.xs[i,t-1,r]
model.cSk = Constraint(model.sI, model.sT_0, model.sR, rule=s_ck_rule)

def s_cl_rule(model, i, r):
	return model.ws[i,0,r]==model.xs[i,0,r]
model.cSl = Constraint(model.sI, model.sR, rule=s_cl_rule)

def s_cm_rule(model, i, r, t):#not finish
	return model.ws[i,t,r]==model.ws[i,t-model.pT.value,r]
model.cSm = Constraint(model.scm, rule=s_cm_rule)

def s_cn_rule(model, i, r, t):#not finish
	return model.ws[i,t,r]==0
model.cSn = Constraint(model.scn, rule=s_cn_rule)

def s_co_rule(model, i, t, r):
	return model.u[i,t,r]-model.v[i,t,r]==\
			model.ws[i,t+model.pLT[i,r],r]-model.ws[i,t,r-1]
model.cSo = Constraint(model.sI, model.sT, model.sR_0, rule=s_co_rule)

def s_cp_rule(model, i, t):
	return model.u[i,t,1]-model.v[i,t,1]==\
			model.ws[i,t+model.pLT[i,1],1]-1
model.cSp = Constraint(model.sI, model.sT, rule=s_cp_rule)

def s_cq_rule(model, i, t,r):
	return model.u[i,t,r]+model.v[i,t,r]<=1
model.cSq = Constraint(model.sI, model.sT, model.sR, rule=s_cq_rule)
"""
"""
######################################
#Objective
######################################
def exp_stage2_rule(model):
	return sum(model.pCCR[i]*model.ws[i,model.pLT[i,1],1]+model.pCPR[i]*(1-model.ws[i,model.pLT[i,1],1])-\
				sum(model.pCPR[i1]*model.x[i1] for i1 in model.sI)-\
				sum((model.pCCR[i2]-model.pCPR[i2])*model.pKesi[i2] for i2 in model.sI)+\
				sum(model.pCCR[i]*(1-0.5*sum(model.u[i,t,r]+model.v[i,t,r] for t in model.sT))-\
					model.pCCR[i]*(1-0.5*(model.xs[i,model.pT,r]+model.xs[i,model.pT,r-1]))+\
					model.pCPR[i]*0.5*sum(model.u[i,t1,r]+model.v[i,t1,r] for t1 in model.sT)
					for r in model.sR_0
				   )
			  for i in model.sI
			  )+\
			sum(model.pd*model.zs[t] for t in model.sT_0)
model.oSub = Objective(rule=exp_stage2_rule,sense=minimize)
