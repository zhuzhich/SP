#  Implemention of extensive form. Solving the problem directly.
#  Author: Zhicheng Zhu
#  Email: zzhu3@lamar.edu
#  Last Update: 02/23/2018


#
##
##

from pyomo.core import *
import math
import random

model = AbstractModel(name="progressive hedging")

#
model.rc = Suffix(direction=Suffix.IMPORT)
model.lrc = Suffix(direction=Suffix.IMPORT)
model.urc = Suffix(direction=Suffix.IMPORT)
model.dual = Suffix(direction=Suffix.IMPORT)


#####################################
#Parameters
#####################################
model.NUMSCEN 	= Param()								#scenarios
model.Scen		= RangeSet(model.NUMSCEN)				
model.prob 		= Param(default=1.0/model.NUMSCEN)		#equal probability
model.pI		= Param()       						#component numbers
model.sI 		= RangeSet(model.pI)
model.pT		= Param()								#time horizon
model.pT_1		= Param()								#
model.sT		= RangeSet(0,model.pT)					#time set [0,pT]
model.sT_End	= RangeSet(0,model.pT-1)				#sT_end = [0,end-1]
model.sT_0		= RangeSet(1,model.pT)					#sT_0=[1,end]
model.pT_ex		= Param()								#Extended time horizon
model.sT_ex		= RangeSet(0,model.pT_ex)				#extended time set [0,pT_ex]
model.sT_exT	= RangeSet(model.pT+1, model.pT_ex)		#set[pT+1, pT_ex]
model.ps		= Param()								#starting time
model.pR 		= Param()								#max individuals
model.sR 		= RangeSet(1,model.pR)
model.sR_End	= RangeSet(1,model.pR-1)
model.sR_0		= RangeSet(2,model.pR)
model.pd		= Param(default=5)						#set-up cost

model.pCPR		= Param(model.sI)  #PR cost	
model.pCCR		= Param(model.sI)	 #CR cost
model.pKesi		= Param(model.sI)  #failure state
model.w_shape	= Param(model.sI) #weibull shape
model.w_scale	= Param(model.sI) #weibull scale
model.pLT		= Param(model.sI, model.sR)

#For the large number of components
#generate parameters.
"""
#PR cost
def pCPR_init(model, i):
	return 1
model.pCPR		= Param(model.sI, initialize=pCPR_init)  #PR cost	

#CR cost
def pCCR_init(model, i):
	return i*2					
model.pCCR		= Param(model.sI, initialize=pCCR_init)	 #CR cost

#kesi
def pKesi_init(model, i):
	if i==1:
		return 1
	else:
		return 0	
model.pKesi		= Param(model.sI, initialize=pKesi_init)  #failure state

#weibull shape parameter
def w_shape_init(model, i):
	return 6
model.w_shape		= Param(model.sI, initialize=w_shape_init) #weibull shape

#weibull scale parameter
def w_scale_init(model, i):
	return i
model.w_scale		= Param(model.sI, initialize=w_scale_init) #weibull scale


#random life time
def pLT_init(model, i, r):
	#random.seed(i**2+r**2)
	if r == 1:
		if model.pKesi[i]==1:
			return 0
	#return i #for test only..
		else:
			LT=int(round((-math.log(random.uniform(0,1)))**\
					(1.0/model.w_shape[i])*model.w_scale[i]))-\
					model.ps
			return max(1,LT)
	else:	
		LT=int(round((-math.log(random.uniform(0,1)))**\
					(1.0/model.w_shape[i])*model.w_scale[i]))
		return max(1,LT)
model.pLT		= Param(model.sI, model.sR, initialize=pLT_init)
if 0:
	for i in model.sI:
		for r in model.sR:
			print ("(%d,%d)=%d" %(i,r,model.pLT[i,r]()))
"""
#set of constraint f
def scf_init(model):
	return ((i,r,t) for i in model.sI for r in model.sR_End for t in model.sT\
			if t<=model.pT-model.pLT[i,r+1])
model.scf		= Set(dimen=3, initialize=scf_init)

#set of constraint g
def scg_init(model):
	return ( i for i in model.sI if model.pLT[i,1]<=model.pT)
model.scg		= Set(dimen=1, initialize=scg_init)

#set of constraint m
def scm_init(model):
	return ((i,r,t) for i in model.sI for r in model.sR for t in model.sT_exT\
			  if t<=model.pT+model.pLT[i,r]-1)
model.scm		= Set(dimen=3, initialize=scm_init)

#set of constraint n
def scn_init(model):
	return ((i,r,t) for i in model.sI for r in model.sR for t in model.sT_exT\
			  if t>=model.pT+model.pLT[i,r])
model.scn		= Set(dimen=3, initialize=scn_init)

######################################
#Variables
######################################
if 1:
#integer
	model.x 		= Var(model.sI, within=Binary)
	model.z 		= Var(within=Binary)	
	model.xs		= Var(model.sI, model.sT, model.sR, within=Binary)
	model.ws		= Var(model.sI, model.sT_ex, model.sR, within=Binary)
	model.zs		= Var(model.sT, within=Binary)
	model.u			= Var(model.sI, model.sT, model.sR, within=Binary)
	model.v			= Var(model.sI, model.sT, model.sR, within=Binary)
else:
#relax
	model.x 		= Var(model.sI, within=NonNegativeReals)
	model.z 		= Var(within=NonNegativeReals)
	model.xs		= Var(model.sI, model.sT, model.sR, within=NonNegativeReals)
	model.ws		= Var(model.sI, model.sT_ex, model.sR, within=NonNegativeReals)
	model.zs		= Var(model.sT, within=NonNegativeReals)
	model.u			= Var(model.sI, model.sT, model.sR, within=NonNegativeReals)
	model.v			= Var(model.sI, model.sT, model.sR, within=NonNegativeReals)
#model.test = Var(bounds=(0,0.5),initialize=1)
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
	return sum(model.ws[i,t,r] for r in model.sR)<=model.zs[t]
model.cSd = Constraint(model.sI, model.sT_0, rule=s_cd_rule)

def s_ce_rule(model, i):
	return model.xs[i,0,1]<=model.zs[0]
model.cSe = Constraint(model.sI, rule=s_ce_rule)

def s_ce1_rule(model):
	return model.z == model.zs[0]
model.cSe1 = Constraint(rule=s_ce1_rule)

def s_cf_rule(model, i, r, t):
	return model.xs[i,t,r]<=model.xs[i,t+model.pLT[i,r+1],r+1]
model.cSf = Constraint(model.scf, rule=s_cf_rule)

def s_cg_rule(model, i):
	return model.xs[i,model.pLT[i,1],1]==1
model.cSg = Constraint(model.scg, rule=s_cg_rule)

def s_ch_rule(model, i, r):
	return model.xs[i,0,r]==0
model.cSh = Constraint(model.sI, model.sR_0,  rule=s_ch_rule)

def s_ci_rule(model, i):
	return model.xs[i,0,1]==model.x[i]
model.cSi = Constraint(model.sI,  rule=s_ci_rule)

def s_cj_rule(model, i):
	return model.x[i] >= model.pKesi[i]
model.cSj = Constraint(model.sI, rule=s_cj_rule)

def s_ck_rule(model, i, t, r):
	return model.ws[i,t,r] == model.xs[i,t,r]-model.xs[i,t-1,r]
model.cSk = Constraint(model.sI, model.sT_0, model.sR,  rule=s_ck_rule)

def s_cl_rule(model, i, r):
	return model.ws[i,0,r]==model.xs[i,0,r]
model.cSl = Constraint(model.sI, model.sR,  rule=s_cl_rule)

def s_cm_rule(model, i, r, t):
	return model.ws[i,t,r]==model.ws[i,t-model.pT.value,r]
model.cSm = Constraint(model.scm, rule=s_cm_rule)

def s_cn_rule(model, i, r, t):
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

def s_cq_rule(model, i, t, r):
	return model.u[i,t,r]+model.v[i,t,r]<=1
model.cSq = Constraint(model.sI, model.sT, model.sR, rule=s_cq_rule)

def s_cr_rule(model, i, t, r):
	return model.u[i,t,r]<=1 #and model.u[i,t,r]>=0
model.cSr = Constraint(model.sI, model.sT, model.sR, rule=s_cr_rule)

def s_cs_rule(model, i, t, r):
	return model.v[i,t,r]<=1 #and model.v[i,t,r]>=0
model.cSs = Constraint(model.sI, model.sT, model.sR, rule=s_cs_rule)

def s_cs1_rule(model, i, r):
	return sum(model.u[i,t,r]+model.v[i,t,r] for t in model.sT) <= 2
model.cSs1 = Constraint(model.sI, model.sR_0, rule=s_cs1_rule)

def s_ct_rule(model, i, t, r):
	return model.xs[i,t,r]<=1 #and model.xs[i,t,r]>=0
model.cSt = Constraint(model.sI, model.sT, model.sR, rule=s_ct_rule)

def s_cu_rule(model, t):
	return model.zs[t]<=1 #and model.zs[t]>=0
model.cSu = Constraint(model.sT, rule=s_cu_rule)

def s_cv_rule(model, i, t, r):
	return model.ws[i,t,r]<=1 #and model.ws[i,t,r]>=0
model.cSv = Constraint(model.sI, model.sT_ex, model.sR, rule=s_cv_rule)

def s_cw_rule(model, i):
	return model.x[i]<=1 
model.cSw = Constraint(model.sI, rule=s_cw_rule)
"""
"""
######################################
#Objective
######################################
def computeFirstStageObj_rule(model):
		return model.pd*model.z+\
			   sum(model.pCPR[i]*model.x[i] for i in model.sI)+\
			   sum((model.pCCR[i]-model.pCPR[i])*model.pKesi[i] for i in model.sI)
model.firstStageObj = Expression(rule=computeFirstStageObj_rule)	
	
def computeSecondStageObj_rule(model):
	model.idr1   = sum(model.pCCR[i]*model.ws[i,model.pLT[i,1],1]+\
					model.pCPR[i]*(1-model.ws[i,model.pLT[i,1],1])\
					for i in model.sI)
	model.stage1 = sum(model.pCPR[i1]*model.x[i1] for i1 in model.sI)+\
				sum((model.pCCR[i2]-model.pCPR[i2])*model.pKesi[i2] for i2 in model.sI)
	model.othIdr = sum(sum(model.pCCR[i]*(1-0.5*sum(model.u[i,t,r]+model.v[i,t,r] for t in model.sT))-\
					model.pCCR[i]*(1-0.5*(model.xs[i,model.pT,r]+model.xs[i,model.pT,r-1]))+\
					model.pCPR[i]*0.5*sum(model.u[i,t1,r]+model.v[i,t1,r] for t1 in model.sT)\
					for r in model.sR_0)-0.5*model.pCPR[i]\
					for i in model.sI)
	model.setup  = sum(model.pd*model.zs[t] for t in model.sT_0)
	return model.idr1 - model.stage1 + model.othIdr + model.setup
model.secondStageObj = Expression(rule=computeSecondStageObj_rule)

def totalObj_rule(model):
	return model.firstStageObj + model.secondStageObj
model.totalObj = Objective(rule=totalObj_rule,sense=minimize)
