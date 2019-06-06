#  Author: Zhicheng Zhu
#  Email: zzhu3@lamar.edu

#  Implemention of extensive form. Whole problem.
#
#Last Update: 03/21/2018
#

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
model.iter  	= Param()

#For the large number of components
#generate parameters.

#PR cost
def pCPR_init(model, i):
	return 1
model.pCPR		= Param(model.sI, initialize=pCPR_init)  #PR cost	

#CR cost
def pCCR_init(model, i):
	random.seed((i-1)*30)
	temp = random.uniform(6,16)
	return round(temp,1)				
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
	random.seed((i-1)*20) ###control the seed     
	temp = random.uniform(4,7)
	return round(temp,1)
model.w_shape		= Param(model.sI, initialize=w_shape_init) #weibull shape

#weibull scale parameter
def w_scale_init(model, i):
	random.seed((i-1)*10) ###control the seed     
	temp = random.uniform(1,8)#(4,11)
	return round(temp,1)
model.w_scale		= Param(model.sI, initialize=w_scale_init) #weibull scale

#random life time
def pLT_init(model, i, r, w):
	random.seed(i-1+r-1+w-1+model.iter)
	if r == 1:
		if model.pKesi[i]==1:
			return 0
	#return i #for test only..
		else:
			ran_num = random.uniform(0,1)
			part1 = math.log(ran_num)
			s_inv = 1.0/model.w_shape[i]
			surv_time = model.ps	
			part2 = (surv_time/model.w_scale[i])**model.w_shape[i];
			part3 = part2 - part1;
			LT1 = round((part3**s_inv)*model.w_scale[i]) - surv_time;
			'''
			ran_num_log = -math.log(ran_num)
			s_inv = 1.0/model.w_shape[i]	
			LT1 = int((ran_num_log**s_inv)*model.w_scale[i]) - model.ps	
			'''
			LT = int(max(1,LT1))
			return LT
	else:	
		ran_num = random.uniform(0,1)
		ran_num_log = -math.log(ran_num)
		s_inv = 1.0/model.w_shape[i]	
		LT1 = round((ran_num_log**s_inv)*model.w_scale[i])				
		LT = int(max(1,LT1))
		return LT
model.pLT		= Param(model.sI, model.sR, model.Scen, initialize=pLT_init)

#set of constraint f
def scf_init(model):
	return ((i,r,t,w) for i in model.sI for r in model.sR_End for t in model.sT\
			for w in model.Scen if t<=model.pT-model.pLT[i,r+1,w])
model.scf		= Set(dimen=4, initialize=scf_init)

#set of constraint g
def scg_init(model):
	return ( (i,w) for i in model.sI for w in model.Scen if model.pLT[i,1,w]<=model.pT)
model.scg		= Set(dimen=2, initialize=scg_init)

#set of constraint m
def scm_init(model):
	return ((i,r,t,w) for i in model.sI for r in model.sR for t in model.sT_exT\
			for w in model.Scen if t<=model.pT+model.pLT[i,r,w]-1)
model.scm		= Set(dimen=4, initialize=scm_init)

#set of constraint n
def scn_init(model):
	return ((i,r,t,w) for i in model.sI for r in model.sR for t in model.sT_exT\
			for w in model.Scen if t>=model.pT+model.pLT[i,r,w])
model.scn		= Set(dimen=4, initialize=scn_init)

######################################
#Variables
######################################
if 1:
#integer
	model.x 		= Var(model.sI, within=Binary)
	model.xs		= Var(model.sI, model.sT, model.sR, model.Scen, within=Binary)
	model.ws		= Var(model.sI, model.sT_ex, model.sR, model.Scen, within=Binary)
	model.zs		= Var(model.sT, model.Scen, within=Binary)
	model.u			= Var(model.sI, model.sT, model.sR, model.Scen, within=Binary)
	model.v			= Var(model.sI, model.sT, model.sR, model.Scen, within=Binary)
else:
#relax
	model.x 		= Var(model.sI, within=NonNegativeReals)
	model.xs		= Var(model.sI, model.sT, model.sR, model.Scen, within=NonNegativeReals)
	model.ws		= Var(model.sI, model.sT_ex, model.sR, model.Scen, within=NonNegativeReals)
	model.zs		= Var(model.sT, model.Scen, within=NonNegativeReals)
	model.u			= Var(model.sI, model.sT, model.sR, model.Scen, within=NonNegativeReals)
	model.v			= Var(model.sI, model.sT, model.sR, model.Scen, within=NonNegativeReals)
#model.test = Var(bounds=(0,0.5),initialize=1)
######################################
#Constraints
######################################

def s_cb_rule(model, i, t, r, w):
	return model.xs[i,t,r,w]<=model.xs[i,t+1,r,w]
model.cSb = Constraint(model.sI, model.sT_End, model.sR, model.Scen, rule=s_cb_rule)

def s_cc_rule(model, i, t, r, w):
	return model.xs[i,t+1,r+1,w]<=model.xs[i,t,r,w]
model.cSc = Constraint(model.sI, model.sT_End, model.sR_End, model.Scen, rule=s_cc_rule)

def s_cd_rule(model, i, t, w):
	#return sum(model.xs[i,t,r]-model.xs[i,t-1,r] for r in model.sR)<=model.zs[t]
	return sum(model.ws[i,t,r,w] for r in model.sR)<=model.zs[t,w]
model.cSd = Constraint(model.sI, model.sT_0, model.Scen, rule=s_cd_rule)

def s_ce_rule(model, i, w):
	return model.xs[i,0,1,w]<=model.zs[0,w]
model.cSe = Constraint(model.sI,model.Scen, rule=s_ce_rule)

#def s_ce1_rule(model):
#	return model.z == model.zs[0]
#model.cSe1 = Constraint(rule=s_ce1_rule)

def s_cf_rule(model, i, r, t,w):
	return model.xs[i,t,r,w]<=model.xs[i,t+model.pLT[i,r+1,w],r+1,w]
model.cSf = Constraint(model.scf, rule=s_cf_rule)

def s_cg_rule(model, i, w):
	return model.xs[i,model.pLT[i,1,w],1,w]==1
model.cSg = Constraint(model.scg, rule=s_cg_rule)

def s_ch_rule(model, i, r, w):
	return model.xs[i,0,r,w]==0
model.cSh = Constraint(model.sI, model.sR_0, model.Scen, rule=s_ch_rule)

def s_ci_rule(model, i, w):
	return model.xs[i,0,1,w]==model.x[i]
model.cSi = Constraint(model.sI, model.Scen, rule=s_ci_rule)

def s_cj_rule(model, i):
	return model.x[i] >= model.pKesi[i]
model.cSj = Constraint(model.sI, rule=s_cj_rule)

def s_ck_rule(model, i, t, r, w):
	return model.ws[i,t,r,w] == model.xs[i,t,r,w]-model.xs[i,t-1,r,w]
model.cSk = Constraint(model.sI, model.sT_0, model.sR, model.Scen, rule=s_ck_rule)

def s_cl_rule(model, i, r, w):
	return model.ws[i,0,r,w]==model.xs[i,0,r,w]
model.cSl = Constraint(model.sI, model.sR, model.Scen, rule=s_cl_rule)

def s_cm_rule(model, i, r, t, w):
	return model.ws[i,t,r,w]==model.ws[i,t-model.pT.value,r,w]
model.cSm = Constraint(model.scm, rule=s_cm_rule)

def s_cn_rule(model, i, r, t, w):
	return model.ws[i,t,r,w]==0
model.cSn = Constraint(model.scn, rule=s_cn_rule)

def s_co_rule(model, i, t, r, w):
	return model.u[i,t,r,w]-model.v[i,t,r,w]==\
			model.ws[i,t+model.pLT[i,r,w],r,w]-model.ws[i,t,r-1,w]
model.cSo = Constraint(model.sI, model.sT, model.sR_0, model.Scen, rule=s_co_rule)

def s_cp_rule(model, i, t, w):
	return model.u[i,t,1,w]-model.v[i,t,1,w]==\
			model.ws[i,t+model.pLT[i,1,w],1,w]-1
model.cSp = Constraint(model.sI, model.sT, model.Scen, rule=s_cp_rule)

def s_cq_rule(model, i, t, r, w):
	return model.u[i,t,r,w]+model.v[i,t,r,w]<=1
model.cSq = Constraint(model.sI, model.sT, model.sR, model.Scen, rule=s_cq_rule)

def s_cr_rule(model, i, t, r, w):
	return model.u[i,t,r,w]<=1 #and model.u[i,t,r]>=0
model.cSr = Constraint(model.sI, model.sT, model.sR, model.Scen, rule=s_cr_rule)

def s_cs_rule(model, i, t, r, w):
	return model.v[i,t,r,w]<=1 #and model.v[i,t,r]>=0
model.cSs = Constraint(model.sI, model.sT, model.sR, model.Scen, rule=s_cs_rule)

def s_cs1_rule(model, i, r, w):
	return sum(model.u[i,t,r,w]+model.v[i,t,r,w] for t in model.sT) <= 2
model.cSs1 = Constraint(model.sI, model.sR_0, model.Scen, rule=s_cs1_rule)

def s_ct_rule(model, i, t, r, w):
	return model.xs[i,t,r,w]<=1 #and model.xs[i,t,r]>=0
model.cSt = Constraint(model.sI, model.sT, model.sR, model.Scen, rule=s_ct_rule)

def s_cu_rule(model, t, w):
	return model.zs[t, w]<=1 #and model.zs[t]>=0
model.cSu = Constraint(model.sT, model.Scen, rule=s_cu_rule)

def s_cv_rule(model, i, t, r, w):
	return model.ws[i,t,r,w]<=1 #and model.ws[i,t,r]>=0
model.cSv = Constraint(model.sI, model.sT_ex, model.sR, model.Scen, rule=s_cv_rule)

def s_cw_rule(model, i):
	return model.x[i]<=1 
model.cSw = Constraint(model.sI, rule=s_cw_rule)
"""
"""
######################################
#Objective
######################################
def exp_cost_rule(model):
	model.idr1   = sum(model.pCCR[i]*model.ws[i,model.pLT[i,1,w],1,w]+\
			 model.pCPR[i]*(1-model.ws[i,model.pLT[i,1,w],1,w])\
			 for i in model.sI for w in model.Scen)
	model.othIdr = sum(sum(model.pCCR[i]*(1-0.5*sum(model.u[i,t,r,w]+model.v[i,t,r,w] for t in model.sT))-\
					model.pCCR[i]*(1-0.5*(model.xs[i,model.pT,r,w]+model.xs[i,model.pT,r-1,w]))+\
					model.pCPR[i]*0.5*sum(model.u[i,t1,r,w]+model.v[i,t1,r,w] for t1 in model.sT)\
					for r in model.sR_0)-0.5*model.pCPR[i]\
					for i in model.sI for w in model.Scen)
	model.setup  = sum(model.pd*model.zs[t,w] for t in model.sT for w in model.Scen)
	return (model.idr1+model.othIdr+model.setup)/model.NUMSCEN
model.objCost = Objective(rule=exp_cost_rule,sense=minimize)
