# Author Zhicheng Zhu
# Email: zzhu3@lamar.edu

# solve extensive form model using Pyomo
# 
# Structure:
# run.py: main script file
# ef.py:  extensive form model
# ef.dat: extensive form data.

#
#Last Update: 03/21/2018
#
from pyutilib.misc import import_file
from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.opt.base import SolverFactory
from pyomo.opt.parallel import SolverManagerFactory
from pyomo.opt.parallel.manager import solve_all_instances
import time
import math
import random
import numpy as np
class system_info():
	def weibull_cdf(self, i, t):
		tmp = float(t)/self.w_scale[i];
		tmp1 = tmp ** self.w_shape[i];
		return 1 - math.exp(-tmp1);
		
	def cond_fail_prob(self, i, ageFrom, ageTo):
		if (ageTo - ageFrom) != 1:
			prob = 0;
		else:
			tmp = self.weibull_cdf(i, ageTo) -  self.weibull_cdf(i, ageFrom);
			tmp1 = 1 -  self.weibull_cdf(i, ageFrom);
			if tmp1 == 0:
				prob = 1;
			else:
				prob = float(tmp)/tmp1;	
		return prob;


	def __init__(self, N, T, inspInterval, cS, age,\
				failState, cPR, cCR, w_shape, w_scale):
		self.nComponents = N;
		self.nStages = T;
		self.inspItvl = inspInterval;
		self.cS = cS;
		self.age = age;
		self.failState = failState;
		self.cPR = cPR;
		self.cCR = cCR;
		self.cS = cS;
		self.w_shape = w_shape;
		self.w_scale = w_scale;

		#output


def create_files(randSeed, sysInfo):

	file_name = "C:\\Users\\Yisha\\Desktop\\ZZC\\DEF_roll\\ef.dat"
	f = open(file_name,"w")


	W = 100;
	T_ex = 50;

	f.write("param NUMSCEN := " + str(W)+";\n\n");
	
	f.write("param prob := " + str(1.0/W)+";\n\n")	;
	
	f.write("param pI := " + str(sysInfo.nComponents)+";\n\n");

	f.write("param pT := " + str(sysInfo.nStages)+";\n\n");

	f.write("param pT_ex := " + str(T_ex)+";\n\n");
	
	f.write("param pd := " + str(sysInfo.cS)+";\n\n");

	f.write("param pR := " + str(sysInfo.nStages+2)+";\n\n")
	
	f.write("param randSeed := " + str(randSeed)+";\n\n")

	text1 = "param ps := ";
	text2 = "param pKesi := ";
	text3 = "param pCCR := ";
	text4 = "param pCPR := ";
	text5 = "param w_shape := ";
	text6 = "param w_scale := ";
	for i in range(sysInfo.nComponents):
		text1 = text1 + str(i+1) + " " + str(sysInfo.age[i]) + " ";
		text2 = text2 + str(i+1) + " " + str(sysInfo.failState[i]) + " ";
		text3 = text3 + str(i+1) + " " + str(sysInfo.cCR[i]) + " ";
		text4 = text4 + str(i+1) + " " + str(sysInfo.cPR[i]) + " ";
		text5 = text5 + str(i+1) + " " + str(sysInfo.w_shape[i]) + " ";
		text6 = text6 + str(i+1) + " " + str(sysInfo.w_scale[i]) + " ";
	
	text1 = text1 + ";\n\n";
	text2 = text2 + ";\n\n";
	text3 = text3 + ";\n\n";
	text4 = text4 + ";\n\n";
	text5 = text5 + ";\n\n";
	text6 = text6 + ";\n\n";
	
	f.write(text1);
	f.write(text2);	
	f.write(text3);	
	f.write(text4);
	f.write(text5);
	f.write(text6);	
	'''
## create lifetimes
	for idx in range(0,w_Lo):
		file_name = dir_local + "\\nodedata\\ScenNode"+str(idx+1)+".dat"
		f = open(file_name,"w")
		f.write("#pLT is two-dimensional [i,r]=[component, individual]\n#[i,*] means all data in row i.\n\n")
		f.write("param pLT := ")
		for idx1 in range(0,I_Lo):
			i = idx1 + 1
			f.write("\n["+str(i)+",*]\n")
			for idx2 in range(0,R_Lo):
				r = idx2 + 1
				random.seed(i+r+idx+1) ###control the seed
				if r == 1:
					if pKesi_Lo[idx1] == 1:
						LT = 0
					else:    
						ran_num = random.uniform(0,1)
						ran_num_log = -math.log(ran_num)
						s_inv = 1.0/w_shape_Lo[idx1]	
						LT1 = int((ran_num_log**s_inv)*w_scale_Lo[idx1]) - s_Lo					
						LT = max(1,LT1)
				else:
					ran_num = random.uniform(0,1)
					ran_num_log = -math.log(ran_num)
					s_inv = 1.0/w_shape_Lo[idx1]
					LT1 = int((ran_num_log**s_inv)*w_scale_Lo[idx1])
					LT = max(1,LT1)
				f.write(str(r)+ " " + " "*(2-len(str(r))) + str(LT)+ "    ")
				#f.write(str(r)+ " " + str(LT)+ "   ")
				if r%10 == 0:
					f.write("\n")
			f.write("\n")
		f.write(";")
		f.close()
'''
	f.close();

##
#start from here!
nComponents = 3;
nStages = 6;		#total stages - 1, t=0, 1, 2, ..., T
inspInterval = 1;
cS = 5;

failState = [0]*nComponents;
initFailSate = [0]*nComponents;
age = [0]*nComponents;
initAge = [0]*nComponents;
cPR = [0]*nComponents;
cCR = [0]*nComponents;
w_shape = [0]*nComponents;
w_scale = [0]*nComponents;

for i in range(nComponents):
	#kesi
	if i==0:
		failState[i] = 1;
		initFailSate[i] = 1;
	#age
	age[i] = 2;
	initAge[i] = 2;
	#cPR
	cPR[i] = 1;
	#cCR
	random.seed(i*30);
	temp = random.uniform(6,16);
	cCR[i] = round(temp,1);			
	#shape
	random.seed(i*20)   
	temp = random.uniform(4,7)
	w_shape[i] = round(temp,1);			
	#scale
	random.seed(i*10)   
	temp = random.uniform(1,8)#(4,11)
	w_scale[i] = round(temp,1);		



sysInfo = system_info(nComponents, nStages, inspInterval, cS, age,\
				failState, cPR, cCR, w_shape, w_scale);


rep = 100;
cost = []	

for randSeed in range(rep):	
	cost_rep = 0;
	sysInfo.nStages = nStages;
	for i in range(sysInfo.nComponents):
		sysInfo.age[i] = initAge[i];
		sysInfo.failState[i] = initFailSate[i];
	print ("rep");
	print ("=====")
	print (randSeed);
	for t in range(nStages+1):
		#print ("t, age, failState, stages");
		#print (t, sysInfo.age, sysInfo.failState, sysInfo.nStages);
		#print ("current cost");
		#print (cost_rep);
		if t == nStages:
			tmp = 0;
			for i in range(sysInfo.nComponents):
				tmp += sysInfo.cCR[i] * sysInfo.failState[i];
			if tmp > 0:
				tmp += sysInfo.cS;
			cost_rep += tmp;
		else:
			create_files(randSeed, sysInfo);
			# initialize the instances.
			ef_mdl = import_file("ef.py").model
			ef_insts=[]
			ef_insts.append(ef_mdl.create_instance(name="ef_instance", filename="ef.dat"))
			solver_manager = SolverManagerFactory("serial")
			#start_time1 = time.clock()
			solve_all_instances(solver_manager, 'cplex', ef_insts);
			# get first stage solution
			x = [];
			tmp = 0;
			for i in ef_insts[0].sI:
				x.append(round(ef_insts[0].x[i](),0));
				tmp += sysInfo.cPR[i-1]*x[i-1];
				tmp += (sysInfo.cCR[i-1] - sysInfo.cPR[i-1])*sysInfo.failState[i-1];
			if tmp > 0:
				tmp += sysInfo.cS;
			cost_rep += tmp;
			#print ("solution");
			#print (x);
			#prepare next stage:
			age = [sysInfo.age[i]*(1-x[i]) for i in range(sysInfo.nComponents)];	#age after maintenance
			failProb = [sysInfo.cond_fail_prob(i, age[i], age[i]+1) for i in range(sysInfo.nComponents)];

			randProb = [random.uniform(0,1) for i in range(sysInfo.nComponents)];
			failState = [int(randProb[i]<failProb[i]) for i in range(sysInfo.nComponents)];
			#print (kesi);
			for i in range(sysInfo.nComponents):
				sysInfo.failState[i] = failState[i];	
				sysInfo.age[i] = age[i] + 1;
				sysInfo.nStages = nStages - t - 1;
	cost.append(cost_rep);
print ("cost");
print (cost);
print ("mean cost");
print (np.mean(cost));
print ("variance cost");
print (np.var(cost));	
print ("max and min");
print (max(cost), min(cost));
#print ("objective value=%f" %ef_insts[0].objCost())
	
#print ("====debug===");
#print (ef_insts[0].w_shape());	
#for i in ef_insts[0].sI:
#	print("w_shape["+str(i)+"]="+str(round(ef_insts[0].w_shape[i], 4)))	
#	print("x["+str(i)+"]="+str(round(ef_insts[0].x[i](), 4)))	
			
			
			
#end_time1 = time.clock()
#print ("solving time = ",str(end_time1-start_time1))
"""
for i in ef_insts[0].sI:
	for t in ef_insts[0].sT:
		print ("")
		for r in ef_insts[0].sR:
			print ("(%d,%d,%d)=(%d,%d)" %(i,t,r,ef_insts[0].u[i,t,r,1](),ef_insts[0].v[i,t,r,1]())),
for i in ef_insts[0].sI:
	print("x["+str(i)+"]="+str(round(ef_insts[0].x[i](), 4)))
print ("objective value=%f" %ef_insts[0].objCost())
"""

