
import math
import random
import numpy as np
#####################################
#class info
#####################################

class heuristic_info():
	def d_print(self):
		print ("currentWin");
		print (self.currentWin);
		print ("shadowWin");
		print (self.shadowWin);
		print ("lastIndTime");
		print (self.lastIndTime);
		print ("workingInd");
		print (self.workingInd);
		print ("sol");
		print (self.sol);
		print ("cost");
		print (self.cost);
	
	def update_i(self, sysInfo, scen, i, t1, repTime):
		# replace i at repTime
		iIdx = self.currentWin[i][1];
		rIdx = self.workingInd[iIdx];
		LT = sysInfo.comInfoAll[iIdx].LT[scen][rIdx];
		repTime = int(repTime);
		if LT == 0:	#fail
			self.sol[iIdx][repTime] = 2;
		else:
			failTime = max(self.lastIndTime[iIdx],0) + LT;
			if repTime == failTime:
				self.sol[iIdx][repTime] = 2;
			elif repTime < failTime:
				self.sol[iIdx][repTime] = 1;
			else:
				self.d_print();
				print (i, iIdx, rIdx, LT, failTime);
				print("error! dangerous error!!");
		# update the first working individual's info
		self.lastIndTime[iIdx] = repTime;
		self.workingInd[iIdx] += 1;
		rIdx = self.workingInd[iIdx];

		LT = sysInfo.comInfoAll[iIdx].LT[scen][rIdx];
		self.shadowWin[i][0] = repTime + \
										max(1, LT-t1);	


	def sort_currentWin(self):
		self.currentWin.sort(key=lambda x:x[0]);
	
	def __init__(self, sysInfo):
		#sorted lifetime,
		#first column is lifetime;
		#second colum is the component index
		self.currentWin = [[0]*2 \
						for i in range(sysInfo.nComponents)];
		self.shadowWin = [[0]*2 \
						for i in range(sysInfo.nComponents)];
		
		
		self.lastIndTime = [-1]*sysInfo.nComponents;
		self.workingInd = [0]*sysInfo.nComponents;
		#global optimal soution
		self.sol = [[0]*sysInfo.nStages \
					for i in range(sysInfo.nComponents)];
		self.cost = 0;

class component_info():
	def weibull_cdf(self, t):
		tmp = float(t)/self.w_scale;
		tmp1 = tmp ** self.w_shape;
		return 1 - math.exp(-tmp1);
		
	def cond_fail_prob(self, ageFrom, ageTo):
		if ageTo - ageFrom != 1:
			prob = 0;
		else:
			tmp = self.weibull_cdf(ageTo) -  self.weibull_cdf(ageFrom);
			tmp1 = 1 -  self.weibull_cdf(ageFrom);
			prob = float(tmp)/tmp1;	
		return prob;
	
	
	def create_LifeTime(self, ranSeed):
		# if self.LT not empty, only update the first individual's lifetime
		# e.g., fail -> 0 else -> 1 or more
		#		
		if len(self.LT) > 0:
			LT_in_arr = np.array(self.LT);
			LT_arr = LT_in_arr[:,range(int(self.initFail==1), self.nIndividuals+ int(self.initFail==1)) ];
			for idx_w in range(self.nScenarios):
				if self.initFail == 1:
					LT_arr[idx_w,0] = 0;	#force to replace the first individual	
				else:
					random.seed(self.idx+idx_w+ranSeed); ###control the seed  				
					ran_num = random.uniform(0,1);
					part1 = math.log(ran_num);
					s_inv = 1.0/self.w_shape;
					surv_time = self.initAge;
					part2 = (surv_time/self.w_scale)**self.w_shape;
					part3 = part2 - part1;
					tmp = int((part3**s_inv)*self.w_scale) - surv_time;	
					##############
					#round - > int
					##############					
					LT_arr[idx_w,0] = int(max(1,tmp));
			self.LT = LT_arr.tolist();
		
		######for the first time######################
		else:
			for idx_w in range(self.nScenarios):
				tmp = [0]*self.nIndividuals;
				for idx2 in range(self.nIndividuals):
					random.seed(self.idx + idx2+idx_w+ranSeed);###control the seed  
					ran_num = random.uniform(0,1);
					if idx2 == 0:
						if self.initFail == 1:
							tmp[0]  = 0	#force to replace the first individual	
						else:
							part1 = math.log(ran_num);
							s_inv = 1.0/self.w_shape;
							surv_time = self.initAge;
							part2 = (surv_time/self.w_scale)**self.w_shape;
							part3 = part2 - part1;
							tmp[idx2] = max(1,round((part3**s_inv)*self.w_scale) - surv_time);	
							tmp[idx2] = int(tmp[idx2]);
					else:
						ran_num_log = -math.log(ran_num)
						s_inv = 1.0/self.w_shape
						LT1 = round((ran_num_log**s_inv)*self.w_scale)
						##############
						#round - > int
						##############
						tmp[idx2] = int(max(1,LT1))	
				#LT_tmp <- tmp
				self.LT.append(tmp);
		
		#print ("LT", self.LT);


	def __init__(self, idx, w_shape, w_scale, initAge,initFail, \
				intvl, cCM, cPM, cS, nScenarios, nIndividuals):
		self.idx = idx;
		self.w_shape = w_shape;
		self.w_scale = w_scale;
		self.initAge = initAge;
		self.initFail = initFail;
		self.intvl = intvl;
		self.cCR = cCM;
		self.cPR = cPM;
		self.cS = cS;
		self.nScenarios = nScenarios;
		self.nIndividuals = nIndividuals;
		self.LT = [];

		
		
#system information
#parameters
class system_info():
	def add_com(self, comInfo):
		self.comInfoAll.append(comInfo);
	def __init__(self, N, T, inspInterval, cS, nScenarios, ranSeed):
		self.nComponents = N;
		self.nStages = T;
		self.inspItvl = inspInterval;
		self.cS = cS;
		self.nScenarios = nScenarios;
		self.ranSeed = ranSeed;
		#self.cInsp = cInsp;
		self.comInfoAll = [];
		#output
		self.solX = [];
		self.N0 = [];
		self.N1 = [];
		self.Nu = [];
		self.Nu_0 = [];
		self.Nu_1 = [];
		self.time = 0;
		self.objValue = 0;
		self.iterInfo = [];

		