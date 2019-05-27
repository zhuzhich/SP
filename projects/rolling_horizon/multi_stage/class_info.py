
import math
#####################################
#class info
#####################################
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
			#print (ageFrom)
			if tmp1 == 0:
				prob = 1;
			else:
				prob = float(tmp)/tmp1;	
		return prob;
	
	def __init__(self, idx, w_shape, w_scale, initAge,initFail, intvl, cCM, cPM, cS):
		self.idx = idx;
		self.w_shape = w_shape;
		self.w_scale = w_scale;
		self.initAge = initAge;
		self.initFail = initFail;
		self.intvl = intvl;
		self.cCR = cCM;
		self.cPR = cPM;
		self.cS = cS;

		
		
#system information
#parameters
class system_info():
	def add_com(self, comInfo):
		self.comInfoAll.append(comInfo);
	def __init__(self, N, T, inspInterval, cS):
		self.nComponents = N;
		self.nStages = T;
		self.inspItvl = inspInterval;
		self.cS = cS;
		#self.cInsp = cInsp;
		self.comInfoAll = [];
		#output
		self.N0 = [];
		self.N1 = [];
		self.Nu = [];
		self.Nu_0 = [];
		self.Nu_1 = [];
		self.time = 0;
		self.objValue = 0;
		self.iterInfo = [];
		