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
import numpy as np

def create_files(iter,W):

	file_name = "C:\\Users\\Yisha\\Desktop\\ZZC\\extensive_form\\ef.dat"
	f = open(file_name,"w")


	#W = 140;
	I = 2;
	T = 49;
	T_ex = 100;
	d = 5;
	f.write("param NUMSCEN := " + str(W)+";\n\n")
	
	f.write("param prob := " + str(1.0/W)+";\n\n")	
	
	f.write("param pI := " + str(I)+";\n\n")

	f.write("param pT := " + str(T)+";\n\n")

	f.write("param pT_ex := " + str(T_ex)+";\n\n")

	f.write("param ps := " + str(2)+";\n\n")

	f.write("param pR := " + str(T+2)+";\n\n")

	f.write("param pd := " + str(d)+";\n\n")
	
	f.write("param iter := " + str(iter)+";\n\n")

	f.close()

cost = [];
start_time1 = time.clock()

for iter in range(100):
	print ("================")
	print ("iter="+str(iter));
	W = 1;
	create_files(iter,W);
	
	# initialize the instances.
	ef_mdl = import_file("ef.py").model
	ef_insts=[]
	ef_insts.append(ef_mdl.create_instance(name="ef_instance", filename="ef.dat"))
	solver_manager = SolverManagerFactory("serial")
	solve_all_instances(solver_manager, 'cplex', ef_insts)
	#print ("solving time = ",str(end_time1-start_time1))
	"""
	for i in ef_insts[0].sI:
		for t in ef_insts[0].sT:
			print ("")
			for r in ef_insts[0].sR:
				print ("(%d,%d,%d)=(%d,%d)" %(i,t,r,ef_insts[0].u[i,t,r,1](),ef_insts[0].v[i,t,r,1]())),
	"""			
	#for i in ef_insts[0].sI:
	#	for w in ef_insts[0].Scen:
	#		tmp = [];
	#		for r in ef_insts[0].sR:
	#			tmp.append(ef_insts[0].pLT[i,r,w]);
	#		print (tmp)

	for i in ef_insts[0].sI:
		#print("w_shape["+str(i)+"]="+str(round(ef_insts[0].w_shape[i], 4)))
		#print("w_scale["+str(i)+"]="+str(round(ef_insts[0].w_scale[i], 4)))
		#print("cCR["+str(i)+"]="+str(round(ef_insts[0].pCCR[i], 4)))	
		#print("LT["+str(i)+"]="+str(round(ef_insts[0].pLT[i,:,:], 4)))			
		print("x["+str(i)+"]="+str(round(ef_insts[0].x[i](), 4)))
		for t in ef_insts[0].sT:
			#print ("-----t------",t)
			tmpU = [];
			tmpV = [];
			tmpWs = [];
			tmpZs = [];
			tmpXs = [];
			for r in ef_insts[0].sR:
				#print ("-----r------",r)
				#print("LT", ef_insts[0].pLT[i,r,1])
				for s in ef_insts[0].Scen:
					#tmpU.append(ef_insts[0].u[i,t,r,s]());
					#tmpV.append(ef_insts[0].v[i,t,r,s]());
					tmpWs.append(ef_insts[0].ws[i,t,r,s]());
					#tmpXs.append(ef_insts[0].xs[i,t,r,s]());
					#tmpZs.append(ef_insts[0].zs[t,s]());
			#print("u",tmpU);
			#print("v",tmpV);
			#if sum(tmpWs) > 0:
				#print("t",t+1);
			#print("xs",tmpXs);
			#print("zs",tmpZs);
	print ("objective value=%f" %ef_insts[0].objCost())
	cost.append(ef_insts[0].objCost());
end_time1 = time.clock()
print ("avg time=%d" %((end_time1 - start_time1)/100.0));
print ("average cost");
print (np.average(cost))
print ("cost", cost)

		
