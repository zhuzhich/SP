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

def create_files(iter):

	file_name = "C:\\Users\\Yisha\\Desktop\\ZZC\\extensive_form\\ef.dat"
	f = open(file_name,"w")


	W = 50;
	I = 3;
	T = 7;
	T_ex = 50;
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
for iter in range(1):
	print ("================")
	print ("iter="+str(iter));
	create_files(iter);
	
	# initialize the instances.
	ef_mdl = import_file("ef.py").model
	ef_insts=[]
	ef_insts.append(ef_mdl.create_instance(name="ef_instance", filename="ef.dat"))
	solver_manager = SolverManagerFactory("serial")
	start_time1 = time.clock()
	solve_all_instances(solver_manager, 'cplex', ef_insts)
	end_time1 = time.clock()
	print ("solving time = ",str(end_time1-start_time1))
	"""
	for i in ef_insts[0].sI:
		for t in ef_insts[0].sT:
			print ("")
			for r in ef_insts[0].sR:
				print ("(%d,%d,%d)=(%d,%d)" %(i,t,r,ef_insts[0].u[i,t,r,1](),ef_insts[0].v[i,t,r,1]())),
	"""			

	for i in ef_insts[0].sI:
		print("w_shape["+str(i)+"]="+str(round(ef_insts[0].w_shape[i], 4)))	
		print("x["+str(i)+"]="+str(round(ef_insts[0].x[i](), 4)))
	print ("objective value=%f" %ef_insts[0].objCost())
		
