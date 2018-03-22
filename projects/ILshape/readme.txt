This code trys to implement integer L-shape method to solve stochastic grouping problem in multi-component system.
Which is the algorithm 2 in paper:
xxxxxxxxx
Code author: Zhicheng Zhu
Email: zzhu3@lamar.edu, Zhich.zhu@gmail.com


Code structure:
===============
./data/:
------- 
	data file for each scenario. Named as ScenNode[x].dat, where [x] is the scenario index, starts from 1.
./linux_run/: 
-------------
	the job script for HPC (high performance computer). Ignore it.
./otherFiles/:
------------
	main_multi: the file inside the multi-cut version of integer L-shape. Copy it to root folder for using.
	master.py (master_Int.py):relaxation (integer) version of master problem implementation
./result/:
----------
	result files
./main_single.py:
-----------------
    main file for single-cut implementation.
./genFile.py:
-------------
	generate scenario-specific data file in ./data
./sub.py (subIP.py):
--------------------
	relaxation (integer) version of sub-problem implementation.

Notice:
=======
1. We are using python-CPLEX to solve master problem by using branch-and-cut.
2. master problem is built by python-CPLEX.
3. B&C is embeded in CPLEX. So master problem is solved as integer problem.
4. Sub-problem is solved via Pyomo package, solver is CPLEX
 