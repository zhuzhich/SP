Deterministic dxtensive form (DEF) to solve stochastic grouping problem in multi-component system.

For paper: xxxxxxxxx

Code author: Zhicheng Zhu
Email: zzhu3@lamar.edu, zhich.zhu@gmail.com
Last Update: 03/21/2018

Code structure:
===============
./run.py:
--------- 
	main script file
./ef.py:
-------- 
	extensive form model
./ef.dat: 
-------
	extensive form data.
./result/:
----------
	result files

Notice:
=======
1. No script to run multiple cases, i.e, different number of component and/or scenarios and/or time horizon.
2. For each experiment, one only need to amend file ./ef.dat for different parameters.
3. "REMEMBER": number of indidividuals changes along with time time horizon in ef.dat, i.e., R = T + 2
4. Pyomo version: 5.5.0, DO NOT USE 5.6+
5. Important paths: lib/sitepackage/pyomo,  pkgs/pyomo-[version], the package is in folder ../../../sp_support
6. CPLEX version that still working: 12.8.00
7. Copy C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\python\3.6\x64_win64\cplex to C:\Users\Zhicheng\Anaconda3\Lib\site-packages;
 