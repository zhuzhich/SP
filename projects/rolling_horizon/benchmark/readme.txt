implement paper:
A dynamic policy for grouping maintenance activities
EJOR 1997, Wildeman et al.

Code author: Zhicheng Zhu
Email: zzhu3@lamar.edu, zhich.zhu@gmail.com
Last Update: 03/21/2018


Code structure:
===============
./main.py:
-----------------
    main file.
./class_info.py:
----------------
	all 4 classes in main.py:
	1. system_parameter: system static information
	2. component_parameters: component static information
	3. component_running: component running information. Dynamically changing in program
	4. system_running: system running information. Dynamically changing in program
./result.txt:
-------------
	results showed in paper: xxxxx
	

Notice:
=======
1. It's using age-based maintenance + perfect maintainenace here.
   Reference paper: block-based maintenance + minimum repair
2. It's a rolling horizon version.
 