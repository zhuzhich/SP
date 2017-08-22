def instance_info(instance):
	print "==================xs===================="
	for i in instance.sI:
		print "\n i=%d" %(i)
		for r in instance.sR:
			print "\n r=%d" %(r)
			for t in instance.sT:
				#print("xs(%d,%d,%d)=%d" %(i,t,r,instance.xs[i,t,r]() )),
				print int(instance.xs[i,t,r]()),
	print("")
	print "==================ws===================="
	for i in instance.sI:
		print "\n i=%d" %(i)
		for r in instance.sR:
			print "\n r=%d" %(r)
			for t in instance.sT_ex:
				print int(instance.ws[i,t,r]()),
	print("")
	print "==================zs===================="
	for t in instance.sT_0:
		print int(instance.zs[t]()),
	print("")
	print "==================(u,v)===================="
	for i in instance.sI:
		print "\n i=%d" %(i)
		for r in instance.sR:
			print "\n r=%d" %(r)
			for t in instance.sT:
				#print("xs(%d,%d,%d)=%d" %(i,t,r,instance.xs[i,t,r]() )),
				print "(%d,%d)" %(int(instance.u[i,t,r]()),int(instance.v[i,t,r]())),
	
	print "\n==================first individual================="
	print sum(instance.pCCR[i]*instance.ws[i,instance.pLT[i,1],1]()+instance.pCPR[i]*(1-instance.ws[i,instance.pLT[i,1],1]())
		for i in instance.sI)
	print "\n==================first stage================="	
	print sum(instance.pCPR[i1]*instance.x[i1]() for i1 in instance.sI)+\
			sum((instance.pCCR[i2]-instance.pCPR[i2])*instance.pKesi[i2] for i2 in instance.sI)
	print "\n==================others================="	
	for i in instance.sI:
		print ("i=%d"%i)
		print sum(instance.pCCR[i]*(1-0.5*sum(instance.u[i,t,r]()+instance.v[i,t,r]() for t in instance.sT))-\
				instance.pCCR[i]*(1-0.5*(instance.xs[i,instance.pT,r]()+instance.xs[i,instance.pT,r-1]))()+\
				instance.pCPR[i]*0.5*sum(instance.u[i,t1,r]()+instance.v[i,t1,r]() for t1 in instance.sT)\
				for r in instance.sR_0)-0.5*instance.pCPR[i]			
	print "\n==================z================="						
	print sum(instance.pd*instance.zs[t]() for t in instance.sT_0)
