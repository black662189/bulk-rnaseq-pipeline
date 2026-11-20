import sys
count ={}
for line in open(sys.argv[1],'r'):
	if line.startswith('gene'):continue
	line =line.strip().split()
	if line[0] not in count:
		count[line[0]]=1
	else:
		count[line[0]]+=1
	if line[1] not in count:
		count[line[1]]=1
	else:
		count[line[1]]+=1
fd={}
for line in open(sys.argv[2],'r'):
	line=line.strip().split('\t')
	fd[line[0]]=line[-5]
for line in count:
	print "%s\t%s\t%s"%(line,count[line],fd[line])

	
	
	
	
	
