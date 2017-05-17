import sys

file1 = sys.argv[1]
file2 = sys.argv[2]
file3 = sys.argv[3]
ids = []
count1 = 0
count2 = 0

#file1 = SAM to obtain reference IDs
#file2 = fastq to match references
#file3 = output 

fdw = open(file3, 'a')

with open(file1, 'rb') as fd1:
	for line in fd1:
		unique = line.split("	")[0]
		ids.append(unique)

with open(file2, 'rb') as fd2:
	count = 0
	line = fd2.readline()
	while line:	
		if(count%4==0):
			if line.split(' ')[0].split('@')[1] in ids:
				fdw.write(line)
				fdw.write(fd2.readline())
				fdw.write(fd2.readline())
				fdw.write(fd2.readline())
				count += 3
		count += 1
		line = fd2.readline()
