import sys

#Prints total number of aligned reads for each species

file1 = sys.argv[1]
ids = {}
with open(file1, 'rb') as fd1:
	for line in fd1:
		unique = line.split("	")[2]
		if unique in ids:
			ids[unique] += 1
		else:
			ids[unique] = 1	

#print pairs
for k,v in sorted(ids.items()):
        print(str(k) + ': ' + str(v))
