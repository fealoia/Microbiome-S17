import matplotlib.pyplot as plt
import sys

#Visualization of aligned reads along DNA sequence

f_name = sys.argv[1]
length = sys.argv[2]

data = []
with open(f_name, 'r') as fd:
        for line in fd:
                data.append(line.split(':')[0]) 

plt.figure(figsize=(12,9))
ax = plt.subpot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

plt.yticks(fontsize=12)
plt.xticks(range(0,length+1),fontsize=12)

plt.xlabel("Read Position",fontsize=14)
plt.ylabel("Non_Concat Only Reads Frequency",fontsize=14)

plt.hist(data,bins=200)
plt.show()		
