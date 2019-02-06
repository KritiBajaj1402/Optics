import matplotlib.pyplot as plt
# import numpy as np

f1 = open("reachabilityDist.txt")
line = f1.readline()
y = list()
while line:
	dist = (float)(line.strip())
	y.append(dist)
	line = f1.readline()
f1.close()
print len(y)

# y.sort()
max_rd = 0
ind_to_be_modified = list()
for i in xrange(len(y)):
	if y[i] < 1e9 - 10 and y[i] > max_rd:
		max_rd = y[i]
	elif y[i] > 1e9 - 10:
		ind_to_be_modified.append(i)

print "len(ind_to_be_modified)",
print len(ind_to_be_modified)

print "max_rd", max_rd

for i in ind_to_be_modified:
	y[i] = 2*max_rd


# y = y[0:5]
# print y

# x = np.arange(len(y))
x = range(len(y))
plt.bar(x, y)
# plt.plot(x, y)
plt.xlabel("Cluster-order of the objects")
plt.ylabel("Reachability distance")
plt.title("Reachability Distance Plot")
plt.show()