import sys

n = 1024
q = 8404993

for j in range(q):
	if(pow(j,n,q) == (q-1)):
		print(j)
		break

