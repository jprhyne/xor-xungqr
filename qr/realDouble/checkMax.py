#!/bin/env python

f = open("results.txt")

maxOrth = 0
maxRepr = 0

f.readline()
for line in f:
    a = line.split()
    maxRepr = max(maxRepr,float(a[0]))
    maxOrth = max(maxOrth,float(a[1]))

print("Maximum representation error: " + str(maxRepr))
print("Maximum orthogonal error: " + str(maxOrth))
f.close()
