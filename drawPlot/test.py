import numpy as np

second = []
third = []
fourth = []


with open('test.csv', 'r') as f:
	lines = f.readlines()
	for line in lines:
		values = line.strip().split()
		second.append(float(values[1]))
		third.append(float(values[2]))
		fourth.append(float(values[3]))

print("2nd : ", second)
print("3rd : ", third)
print("4th : ", fourth)
