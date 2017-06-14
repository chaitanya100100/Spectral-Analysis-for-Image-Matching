import os
import sys

file_list = os.listdir('.')

for f in file_list:
	if ".mat" in f:
		k = f[:-4]
		print k+"/"+f
		os.mkdir(k)
		os.rename(f, k+"/"+f)

