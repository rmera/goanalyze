#/usr/bin/env python

from os.path import expanduser
import sys, os

print("This is the goAnalyze! install system")

version="v0.1.0"

OS=["Linux","Darwin","Windows"]

packed_files=["LICENSE","README","INSTALL.py","gochem.py","goanalyze.py","gomd.py"]

base=os.getcwd()  #all of it are belong to us

os.system("mkdir packed")


for i in  OS:
	base=os.getcwd()
	dirname="packed/"+i
	os.system("mkdir %s"%(dirname))
	os.system("cp -R bin_%s %s"%(i,dirname))
	for j in packed_files:
		os.system("cp %s %s/"%(j,dirname))
	os.chdir(dirname)
	os.system("tar cvzf ../goAnalyze-%s-%s.tgz *"%(version,i))
	os.chdir(base)

