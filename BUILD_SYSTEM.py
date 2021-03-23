#/usr/bin/env python


import sys, os

print("This is the goAnalyze! build system")
print("If you are a user you probably don't need this. Use the INSTALL.py script, instead")

dirs=["gorama","goreduce","goopt","gotraj","gomd","goshape"]
OS=["linux","darwin","windows"]

#OS=["linux"]
#dirs=["gomd"]

base=os.getcwd()  #all of it are belong to us


for j in OS:
	os.system("mkdir bin_%s"%j.title())

for i in dirs:
	os.chdir(i)
	for j in OS:
		print("GOOS=%s GOARCH=386 go build"%j)
		binname=i
		if j=="windows":
			binname=i+".exe"
		os.system("GOOS=%s GOARCH=386 go build"%j)
		os.system("mv %s ../bin_%s/"%(binname,j.title()))
	os.chdir(base)
	
