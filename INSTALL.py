#/usr/bin/env python

from os.path import expanduser
import sys, os

print("This is the goAnalyze! install system")

if len(sys.argv)<1 or not sys.argv[1] in ["linux","darwin"]:
	print("USE: python INSTALL.py OS_NAME\n OS_NAME must be either 'linux' or 'darwin'")


OS=sys.argv[1]

base=os.getcwd()  #all of it are belong to us

#You need to have bash as your shell, which is, apparently,
#not the default in MacOS. Sorry about that, I don't know how zsh 
#configuration works.
if OS in ["linux", "darwin"]:
	home= expanduser("~")
	os.chdir(home)
	#shell settings first. I assume Bash. If you use zsh or something, 
	# you'll have to edit this, or just set things by hand, you hipster :-D
	bash=open(".bashrc","a")
	if os.path.isfile(".zprofile"):
	    	bash.close()
		bash=open(".zprofile","a") #not tested!
	bash.write("\n# goAnalyze stuff\n")
	bash.write("export PATH=%s/bin_%s:$PATH\n"%(base,OS.title()))
	bash.write("export PYTHONPATH=%s:$PYTHONPATH #needed for the gochem.py library\n"%(base))
	bash.close()
	os.system("source ~/.bashrc")
	pymol=open(".pymolrc","a")
	pymol.write("load %s/goanalyze.py\n"%(base))
	pymol.write("load %s/gomd.py\n"%(base))
	pymol.close()
	os.chdir(base)


print("Installation ready. Now goAnalyze! (sorry about that :-) ")
