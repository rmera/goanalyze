#!/usr/bin/python

# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# ----------------------------------------------------------------------
# This PyMOL Plugin is Copyright (C) 2013-2021 by Raul Mera-Adasme
# 
#						All Rights Reserved
# 
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
# 
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ------------------------------


##This work is dedicated to the long life of the Ven. Khempo Phuntzok Tenzin Rinpoche.##


import array
import os
import json

import gochem
from chempy.models import Indexed
from chempy import Bond, Atom
from pymol import cmd
from subprocess import Popen, PIPE

#needed to display Ramachandran plots.
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def gocap(sel):
	if sel=="help":
		print("gocap selection")
		print("Caps the N-terminal and the C-terminal dangling bonds of a subsequence, if needed")
		print("A new object is produced that can be then fed to goopt")
		print("a 'tofix' selection with atoms recommended to be fixed in an optimization is also produced for convenience")
		print("If gocap is used on several different selections, tofix will accumulate the atoms to be fixed on every selection used")
		print("So it can be fed directly to goopt.")
		return
	selo=cmd.get_model(sel,state=-1)
	cmd.delete("tmV")
	if selo.atom[0].chain!="":
		cmd.select("tmpV","chain "+selo.chain)
	else:
		cmd.select("tmpV","polymer.protein")
	last=cmd.get_model("tmpV",-1).atom[-1].resi_number
	cmd.delete("tmpV")
	firstresi=selo.atom[0].resi_number
	lastresi=selo.atom[-1].resi_number
	pattern=sel
	if firstresi!=1:
		pattern="(resi %d and name C) or "%(firstresi-1)+pattern
	if lastresi!=last:
		pattern=pattern+" or (resi %d and name N)"%(lastresi+1)
	cmd.select("tmpV",pattern)
	mod=cmd.get_model("tmpV",-1)
	if firstresi!=1:
		mod.atom[0].resi=firstresi
		mod.atom[0].name="HT"
		mod.atom[0].symbol="H"
	if lastresi!=last:
		mod.atom[-1].resi=lastresi
		mod.atom[-1].name="HT"
		mod.atom[-1].symbol="H"
	cmd.delete("tmpV")
	allnames=cmd.get_names("all")
	name=sel+"__cap" #trying to choose a a suffix that is unlikely to be in another selection.
	cont=0
	while name in allnames:
		name=name+str(cont)
		cont+=1
	cmd.load_model(mod,name)
	if "tofix" in allnames:
		print("tofix or (%s and resi %d and name N) or (%s and resi %d and name C)"%(name,firstresi,name,lastresi))
		cmd.select("tofix", "tofix or (%s and resi %d and name N) or (%s and resi %d and name C)"%(name,firstresi,name,lastresi))
	else:
	#	print("(%s_c and resi %d and name N) or (%s_c and resi %d and name C)"%(sel,firstresi,sel,lastresi))
		cmd.select("tofix","(%s and resi %d and name N) or (%s and resi %d and name C)"%(name,firstresi,name,lastresi))
	

		


def goopt(sel, quality="high",fixed=None,charge=0,multi=1,solvent="water",dryrun=False):
	if sel=="help":
		print("goopt selection quality=high, fixed=None, charge=0, multi=1,solvent=''\n")
		print("Optimizes the selection with xtb (must be installed and in the PATH)")
		print("The optimized structure is loaded into PyMOL")
		print("Selections can't have dangling bonds! if you want to optimize a sub-sequence in a protein")
		print("Use the gocap command to cap the dangling bonds first")
		print("fixed is a selections with the atoms to be fixed/constrained during the optimization")
		print("gocap provides a 'tofix' selection you can use here")
		print("Only one selection is allowed, so put together all your selections into one")
		print("Remember to get the charge right!")
		print("3 qualities for the optimization are supported. 'high' (but slower) 'medium' and 'slow'")
		print("3 possible solvents are supported. water, protein (actually modeled as dicloromethane, e=5, and no solvent (default)")
		print("goopt always takes the first state for each selection. If you wish a different state, copy it as a new object (which gocap does)")
		return
	qual={"vhigh":"gfn2","high":"gfn0","medium":"gfnff"}
	diel={"water":80,"protein":5,"":-1}
	proc = Popen("goopt", shell=True, stdin=PIPE, stdout=PIPE)
	fix=[]
	#If we get a selection with atoms to be fixed
	#we can't just pick those indexes and give them to Go.
	#we need to find the indexes in the large selection that have the same
	#id as the atoms in the fixed selection, and give _those_ indexes to Go.
	if fixed:
		mod=cmd.get_model(sel)
		tofix=cmd.get_model(fixed)
		fixids=[]
		for i in tofix.atom:
			fixids.append(i.id)
		for i,v in enumerate(mod.atom):
			if v.id in fixids:
				fix.append(i)
	opt=gochem.Options()
	method=qual[quality]
	opt.AddStringOptions([method])
	opt.AddBoolOptions([dryrun])
	opt.AddIntOptions([int(charge),int(multi)])
	opt.AddIntOptions(fix)
	opt.AddFloatOptions([float(diel[solvent])])
	stat={"first":-1,"last":0,"skip":1}
	#Ojects created from Chempy and objects read by PyMOL are somehow treated differently by the program.
	#For a PDB its OK to get the current "-1" state, but for the chempy-ones, it causes errors later. So,
	#I just put this completely ad-hoc fix. IF gocap created the object (thus, it's a Chempy object) we 
	#ask for the first state, not the current. These objects have one state, so it's the same, but it prevents an
	#error later on.
	if sel.endswith("__cap"):
		stat={"first":0,"last":0,"skip":1}
	gochem.SendgoChem(proc,selnames=[sel],states=[stat],options=opt)
	info = gochem.get_info(proc)
	mod=gochem.get_model(proc,info,0)
	gochem.LoadNew(mod,sel+"_opt")



def gotraj(sel,filename="unlikelytobearealname.exe",skip=1,extension=None):
	if sel=="help":
		print("gotraj selection traj_filename skip\n")
		print("The trajectory file will be loaded on the selection, in a new object")
		print("A frame will be read every skip frames (default=1, all frames are read)")
		return
	ed = {"crd":"amber","trj":"amber","xyz":"xyz","xtc":"xtc","dcd":"dcd","pdb":"pdb","xtc":"xtc"}
	proc = Popen("gotraj", shell=True, stdin=PIPE, stdout=PIPE)
	if filename=="unlikelytobearealname.exe":
		raise "usage: gotraj selection traj_file_name\n (you must give the name of the trajectory to be read)"
	if not extension:
		ext=filename.split(".")[-1]
		extension=ed[ext] #if this doesn't work you get a well deserved exception
	if extension=="xtc":
		print("this extension is probably not supported!") #you'll get a crash later 
	opt=gochem.Options()
	opt.AddStringOptions([filename,extension])
	opt.AddIntOptions([int(skip)])
	gochem.SendgoChem(proc,selnames=[sel],options=opt,sort=True)
	info = gochem.get_info(proc)
	mod,states=gochem.get_model(proc,info,0)
	gochem.LoadNew(mod,sel+"_Traj")
	while states.GetState(mod):
		cmd.load_model(mod,sel+"_Traj") #we don't want to re-calculate bonds, so we call this function instead.




#Uses goChem and Reduce to protonate a selection, and display the protonated version
#as a new model. Both Reduce and the goChem part of the plugin need to be in the PATH.
def goreduce(sel):
	if sel=="help":
		print("goreduce selectio\n")
		print("The selection will be used in the _current_ state")
		print("The selection (must be a protein!) will be protonated with") 
		print("Reduce, which must be installed and in the PATH)")
		print("The protonated structure will be shown in PyMOL")
		return
	proc = Popen("goreduce", shell=True, stdin=PIPE, stdout=PIPE)
	gochem.SendgoChem(proc,selnames=[sel])
	info = gochem.get_info(proc)
	mod=gochem.get_model(proc,info,0)
	gochem.LoadNew(mod,sel+"_H")




#this reads 2 structures, calculates RMSD for each atom in backbone and assigns that value to b-factor of the atom. It sets b-factors for all other atoms to 0
def gorama(*args,**kwargs):
	if args[0]=="help":
		print("gorama selection1, selection2, ..., selectionN")
		print("all selections will be used in their _current_ state")
		print("Each selection will be colored differently in the plot, from blue (first) to red (last)")
		print("If one selection is given the colors will go from blue to red with increasing residue number.")
		return
	seles=args
	proc = Popen("gorama", shell=True, stdin=PIPE)
	opt=gochem.Options()
	number=0
	name="_".join(args)+"_Rama_%d.png"%number
	while os.path.isfile(name):
		number+=1
		name="_".join(args)+"_Rama_%d.png"%number
	opt.AddStringOptions([name.replace(".png","")])
	gochem.SendgoChem(proc,seles,options=opt)
	if  proc.wait() != 0:
		print("There were some errors")	
	print("The points are colored by residue")
	print("First residues in the (sub)-sequence are in blue")
	print("Last residues, in red")
	img = mpimg.imread("./"+name)
	fig = plt.imshow(img)		
	fig.set_cmap('hot')
	fig.axes.get_xaxis().set_visible(False)
	fig.axes.get_yaxis().set_visible(False)
	plt.show()



def goshape(*args,**kwargs):
	if args[0]=="help":
		print("goshape selection1, selection2, ..., selectionN")
		print("Obtains shape indicators (percentage of linearity/elongation and percentage of")
		print("planarity, both based on the eigenvalues of the moment of inertia tensor, for each")
		print("selection given")
		print("all selections will be used in their _current_ state")
		print("The indicators are based on those proposed by")
		print("Taylor et al., .(1983), J Mol Graph, 1, 3")
		return
	seles=args
	proc = Popen("goshape", shell=True, stdin=PIPE,stdout=PIPE)
	gochem.SendgoChem(proc,seles)
	info = gochem.get_info(proc)
	lincirc=info["FloatInfo"]
	for i,v in enumerate(seles):
		print("%s: Elongation %4.1f %% Planarity %4.1f %%"%(v,lincirc[i][0],lincirc[i][1]))



def gohelp(*args,**kwargs):
		print("goModel! Available commands and their basic use:")
		print("gocap selection")
		print("goopt selection quality=high, fixed=None, charge=0, multi=1,solvent=''")
		print("gorama selection1 selection2 ... selectionN")
		print("goreduce selection")
		print("goshape selection1, selection2, ..., selectionN")
		print("gotraj selection traj_filename skip")
		print("gomd task, sel1, selection2 ... (trajectory analysis)")
		print("for more info on each command use:")
		print("command help")
		return

def gomodel(*args,**kwargs):
	print("goModel! is a set of tools for molecular modelling on PyMOL, written in Go and Python")
	gohelp()


cmd.extend("gorama", gorama)
cmd.extend("goreduce", goreduce)
cmd.extend("gotraj", gotraj)
cmd.extend("goopt", goopt)
cmd.extend("gocap", gocap)
cmd.extend("goshape",goshape)

#Die hilfe
cmd.extend("gohelp", gohelp)
cmd.extend("gomodel", gomodel)


