#!/usr/bin/python

# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# ----------------------------------------------------------------------
# This PyMOL Plugin is Copyright (C) 2013 by Raul Mera-Adasme
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
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR-
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ------------------------------


##This work is dedicated to the long life of the Ven. Khempo Phuntzok Tenzin Rinpoche.##

	
from chempy.models import Indexed
from chempy import Bond, Atom
from pymol import cmd
import json


BONDNOHMAX=1.9
BONDHMAX=1.45
SSBONDMAX=2.3



## The options have to be given by completing and json-marshalling the follwing dictionary




#states is a list of the states wanted, starting from zero
#each element in "states" has to mark the first, last and skip for that selection.
#if "last" is -1 we assume you want up to the last state:
#you can choose if wait for the Go program to end, or (default) not.
def SendgoChem(proc,selnames=["sele"],states=[],options=None,wait=False,sort=False):
	if not options:
		options=Options()
	sels=[]
	options._o["SelNames"]=selnames
	for i,v in enumerate(selnames):
		if not states or len(states)<len(selnames):
			states.append({"first":-1,"last":0,"skip":1})  #if first==-1 only the current state will be read.
		sels.append([])
		cont=-1
		while True:
			if states[i]["first"]>cont:
				cont+=1
				continue
			if (cont)%states[i]["skip"]!=0:
				cont+=1
				continue
			if cont>states[i]["last"] and states[i]["last"]>=0:
				break
			sels[-1].append(cmd.get_model(v,cont))
			#if you are going to read a trajectory, you have to do this, for your selection atom's and the trajs' atoms to match, because pymol,
			# for some reason, re-orders the atoms in the object
			if sort:
				sels[-1][-1].atom=sorted(sels[-1][-1].atom,key=lambda x:x.id)
			cont+=1
			if sels[-1][-1].atom==[]:
				del sels[-1][-1]
				break
			if states[i]["first"]==-1:
				break
	APS=[]
	SPS=[]
	for i in sels:
		APS.append(len(i[0].atom))
		SPS.append(len(i))
	options._o["AtomsPerSel"]=APS
	options._o["StatesPerSel"]=SPS
	#print(options ._o) #######################
	#We have all the data now. We just have to send it.
	#The options first, as this is what the Go side will be expecting
	options.Transmit(proc)
	#Now the molecule data.
	#I freely admit that these nested loops are horrible.
	for i,v in enumerate(sels):
		for j,w in enumerate(v):
			for k in w.atom:
				atom,coords=Atom2gcRef(k)
				if j==0: #i.e. this is the first (or only) state in the selection.
					proc.stdin.write((atom+"\n").encode('utf-8'))
					proc.stdin.write((coords+"\n").encode('utf-8'))
				else:
					proc.stdin.write((coords+"\n").encode('utf-8')) #This is a silly way to do it, I'm trying to debug this function
	proc.stdin.close()
	if  wait and proc.wait() != 0:
		print("There were some errors")



#Transforms an atom into a JSON structure that goChem can decode
def Atom2gcRef(i):
	Ref=json.dumps({"Name":i.name,"Id":i.id,"Molname":i.resn,"Symbol":i.symbol,"Molid":int(i.resi_number),"Chain":i.chain})
	coords=json.dumps({"Coords":i.coord})
	return Ref,coords
	



class Options:
	def __init__(self):
		self._o={"SelNames":[],"AtomsPerSel":[],"StatesPerSel":[],"StringOptions":[],"IntOptions":[],"BoolOptions":[],"FloatOptions":[]}
	def AddSelName(self,name):
		if isinstance(name,str):
			self._o["SelNames"].append(name)
		else:
			return False
		return True
	def AddAtomsPerSel(self,atoms):
		if isinstance(atoms,int):
			self._o["AtomsPerSel"].append(atoms)
		else:
			return False
		return True
	def AddStatesPerSel(self,states):
		if isinstance(states,int):
			self._o["AtomsPerSel"].append(states)
		else:
			return False
		return True
	def AddStringOptions(self,strings):
		for i in strings:
			if not isinstance(i,str):
				return False
		self._o["StringOptions"].append(strings)
		return True
	def AddIntOptions(self,ints):
		for i in ints:
			if not isinstance(i,int):
				return False
		self._o["IntOptions"].append(ints)
		return True
	def AddBoolOptions(self,bools):
		for i in bools:
			if not isinstance(i,bool):
				return False
		self._o["BoolOptions"].append(bools)
		return True
	def AddFloatOptions(self,floats):
		for i in floats:
			if not isinstance(i,float):
				return False
		self._o["FloatOptions"].append(floats)
		return True

	def Transmit(self,proc):
		totransmit=json.dumps(self._o)
		proc.stdin.write((totransmit+"\n").encode('utf-8'))










#returns dictionary with information returned by the go program
def get_info(proc):
	first=False
	v=proc.stdout.readline().decode("utf-8")
	info=json.loads(v)
	return info


#puts coordinates returned by the go program into the model model. Coordinates of atom with names in the 
#list name are the only one replaced, or the ones not replaced depending on whether included is True.
#proc is the process object for the go program, atom is the lenght of the model.
#The residues whose IDs are in the list resids are alwayes excluded.
def get_coords(proc,model,names,resids,included,info,number):
	atomsread=0
	atoms=info["AtomsPerMolecule"][number]
	rcoords=True
	rbfactors=False
	rss=False
	first=True
	while(True):
		v=proc.stdout.readline().decode("utf-8")
	#	print "como las w", number, v, atoms
		if "Coords" in v and not "Molname" in v and rcoords:
			coords=json.loads(v)
			if included:
				if model.atom[atomsread].name in names and not model.atom[atomsread].resi in resids:
					model.atom[atomsread].coord=coords["Coords"]
			else:
				if not model.atom[atomsread].name in names and not model.atom[atomsread].resi in resids:
					#print "YEEEY", number
					model.atom[atomsread].coord=coords["Coords"]
			atomsread=atomsread+1
			if atomsread==atoms:
				rcoords=False
				atomsread=0
				if  info["Bfactors"]:
					rbfactors=True
				if info["SS"]:
					rss=True
			continue
		#In many cases this part will not be needed
		if "Bfactors" in v and rbfactors:
			bf=json.loads(v)
			model.atom[atomsread].b=bf["Bfactors"]
			atomsread=atomsread+1
			if atomsread==atoms-1:
				atomsread=0
				rbfactors=False
				if info["SS"]:
					rss=True
			continue
		#This one should be needed only seldom
		if "SS" in v and rss:
			SS=json.loads(v)
			model.atom[atomsread].ss=SS["SS"]
			++atomsread
			if atomsread==atoms-1:
				atomsread=0
				rss=False
			continue
		break
#		print "me fui con una deuda de 500"
	return model




#replace coordinates and b-factors in vmodel (has to be a chempy.Indexed object)
#with the ones sent from the Go Program.
def GetgoChem(proc,vmodel, info,number):
	vmodel=Indexed()
	atoms=info["AtomsPerMolecule"][number]
	atomsread=0
	rcoords=True
	rbfactors=False
	rss=False
	first=True
	while(True):
		v=proc.stdout.readline().decode("utf-8")
		if "Coords" in v and not "Molname" in v and rcoords:
			coords=json.loads(v)
			vmodel.atom[atomsread].coord=coords["Coords"]
			atomsread=atomsread+1
			if atomsread==atoms:
				rcoords=False
				atomsread=0
				if  info["Bfactors"]:
					rbfactors=True
				if info["SS"]:
					rss=True
			continue
		#In many cases this part will not be needed
		if "Bfactors" in v:
			bf=json.loads(v)
			vmodel.atom[atomsread].b=bf["Bfactors"]
			atomsread=atomsread+1
			if atomsread==atoms:
				atomsread=0
				rbfactors=False
				if info["SS"]:
					rss=True
			continue
		#This one should be needed only seldomly
		if "SS" in v:
			SS=json.loads(v)
			vmodel.atom[atomsread].ss=SS["SS"]
			++atomsread
			if atomsread==atoms:
				atomsread=0
				rss=False
			continue
#		print "me fui con una deuda de 500"
		break




resncodes={"ALA":"A","GLY":"G","VAL":"V","LEU":"L","PHE":"F","MET":"M","TRP":"W","PRO":"P","ILE":"I","LYS":"K","ARG":"R","HIS":"H","GLU":"E","ASP":"E","GLN":"Q","ASN":"N","SER":"S","THR":"R","TYR":"Y","CYS":"C"}
#adds 1-letter codes to all atoms in the Indexed object. For _some_ reason pymol seems to remove them
def add_1lett(model):
	for i,v in enumerate(model.atom):
		model.atom[i].resn_code=resncodes[v.resn]




#similar to get_gochem but it doesnt take a reference model, hence, all the data is taken from the output of the go program.
#Proc is the gochem process object, atoms is the number of atoms that should be added.
#Notice that PyMOL 1.6 will NOT guess bonds for objects created with these functions, so lines, sticks, cartoons representations
#willl be NOT available. 
def GetModel(proc, info,number):
	vmodel=Indexed()
	atoms=info["AtomsPerMolecule"][number]
	states=info["FramesPerMolecule"][number]
	statecontainer=States(atoms)
	atomsread=0
	statesread=0
	ratoms=True
	rcoords=False
	rbfactors=False
	rss=False
	first=True
	while(True):
		v=proc.stdout.readline().decode("utf-8") 
		if "Molname" in v and ratoms:
			ad=json.loads(v)
			at=Atom()
			at.name=ad["Name"]
			at.symbol=ad["Symbol"]
			at.chain=ad["Chain"]
			at.id=ad["ID"]
			at.resi_number=ad["MolID"]
			at.resn=ad["Molname"]
			at.hetatm=1
			if at.resn in resncodes.keys():
				at.hetatm=0
				at.resn_code=resncodes[ad["Molname"]]
			vmodel.atom.append(at)
			atomsread=atomsread+1
			if atomsread==atoms:
				ratoms=False
				rcoords=True
				atomsread=0
			continue
		if "Coords" in v and not "Molname" in v and rcoords:
			coords=json.loads(v)
			if statesread==0:
				vmodel.atom[atomsread].coord=coords["Coords"]
			else:
				statecontainer.NewCoords(coords["Coords"])
			#	print("coords added to state!!!") #############3
			atomsread=atomsread+1
			if atomsread==atoms:
				statesread+=1
				atomsread=0
				if statesread==states:
					rcoords=False
					if  info["Bfactors"]:
						rbfactors=True
					if info["SS"]:
						rss=True
				else:
					statecontainer.NewState()
			#		print("NEW STATE!") ########3

			continue
		#In many cases this part will not be needed
		if "Bfactors" in v:
			bf=json.loads(v)
			vmodel.atom[atomsread].b=bf["Bfactors"]
			atomsread=atomsread+1
			if atomsread==atoms:
				atomsread=0
				rbfactors=False
				if info["SS"]:
					rss=True
			continue
		#This one should be needed only seldom
		if "SS" in v:
			SS=json.loads(v)
			vmodel.atom[atomsread].ss=SS["SS"]
			++atomsread
			if atomsread==atoms:
				atomsread=0
				rss=False
			continue
#		print "me fui con una deuda de 500"
		break
	#i don't like this, but I don't want to break my other plugins by changing the signature of the function
	#for them. Also, most pugins will only get one state anyway.
	if statecontainer.States()<1:
		return vmodel
	else:
		return (vmodel, statecontainer)

#this class will keep sets of coordinates for the different states of a molecule
#except for the first, which will be just loaded onto the chempy object.
class States:
	def __init__(self,natoms):
		self.C=[]
		self.natoms=natoms
		self.cs=0 #current state
		self.cc=0 #current coord
		self.statetoload=0
	def NewState(self):
		self.C.append([])
	def NewCoords(self,c):
		self.C[-1].append(c)
	def Coord(self):
		if self.cc+1<self.natoms:
			self.cc+=1
			return self.C[self.cs][self.cc]
		else:
			if self.cs+1>=len(self.C):
				return None
			self.cs+=1
			self.cc=0
			return self.C[self.cs][self.cc]
	def GetState(self,model):
		if self.statetoload>=len(self.C):
			return False #we just ran out of states, not really an error.
		if len(self.C[self.statetoload])!=len(model.atom):
			#something is quite wrong here. The States object is not yet fully prepared, or the data behind it is wrong. 
			#Either way, the program should crash.
			raise AttributeError("State %d has %d atoms, while the model has %d atoms"%(self.statetoload,len(self.C[self.statetoload]),len(model.atoms)))
		for i,v in enumerate(self.C[self.statetoload]):
			model.atom[i].coord=v
		self.statetoload+=1
		return True
	def States(self):
		return len(self.C)


		

				


get_model = GetModel #compatibility with previous versions



##A bit low level stuff
import chempy.neighbor
import chempy.cpv

def bondseeker(model, nbr,excludeH=True,ceil=BONDNOHMAX):
	pairs=[]
	for i,at in enumerate(model.atom):
		if at.symbol=="H" and excludeH:
			continue
		if at.symbol!="H" and not excludeH:
			continue
		lst = nbr.get_neighbors(at.coord)
		for b in lst:
			at2 = model.atom[b]
			if at2.symbol=='H': #we can always exclude these, since we don't expect to find H-H bonds. If you want to add H2 molecules, you are not in luck, I guess
				continue
			if excludeH and [b,i] in pairs: #we can save some comparisons because if we are not excluding H, we are sure that at is an H, and at2 isn't, so the pair [at2,at1] can never be added.
				continue
			dst = chempy.cpv.distance(at.coord,at2.coord)
			if dst>ceil:
				if not (at.symbol=="S" and at2.symbol=="S" and dst<SSBONDMAX): #allow for SS bonds.
					continue
			pairs.append([i,b])
			bnd=Bond()
			bnd.index = [i,b]
			bnd.order = 1 #yeah, not going to do the proper bond order.
			model.bond.append(bnd)



#A very simple and not-very-good bond-adding function.
#Based on the original chempy proteins.add_bonds,
#which was written, I believe, by the great Warren L. DeLano.
#This is meant for proteins. I have not idea if these cutoffs are reasonable for P-O bonds for instance.
def add_bonds(model):
	#yeah, not very sophisticated. sorry
	crd = model.get_coord_list()
	#start with the heavy atoms
	nbr = chempy.neighbor.Neighbor(crd,BONDNOHMAX)
	bondseeker(model,nbr,excludeH=True,ceil=BONDNOHMAX)
	nbr = chempy.neighbor.Neighbor(crd,BONDHMAX)
	bondseeker(model,nbr,excludeH=False,ceil=BONDHMAX)



#loads a newly-built chempy Indexed object, i.e. one without bonds.
def LoadNew(model,name, template=None,PDBstr=False,state=0,discrete=0):
	if PDBstr:
		tmpname=name+"_tmp"
		cmd.load_model(model,tmpname,discrete=discrete,zoom=1)
		modR=cmd.get_pdbstr(tmpname)
		cmd.read_pdbstr(modR,name)
		cmd.delete(tmpname)
		return
	if template:
		model.bond=template.bond
	else:
		add_bonds(model)
	cmd.load_model(model,name,state=state,discrete=discrete)






