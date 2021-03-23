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
import math
import gochem
from chempy.models import Indexed
from chempy import Bond, Atom
from pymol import cmd
from subprocess import Popen, PIPE

#needed to display Ramachandran plots.
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg


#takes hue (0-360), v and s (0-1), returns r,g,b (0-255)
def iHVS2RGB(h, v, s):
	i =0.0
	f = 0.0
	p = 0.0
	q = 0.0
	t = 0.0
	r =0.0
	g=0.0
	b=0.0
	maxcolor = 255.0
	conversion = maxcolor * v
	if s == 0.0:
		return int(conversion), int(conversion), int(conversion)
	h = h / 60
	i = math.floor(h)
	f = h - i
	p = v * (1 - s)
	q = v * (1 - s*f)
	t = v * (1 - s*(1-f))
	if int(i)==0:
		r = v
		g = t
		b = p
	elif int(i)==1: 
		r = q
		g = v
		b = p
	elif int(i)==2:
		r = p
		g = v
		b = t
	elif int(i)==3:
		r = p
		g = q
		b = v
	elif int(i)==4:
		r = t
		g = p
		b = v
	else: #case 5
		r = v
		g = p
		b = q
	r = r * conversion
	g = g * conversion
	b = b * conversion
	return int(r), int(g), int(b)


#colors assigns an color to the number key  depending
#on its place in a linear interpolation from 0 to steps
#the colors go by decreasing hue with step, from 0 to 240
#we stop at a hue of 240 to avoid violet colors that can be
#confusing
#meaning that keys closer to 0 will be blue, while keys close
#to steps will be red.
def colors(key, steps):
	MAXHUE  = 240.0
	norm  = MAXHUE / float(steps)
	hp = float((float(key) * norm))
	h = MAXHUE - hp #we invert the order
	s = 1.0
	v = 1.0
	r, g, b = iHVS2RGB(h, v, s)
	return r/255, g/255, b/255 #apparently this is what matplotlib wants



#this reads 2 structures, calculates RMSD for each atom in backbone and assigns that value to b-factor of the atom. It sets b-factors for all other atoms to 0
def gomd(*args,**kwargs):
	if args[0]=="help":
		print("gomd: gomodel! MD analysis module")
		print("gomd taks selection1, selection2, ..., selectionN")
		print("available tasks are:")
		print("gomd Distances, sel1a, sel1b, sel2a, sel2b ... \n Plots the distances between pairs of selections along the simulation\n")
		print("gomd Elongation, sel1, sel2, sel3, ... \n Plots the % of elongation for the selections given\n")
		print("gomd Ramachandran, sel1\n  draws the Ramachandran plot for a selection along the trajectory (each frame is colores differently. Takes one selection\n" )
		print("gomd RMSD, ref1, sel1, ref2, sel2 ...\n plot the rmsd of a selection along a simulations, against a fixed reference\n")
		print("gomd RMSF, sel1\n plots the per-atom RMSF for all atoms in the selection. Takes one selection only\n" )
		print("gomd Super  reference_selection, test_selection, move_selection\n Obtains transformations to superimpose each frame of test_selection to reference_selection and applies the transformations to the different frames of move_selection")
		return
	task=args[0].upper()
	if task=="DISTANCE":
		task=="DISTANCES"
	seles=args[1:]
	print(task,"*",seles)
	proc = Popen("gomd", shell=True, stdin=PIPE, stdout=PIPE)
	opt=gochem.Options()
	number=0
	opt.AddStringOptions([task])
	skip=1
	if "skip" in kwargs.keys():
		skip=int(kwargs["skip"])
	opt.AddIntOptions([1]) #the skip is handled in Python. Go just getst the already trimmed trajectory, so it never has to skip
	statesdic=[]
	for i in  seles:
		statesdic.append({"first":0,"last":-1,"skip":skip})
	gochem.SendgoChem(proc,seles,states=statesdic,options=opt)
	info = gochem.get_info(proc)
	if task=="SUPER":
		mod,states=gochem.get_model(proc,info,0)
		gochem.LoadNew(mod,seles[2]+"_sup")
		while states.GetState(mod):
			cmd.load_model(mod,seles[2]+"_sup") #we don't want to re-calculate bonds, so we call this function instead.
		return
	results=info["FloatInfo"]
	if task=="RAMACHANDRAN":
		ramas=[]
		for i in results:
			ramas.append({"phi":[],"psi":[]})
			for j,v in enumerate(i):
				if j%2==0: 
					ramas[-1]["phi"].append(v)

				else: 
					ramas[-1]["psi"].append(v)
		fig, ax = plt.subplots()
		cmap=plt.get_cmap("jet")#"Spectral")
		plt.xlabel('Phi (deg)')
		plt.ylabel('Psi (deg)')
		axes = plt.gca()
		axes.set_xlim([-180,180])
		axes.set_ylim([-180,180])
		for i,v in enumerate(ramas):
			phi=np.array(v["phi"])
			psi=np.array(v["psi"])
			r,g,b=colors(i,len(ramas))
			alpha=1.0
			plt.scatter(phi,psi,color=(r,g,b,alpha),marker="+")
		plt.show()
		return
	if task=="RMSF":	
		fig, ax = plt.subplots()
		r=list(range(1, len(results[0])+1))
		r=np.array(r)
		rmsf=np.array(results[0])
		print(r,rmsf)
		ax.plot(r, rmsf,".-")
		ax.set(xlabel='Residue/Atom ', ylabel='A',title=task)
		plt.show()
		return
	#For distances I'll have half as many results as I have selections, since a distance is obtained for 2 selections
	#so I re-structure the seles list to contain half as many elements, each with the list of the 2 selections measured.
	if task=="DISTANCES":
		s=[]
		for i,v in enumerate(seles):
			if i%2!=0:
				continue
			s.append(v+"-"+seles[i+1])
		seles=s
	#A similar thing happens with RMSD, the first selection is always the reference, so we take it away from the names.
	if task=="RMSD":
		s=[]
		for i,v in enumerate(seles):
			if i%2==0:
				continue
			s.append(v)
		seles=s
	glyphs=["b-","r-","g-","k-","c-","m-","k--","b^-","ro-","g.-","c:"]
	runav=int(len(results[0])/50)
	dorunav=True
	runavlegend="(%d -Running average)"%(runav)
	unit="A"
	if task in ["PLANARITY","ELONGATION"]:
	    unit="%"
	if "runav" in kwargs.keys() and kwargs["runav"].upper()=="FALSE":
		dorunav=False
		runavlegend=""
	fig, ax = plt.subplots()
	ax.set(xlabel="Frame %s"%(runavlegend), ylabel=task+' '+unit,title=task.title())
	for i,y in enumerate(results):
		if i>len(glyphs)-1:
			"No more glyphs to plot"
			continue
		if dorunav:
			y=np.convolve(y, np.ones((runav,))/runav, mode='valid')
		xs=list(range(1, len(y)+1))
		ax.plot(xs,y,glyphs[i],label=seles[i])
	ax.legend(loc='upper right')
	plt.show()





	
cmd.extend("gomd", gomd)


