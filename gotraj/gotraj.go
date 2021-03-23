package main

// Copyright Notice
// ================
//
// The PyMOL Plugin source code in this file is copyrighted, but you can
// freely use and copy it as long as you don't change or remove any of
// the copyright notices.
//
// ----------------------------------------------------------------------
// This PyMOL Plugin is Copyright (C) 2013 by Raul Mera-Adasme
//
//						All Rights Reserved
//
// Permission to use, copy, modify, distribute, and distribute modified
// versions of this software and its documentation for any purpose and
// without fee is hereby granted, provided that the above copyright
// notice appear in all copies and that both the copyright notice and
// this permission notice appear in supporting documentation, and that
// the name(s) of the author(s) not be used in advertising or publicity
// pertaining to distribution of the software without specific, written
// prior permission.
//
// THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
// INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
// NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
// CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
// USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
// OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
// PERFORMANCE OF THIS SOFTWARE.
// ------------------------------

import (
	"bufio"
	"fmt"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/amberold"
	"github.com/rmera/gochem/chemjson"
	"github.com/rmera/gochem/dcd"
	v3 "github.com/rmera/gochem/v3"

	//	"github.com/rmera/scu"
	"log"
	"os"
)

//For most plugins you dont really need to report errors like this, just panicking or using log.fatal should be enough.
//Here I use JSON errors as a test. They would be useful if you want to implement some recovery behavior.

func main() {
	//This is the part that collects all the data from PyMOL, with all  the proper error checking.
	stdin := bufio.NewReader(os.Stdin)
	options, err := chemjson.DecodeOptions(stdin)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	mol, _, err := chemjson.DecodeMolecule(stdin, options.AtomsPerSel[0], 1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	Coords := []*v3.Matrix{}
	trajname := options.StringOptions[0][0]
	format := options.StringOptions[0][1]
	skip := options.IntOptions[0][0]
	fmt.Fprintln(os.Stderr, options) ///////////////
	var traj chem.Traj
	var err1 error
	var mat *v3.Matrix

	switch format {

	case "amber":
		traj, err1 = amberold.New(trajname, mol.Len(), false)
		if err1 != nil {
			log.Fatal(err1)
		}
	case "dcd":
		traj, err1 = dcd.New(trajname)
		if err1 != nil {
			log.Fatal(err1)
		}

	case "xyz":
		traj, err1 = chem.XYZFileRead(trajname)
		if err1 != nil {
			log.Fatal(err1)
		}

	case "pdb":
		traj, err1 = chem.PDBFileRead(trajname, false)
		if err1 != nil {
			log.Fatal(err1)
		}
	}

	for i := 0; ; i++ {
		mat = nil
		if i%skip == 0 {
			mat = v3.Zeros(mol.Len())
		}
		err := traj.Next(mat)
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				break
			}
			log.Fatal(err)
		}
		if mat != nil {
			Coords = append(Coords, mat)
		}
	}

	info := new(chemjson.Info)
	info.Molecules = 1
	info.FramesPerMolecule = []int{len(Coords)}
	info.AtomsPerMolecule = []int{mol.Len()}
	if err2 := info.Send(os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}

	if err2 := chemjson.SendMolecule(mol, Coords, nil, nil, os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}
	//	fmt.Fprintln(os.Stderr, "goChem did it's duty!")

}
