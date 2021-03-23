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
//                        All Rights Reserved
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
	"github.com/rmera/gochem/chemjson"

	//	"github.com/rmera/scu"
	"log"
	"os"
	"strings"
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
	mol, coordarray, err := chemjson.DecodeMolecule(stdin, options.AtomsPerSel[0], 1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	coords := coordarray[0]
	var rep *os.File
	//The program itself
	//	fmt.Fprintln(os.Stderr, "We'll start the program now") //////////

	for {
		reportname := strings.Join([]string{"reduce_report", options.SelNames[0], "log"}, ".")
		rep, err2 := os.Create(reportname)
		if err2 != nil {
			break
		}
		defer rep.Close()
		break
	}
	var build int
	if len(options.IntOptions) == 0 || len(options.IntOptions[0]) == 0 {
		build = 2
	} else {
		build = options.IntOptions[0][0]
	}
	newmol, err2 := chem.Reduce(mol, coords, build, rep, "reduce")
	//	fmt.Fprintln(os.Stderr, "Do we have a molecule?", newmol.Len()) /////////////////////

	if err2 != nil {
		//		if err,ok :=err2.(chem.Error); ok{
		//			fmt.Println(err.Decorate(""))
		//		}
		//		fmt.Println(err2)//////////
		fmt.Fprint(os.Stderr, chemjson.NewError("process", "chem.Reduce", err2)) //Not always fatal.
		if newmol == nil {                                                       // !strings.Contains(err2.Error(),"invalid argument"){
			log.Fatal(err2)
		}
	}
	//Start transfering data back
	//	fmt.Fprintln(os.Stderr, "We start transmitting data") /////////////

	info := new(chemjson.Info)
	info.Molecules = 1
	info.FramesPerMolecule = []int{1}
	info.AtomsPerMolecule = []int{newmol.Len()}
	if err2 := info.Send(os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}
	//	fmt.Fprint(os.Stdout,mar)
	//	fmt.Fprint(os.Stdout,"\n")
	//	fmt.Fprintln(os.Stderr, "Sent the info, now the molecule") ///////////////

	if err2 := chemjson.SendMolecule(newmol, newmol.Coords, nil, nil, os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}
	//	fmt.Fprintln(os.Stderr, "goChem did it's duty!")

}
