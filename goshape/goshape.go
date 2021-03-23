package main

/*
 * goRama.
 *
 * Copyright 2013 Raul Mera <devel@gochem.org>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

import (
	"bufio"
	"fmt"
	"log"
	"os"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/chemjson"
	v3 "github.com/rmera/gochem/v3"
)

//For most plugins you dont really need to report errors like this, just panicking or using log.fatal should be enough.
//Here I use JSON errors as a test. They would be useful if you want to implement some recovery behavior.

//A map for assigning mass to elements.
//Note that just common "bio-elements" are present
var symbolMass = map[string]float64{
	"H":  1.0,
	"C":  12.01,
	"O":  16.00,
	"N":  14.01,
	"P":  30.97,
	"S":  32.06,
	"Se": 78.96,
	"K":  39.1,
	"Ca": 40.08,
	"Mg": 24.30,
	"Cl": 35.45,
	"Na": 22.99,
	"Cu": 63.55,
	"Zn": 65.38,
	"Co": 58.93,
	"Fe": 55.84,
	"Mn": 54.94,
	"Si": 28.08,
	"Be": 9.012,
	"F":  18.998,
}

func main() {
	//This is the part that collects all the data from PyMOL, with all  the proper error checking.
	stdin := bufio.NewReader(os.Stdin)
	options, errj := chemjson.DecodeOptions(stdin)
	if errj != nil {
		fmt.Fprint(os.Stderr, errj.Marshal())
		log.Fatal(errj)
	}
	mols := make([]*chem.Topology, 0, len(options.SelNames))
	coordset := make([]*v3.Matrix, 0, len(options.SelNames))
	for k := range options.SelNames {
		mol, coords, errj := chemjson.DecodeMolecule(stdin, options.AtomsPerSel[k], 1)
		if errj != nil {
			fmt.Fprint(os.Stderr, errj.Marshal())
			log.Fatal(errj)
		}
		mol.FillMasses()
		mols = append(mols, mol)
		coordset = append(coordset, coords[0])
	}

	lincirc := make([][]float64, 0, len(mols))
	const epsilon float64 = 0.00001
	for i, v := range mols {
		lin, circ, err := chem.EasyShape(coordset[i], epsilon, v)
		if err != nil {
			log.Fatal(err)
		}
		lincirc = append(lincirc, []float64{lin, circ})
	}
	info := new(chemjson.Info)
	info.FloatInfo = lincirc
	fmt.Fprintln(os.Stderr, info.FloatInfo) ///////////////
	if err2 := info.Send(os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}
	fmt.Fprintln(os.Stderr, "goshape did its duty!") //////////////
}
