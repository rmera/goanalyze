package main

/*
 * goQM.
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
	"runtime"

	"github.com/rmera/gochem/chemjson"
	"github.com/rmera/gochem/qm"
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
	mol, coordarray, err := chemjson.DecodeMolecule(stdin, options.AtomsPerSel[0], 1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	coords := coordarray[0]
	dielectric := options.FloatOptions[0][0]
	charge := options.IntOptions[0][0]
	multi := options.IntOptions[0][1]
	fixed := options.IntOptions[1]
	method := options.StringOptions[0][0]
	log.Println("fixed", fixed) //////////////////
	mol.SetCharge(charge)
	mol.SetMulti(multi)

	calc := new(qm.Calc)
	calc.Job = qm.Job{Opti: true}
	if len(fixed) > 1 {
		calc.CConstraints = fixed
	}
	if dielectric > 0 {
		calc.Dielectric = dielectric
	}
	calc.Method = method
	xtb := qm.NewXTBHandle()
	xtb.SetnCPU(runtime.NumCPU())

	xtb.BuildInput(coords, mol, calc)
	fmt.Fprint(os.Stderr, options.BoolOptions)
	if options.BoolOptions[0][0] {
		return //Dry run
	}
	if err2 := xtb.Run(true); err != nil {
		log.Fatal(err2.Error())
	}
	//Now we ran the calculation, we must retrive the geometry and divide the coordinates among the original selections.
	var newcoords *v3.Matrix
	var err2 error
	newcoords, err2 = xtb.OptimizedGeometry(mol)
	if err2 != nil {
		log.Fatal(err2.Error())
	}
	//	energy, err2 := xtb.Energy()
	//	if err2 != nil {
	//		log.Fatal(err2.Error())
	//	}
	//Start transfering data back
	info := new(chemjson.Info)
	info.Molecules = 1
	info.FramesPerMolecule = []int{1}
	info.AtomsPerMolecule = []int{mol.Len()}
	//	info.Energies = []float64{energy}
	if err2 := info.Send(os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}

	if err2 := chemjson.SendMolecule(mol, []*v3.Matrix{newcoords}, nil, nil, os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}

}
