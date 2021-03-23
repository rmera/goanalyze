/*
To the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche.

goMD, a little tool for the analysis of MD trajectories.


This program makes extensive use of the goChem Computational Chemistry library.
If you use this program, we kindly ask you support it by to citing the library as:

R. Mera-Adasme, G. Savasci and J. Pesonen, "goChem, a library for computational chemistry", http://www.gochem.org.


LICENSE

Copyright (c) 2017 Raul Mera <rmera{at}usachDOTcl>


This program, including its documentation,
is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2.0 as
published by the Free Software Foundation.

This program and its documentation is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General
Public License along with this program.  If not, see
<http://www.gnu.org/licenses/>.



*/

package main

import (
	"fmt"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	//	"sort"
)

//This is not really analysis, but a took to re-superimpose your trajectory to a different subset of the reference.
//you give the test, reference and move matrices. test is a subset of move, test and reference have the same amount
//of atoms. This will, for each frame, get the rotation+translation matrices to superimpose test and ref, and apply
//them to move.
func Super(mol *chem.Molecule) (func(coord []*v3.Matrix) ([]float64, error), func() []*v3.Matrix) {
	supertraj := make([]*v3.Matrix, 0, 500)
	ref := mol.Coords[0]
	//note that only 2 trajectories should be passes to ret (i.e. to mdan)
	ret := func(coords []*v3.Matrix) ([]float64, error) {
		test := coords[0]
		move := coords[1]
		_, rotation, trans1, trans2, err1 := chem.RotatorTranslatorToSuper(test, ref)
		if err1 != nil {
			return nil, fmt.Errorf("Super: " + err1.Error())
		}
		toappend := v3.Zeros(move.NVecs())
		toappend.AddVec(move, trans1)
		toappend.Mul(toappend, rotation)
		toappend.AddVec(toappend, trans2)
		supertraj = append(supertraj, toappend)
		return []float64{0.0}, nil //meaningless
	}
	final := func() []*v3.Matrix {
		return supertraj
	}
	return ret, final
}
