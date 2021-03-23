/*
To the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche.

goMD a little tool for the analysis of MD trajectories.


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
	"bufio"
	"fmt"
	"log"
	"math"
	"os"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/chemjson"
	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/mat"
	//	"sort"
)

////use:  program [-skip=number -begin=number2] Task pdbname xtcname skip sel1 sel2 .... selN. Some tasks may require that N is odd that n is even.
func main() {

	//This is the part that collects all the data from PyMOL, with all  the proper error checking.
	stdin := bufio.NewReader(os.Stdin)
	options, errj := chemjson.DecodeOptions(stdin)
	if errj != nil {
		fmt.Fprint(os.Stderr, errj.Marshal())
		log.Fatal(errj)
	}
	mols := make([]*chem.Molecule, 0, len(options.SelNames))
	for k := range options.SelNames {
		mol, coords, errj := chemjson.DecodeMolecule(stdin, options.AtomsPerSel[k], options.StatesPerSel[k])
		if errj != nil {
			fmt.Fprint(os.Stderr, errj.Marshal())
			log.Fatal(errj)
		}
		mol.FillMasses()
		tmpmol, err := chem.NewMolecule(coords, mol, nil)
		if err != nil {
			log.Fatal(err)
		}
		mols = append(mols, tmpmol)
	}
	trajs := make([]chem.Traj, 0, len(mols))
	for _, v := range mols {
		trajs = append(trajs, v)
	}
	task := options.StringOptions[0][0]
	skip := options.IntOptions[0][0]
	var f func([]*v3.Matrix) ([]float64, error)
	var g func() []float64
	var h func() []*v3.Matrix
	var results [][]float64
	info := new(chemjson.Info)
	switch task {
	case "DISTANCES":
		f = Distance(mols)
		results = mdan(trajs, f, skip)
		results = transposeResults(results)
	case "RMSD":
		f = RMSD(mols)
		trajs2 := make([]chem.Traj, 0, len(mols)/2)
		//the trajs elements 0, 2, 4 and so on are, for the RMSD calculations, the references, they are not to be included as trajectories.
		for i, v := range trajs {
			if i%2 == 0 {
				continue
			}
			trajs2 = append(trajs2, v)
		}
		results = mdan(trajs2, f, skip)
		results = transposeResults(results)
	case "RAMACHANDRAN":
		f, _ := Ramachandran(mols[0], 0)
		results = mdan(trajs, f, skip)

	case "RMSF":
		f, g = RMSF(mols)
		mdan(trajs, f, skip)
		results = append(results, make([]float64, 0, 0))
		results[0] = g()
	//I should get 3 mols here: test templa and move. Move is the one we will use to send back
	case "SUPER":
		f, h = Super(mols[0])
		mdan(trajs[1:], f, skip) //the first one is the traj
		newcoords := h()
		info.Molecules = 1
		info.FramesPerMolecule = []int{len(newcoords)}
		info.AtomsPerMolecule = []int{mols[2].Len()}
		if err2 := info.Send(os.Stdout); err2 != nil {
			log.Fatal(err2)
		}
		if err2 := chemjson.SendMolecule(mols[2], newcoords, nil, nil, os.Stdout); err2 != nil {
			log.Fatal(err2)
		}
		return

	case "PLANARITY":
		f = Shape(mols, true)
		results = mdan(trajs, f, skip)
		results = transposeResults(results)
	case "ELONGATION":
		f = Shape(mols, false)
		results = mdan(trajs, f, skip)
		results = transposeResults(results)
	default:
		log.Fatal("Task parameter invalid or not present: " + task)
	}
	info.Molecules = 0
	info.FramesPerMolecule = []int{0}
	info.AtomsPerMolecule = []int{0}
	info.FloatInfo = results
	if err2 := info.Send(os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}
}

/********General helper functions************/

//takes slice of N slices of M floats, and returns a slice of M slices of N floats each.
func transposeResults(inp [][]float64) [][]float64 {
	ret := make([][]float64, len(inp[0])) //all subslices in inp should be of the same size!
	for i := range ret {
		ret[i] = make([]float64, 0, len(inp))
	}
	for _, v := range inp {
		for i := range v {
			ret[i] = append(ret[i], v[i])
		}
	}
	return ret
}

//CenterOfMass returns the center of mass the atoms represented by the coordinates in geometry
//and the masses in mass, and an error. If mass is nil, it calculates the geometric center
func com(geometry *v3.Matrix, mass *mat.Dense) (*v3.Matrix, error) {
	if geometry == nil {
		return nil, fmt.Errorf("nil matrix to get the center of mass")
	}
	gr, _ := geometry.Dims()
	if mass == nil { //just obtain the geometric center
		tmp := ones(gr)
		mass = mat.NewDense(gr, 1, tmp) //gnOnes(gr, 1)
	}
	tmp2 := ones(gr)
	gnOnesvector := mat.NewDense(1, gr, tmp2)
	ref := v3.Zeros(gr)
	ref.ScaleByCol(geometry, mass)
	//	fmt.Println("ref", ref) ///////////////////
	ref2 := v3.Zeros(1)
	ref2.Mul(gnOnesvector, ref)
	ref2.Scale(1.0/mat.Sum(mass), ref2)
	return ref2, nil
}

//returns a flat64 slice of the size requested filed with ones
func ones(size int) []float64 {
	slice := make([]float64, size, size)
	for k, _ := range slice {
		slice[k] = 1.0
	}
	return slice
}

/*****************RMSD family ***********/
//RMSD returns a function that will calculate the RMSD of as many selections as requested from a given set of coordinates against the coordinates
//in the mol object.
func RMSD(mols []*chem.Molecule) func([]*v3.Matrix) ([]float64, error) { //	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(mols)
	if argslen < 1 {
		panic("RMSD: Not enough arguments, need at least one!")
	}
	refs := make([]*v3.Matrix, 0, len(mols))
	temps := make([]*v3.Matrix, 0, len(mols))
	for i, v := range mols {
		if i%2 == 0 {
			tr := v3.Zeros(v.Len())
			tr.Copy(v.Coords[0]) //the refs are already correctly filled
			refs = append(refs, tr)
			continue
		}
		ttemp := v3.Zeros(v.Len())
		temps = append(temps, ttemp)

	}
	//we expect the []coords to be without the reference structures, i.e. to be half as long as the []mols that was passed to RMSD
	ret := func(coords []*v3.Matrix) ([]float64, error) {
		RMSDs := make([]float64, 0, len(mols))
		for i, v := range coords {
			rmsd, err := memRMSD(v, refs[i], temps[i])
			if err != nil {
				return nil, err
			}
			RMSDs = append(RMSDs, rmsd)
		}
		return RMSDs, nil
	}
	return ret
}

//smemRMSD calculates the RMSD between test and template, considering only the atoms
//present in the testlst and templalst for each object, respectively.
//It does not superimpose the objects.
//To save memory, it asks for the temporary matrix it needs to be supplied:
//tmp must be Nx3 where N is the number
//of elements in testlst and templalst
func memRMSD(ctest, ctempla, tmp *v3.Matrix) (float64, error) {
	if ctest.NVecs() != ctempla.NVecs() || tmp.NVecs() != ctest.NVecs() {
		return -1, fmt.Errorf("memRMSD: Ill formed matrices for memRMSD calculation")
	}
	tmp.Sub(ctest, ctempla)
	rmsd := tmp.Norm(2)
	return rmsd / math.Sqrt(float64(ctest.NVecs())), nil

}

/*******RMSF functions Family***************/

//RMSF returns a function that will process several frames of a trajectory for a structure, obtain the RMSF for
//all atoms in the structure, with the molecule passed to RMSF as a reference.
func RMSF(mols []*chem.Molecule) (func(coord []*v3.Matrix) ([]float64, error), func() []float64) {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(mols)
	if argslen < 1 {
		panic("RMSF: Not enough arguments, need at least one!")
	}

	frames := 0.0

	accc := v3.Zeros(mols[0].Len()) //Accumulate the coordiantes
	d := v3.Zeros(mols[0].Len())    //Accumulate the coordiantes
	acc := v3.Zeros(mols[0].Len())  //Accumulate the coordiantes
	acc2 := v3.Zeros(mols[0].Len()) //accumulate the squares of the coordinates
	//only the first coordinate set will be considered.
	ret := func(coords []*v3.Matrix) ([]float64, error) {
		fmt.Fprintln(os.Stderr, "llaman a esta wea?", coords) ///////////
		frames++
		c := coords[0]
		acc.Add(acc, c)
		d.Dense.MulElem(c, c)
		acc2.Add(acc2, d)
		//	fmt.Fprintln(os.Stderr, "weaitas0", acc, acc2, c, d) ///////////
		return []float64{0.0}, nil //Dummy output
	}
	final := func() []float64 {
		variance := v3.Zeros(mols[0].Len())
		//		fmt.Fprintln(os.Stderr, "weaitas", acc, acc2) ///////////
		acc2.Scale(1/frames, acc2)
		acc.Scale(1/frames, acc)
		//		fmt.Fprintln(os.Stderr, "weaitas2", acc, acc2) ///////////
		//	accc.CloneFrom(acc)
		accc.Dense.MulElem(acc, acc)
		variance.Sub(acc2, accc)
		res := make([]float64, 0, mols[0].Len())
		var v *v3.Matrix
		for i := 0; i < variance.NVecs(); i++ {
			v = variance.VecView(i)
			for j := 0; j < 3; j++ {
				v.Set(0, j, math.Sqrt(v.At(0, 1)))
			}
			res = append(res, v.Norm(2))

		}
		//		fmt.Fprintln(os.Stderr, "RMSFFFF", res) ///////////
		return res

	}
	return ret, final
}

//********The Distance family functions**********//

func Distance(mols []*chem.Molecule) func([]*v3.Matrix) ([]float64, error) {
	//	fmt.Println("Use: MDan distance sel1 sel2...")
	argslen := len(mols)
	if (argslen)%2 != 0 {
		panic("Distance: The number of selections must be even!")
	}
	masses := make([]*mat.Dense, 0, len(mols))
	for _, v := range mols {
		v.FillMasses()
		m, err := v.Masses()
		if err != nil {
			log.Printf("Distances: %s Will use the geometric center for this molecule", err.Error()) //one of the few cases where I actually handle an error!
			masses = append(masses, nil)
		} else {
			masses = append(masses, mat.NewDense(v.Len(), 1, m))
		}
	}
	var vec1, vec2 *v3.Matrix
	distvec := v3.Zeros(1) //the distance vector
	var err error
	ret := func(coords []*v3.Matrix) ([]float64, error) {
		distances := make([]float64, 0, len(coords)/2)
		for i := 0; i < len(coords); i = i + 2 { //we process them by pairs
			if len(coords) != argslen {
				return nil, fmt.Errorf("Wrong number of coordinate sets! is: %d should be %d", len(coords), argslen)
			}

			//	println("get to the chooopaaaaa!")

			vec1, err = com(coords[i], masses[i])
			if err != nil {
				return nil, fmt.Errorf("Distance: coudn't get the COM for selection %d: %s", i, err.Error())
			}
			vec2, err = com(coords[i+1], masses[i+1])
			if err != nil {
				return nil, fmt.Errorf("Distance: coudn't get the COM for selection %d : %s", i+1, err.Error())
			}
			distvec.Sub(vec2, vec1)
			distances = append(distances, distvec.Norm(2))
		}
		return distances, nil
	}
	return ret
}

//If you use skip_first, you have to make sure the slice of coords has enough slots!
func get_trajs(trajs []chem.Traj, coords []*v3.Matrix, read, skip_first bool) error {
	var err error
	skip := 0
	if skip_first {
		skip = 1
	}
	for i, v := range trajs {
		if read {
			err = v.Next(coords[i+skip])
		} else {
			err = v.Next(nil)
		}
		if err != nil {
			return err
		}
	}
	return nil

}

//mdan takes a topology, a trajectory object and a function that must take a set of coordinates
//and a topology and returns a slice of floats. It applies the function to each snapshot of the trajectory.
//It then, for each snapshot, prints a line with the traj number as first field and the numbers in the returned
//slice as second to N fields, with the fields separated by spaces.
//the passed function should be a closure with everything necessary to obtain the desired data from each frame
//of the trajectory.
func mdan(trajs []chem.Traj, f func([]*v3.Matrix) ([]float64, error), skip int) [][]float64 {
	read := true
	data := make([][]float64, 0, 500) //I'm just throwing a wild guess with the 500
	//this serves no purpose now, but one could use it to allow the caller giving an additional matrix, say a reference to add to the trajs
	skip_first := false
	coordset := make([]*v3.Matrix, 0, len(trajs)+1)
	for _, v := range trajs {
		coordset = append(coordset, v3.Zeros(v.Len()))
	}
	for i := 0; ; i++ { //infinite loop, we only break out of it by using "break"  //modified for profiling
		if i%skip != 0 {
			read = false
		} else {
			read = true
		}
		err := get_trajs(trajs, coordset, read, skip_first)
		if err != nil {
			_, ok := err.(chem.LastFrameError)
			if ok || err.Error() == "EOF" {
				break //We processed all frames and are ready, not a real error.

			} else {
				panic(err.Error())
			}
		}
		if !read { //not so nice check for this twice
			continue
		}
		//The important part
		numbers, err := f(coordset)
		data = append(data, numbers)
		if err != nil {
			panic(fmt.Sprintf("Error: %s in frame %d", err.Error(), i))
		}
	}
	return data
}
