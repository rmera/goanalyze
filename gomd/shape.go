package main

import (
	"fmt"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

//Shape returns a function that will calculate shape indicators (planarity and elongation, returned in that order) on as many selections as requested
func Shape(mols []*chem.Molecule, planarity bool) func(coord []*v3.Matrix) ([]float64, error) {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(mols)
	epsilon := 0.0001
	if argslen < 1 {
		panic("Shape: Not enough arguments, need at least one!")
	}
	ret := func(coords []*v3.Matrix) ([]float64, error) {
		shapes := make([]float64, 0, len(mols))
		for i, v := range coords {
			lin, circ, err := chem.EasyShape(v, epsilon, mols[i])
			if err != nil {
				return nil, fmt.Errorf("Shape: " + err.Error())
			}
			if planarity {
				shapes = append(shapes, circ)
			} else {
				shapes = append(shapes, lin)
			}
		}
		return shapes, nil
	}
	return ret
}
