package main

import (
	"fmt"
	"log"
	"math"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

//Ramachandran returns a Ramachandran-plotting function, which does the work of obtaining the angles for each
//frame in a trajectory, and a "final" function that puts everything together at the end, returning a [][][]float64 that contains the 2 angles
//for all residues, for all frames, and a [][]float64 that contains the r,g,b values to color each frame differently.
func Ramachandran(mol *chem.Molecule, frames int) (func(coord []*v3.Matrix) ([]float64, error), func() [][]float64) {
	RGBs := make([][]float64, 0, frames)
	ramas := make([][]float64, 0, frames)
	ramasets := make([]chem.RamaSet, mol.Len()/10)
	ramasets, err2 := chem.RamaList(mol, "", []int{0, -1})
	if err2 != nil {
		log.Fatal(err2)
	}
	callnumber := 0
	var r, g, b float64
	ret := func(coord []*v3.Matrix) ([]float64, error) {
		angles, err := chem.RamaCalc(coord[0], ramasets)
		ra := make([]float64, 0, len(angles)*2)
		for _, v := range angles {
			ra = append(ra, v[0])
			ra = append(ra, v[1])
		}
		if err != nil {
			return nil, fmt.Errorf("Ramachandran f: " + err.Error())
		}
		ramas = append(ramas, ra)
		if frames != 0 {
			r, g, b = colors(callnumber, frames)
			RGBs = append(RGBs, []float64{r, g, b})
		}
		callnumber++
		return ra, nil //dummy value
	}
	final := func() [][]float64 {
		return RGBs
	}
	return ret, final
}

//Nice trick from icza on StackOverflow
//https://stackoverflow.com/questions/39544571/golang-round-to-nearest-0-05
//Note that id doesn't consider that when the rounded digit is five, numbers are promoted
//to the nearest even number. Still it's more than good enough for our purposes.
func round(x float64) float64 {
	return float64(int64(x + 0.5))
}

//These are the same goChem/chem functions. I will eventually make them
//exported in goChem and then replace these for calls to the original ones.
func colors(key, steps int) (r, g, b float64) {
	norm := 260.0 / float64(steps)
	hp := float64((float64(key) * norm) + 20.0)
	var h float64
	if hp < 55 {
		h = hp - 20.0
	} else {
		h = hp + 20.0
	}
	//	fmt.Println("HUE", h, hp)
	s := 1.0
	v := 1.0
	r, g, b = iHVS2RGB(h, v, s)
	return round(r), round(g), round(b) //slight change to the original to return rounded floats
}

//Here there is a change from the original gochem function to
//return float64 instead of uint8
//takes hue (0-360), v and s (0-1), returns r,g,b (0-255)
func iHVS2RGB(h, v, s float64) (float64, float64, float64) {
	var i, f, p, q, t float64
	var r, g, b float64
	maxcolor := 255.0
	conversion := maxcolor * v
	if s == 0.0 {
		return (conversion), (conversion), (conversion)
	}
	//conversion:=math.Sqrt(3*math.Pow(maxcolor,2))*v
	h = h / 60
	i = math.Floor(h)
	f = h - i
	p = v * (1 - s)
	q = v * (1 - s*f)
	t = v * (1 - s*(1-f))
	switch int(i) {
	case 0:
		r = v
		g = t
		b = p
	case 1:
		r = q
		g = v
		b = p
	case 2:
		r = p
		g = v
		b = t
	case 3:
		r = p
		g = q
		b = v
	case 4:
		r = t
		g = p
		b = v
	default: //case 5
		r = v
		g = p
		b = q
	}

	r = r * conversion
	g = g * conversion
	b = b * conversion
	return (r), (g), (b)
}
