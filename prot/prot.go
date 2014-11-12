package prot

import (
	"fmt"
	"math"

	"github.com/jgcarvalho/PhageAnalysis/pep"
)

type Protein struct {
	Name        string
	Seq         string
	Score       []float64
	RandomScore []float64
	Total       float64
	Length      int
}

type Proteins []Protein

func match(s1, s2 string) []float64 {
	if len(s1) != len(s2) {
		panic("MATCH com strings de tamanho diferente")
	}
	max := len(s1)
	x := 0
	score := make([]float64, max)
	for i := 0; i < max; i++ {
		if s1[i] == s2[i] {
			score[i] = 1
			x++
		}
	}
	// Testar value^2 e outros expoentes
	// value := (float64(x) - 1.0) / (float64(max) - 1.0)
	value := math.Pow((float64(x)-1.0)/(float64(max)-1.0), 2)
	for i := 0; i < max; i++ {
		score[i] *= value
	}
	return score
}
func sum(score, s []float64, freq int) {
	if len(score) != len(s) {
		panic("SUM parametros com tamanho diferente")
	}
	for i := 0; i < len(score); i++ {
		score[i] += s[i] * float64(freq)
	}
}

func (pro *Protein) calcScore(p []pep.Peptide) []float64 {
	score := make([]float64, len(pro.Seq))
	l := len(p[0].Seq)
	ps := make([]float64, l)
	for i := 0; i < len(pro.Seq)-(l-1); i++ {
		for j := 0; j < len(p); j++ {
			ps = match(pro.Seq[i:i+l], p[j].Seq)
			sum(score[i:i+l], ps, p[j].Freq)
		}
	}
	return score
}

func Teste(pro Protein, p []pep.Peptide) {
	// fmt.Print(match(s1, s2))
	score := pro.calcScore(p)
	for i := 0; i < len(score)-1; i++ {
		fmt.Printf("%f, ", score[i])
	}
	fmt.Printf("%f\n", score[len(score)-1])
}
