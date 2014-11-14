package prot

import (
	"fmt"
	"io/ioutil"
	"math"

	"code.google.com/p/plotinum/plot"
	"code.google.com/p/plotinum/plotter"
	"code.google.com/p/plotinum/plotutil"
	"github.com/jgcarvalho/PhageAnalysis/pep"
)

type Protein struct {
	Name        string
	Seq         string
	Score       []float64
	RandomScore []float64
	DiffScore   []float64
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
	pro.Length = len(pro.Seq)
	score := make([]float64, pro.Length)
	l := len(p[0].Seq)
	ps := make([]float64, l)
	for i := 0; i < pro.Length-(l-1); i++ {
		for j := 0; j < len(p); j++ {
			ps = match(pro.Seq[i:i+l], p[j].Seq)
			sum(score[i:i+l], ps, p[j].Freq)
		}
	}
	return score
}

func (pro *Protein) plot() {
	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.Title.Text = pro.Name
	p.X.Label.Text = "Residue"
	p.Y.Label.Text = "Score"
	ptsdiff := make(plotter.XYs, len(pro.Score))
	ptsscore := make(plotter.XYs, len(pro.Score))
	ptsrandom := make(plotter.XYs, len(pro.Score))
	for i := 0; i < len(ptsdiff); i++ {
		ptsdiff[i].X, ptsscore[i].X, ptsrandom[i].X = float64(i), float64(i), float64(i)
		ptsdiff[i].Y, ptsscore[i].Y, ptsrandom[i].Y = pro.DiffScore[i], pro.Score[i], pro.RandomScore[i]
	}

	err = plotutil.AddLinePoints(p,
		"Difference", ptsdiff,
		"Selected", ptsscore,
		"Random", ptsrandom)
	if err != nil {
		panic(err)
	}
	if err := p.Save(16, 8, pro.Name+".svg"); err != nil {
		panic(err)
	}
}

func (pro *Protein) save() {
	out := fmt.Sprintf("%s Score Total: %f Length: %d\n\n", pro.Name, pro.Total, pro.Length)
	out += fmt.Sprintf("Score, Random Score, DiffScore\n")
	for i := 0; i < len(pro.Score); i++ {
		out += fmt.Sprintf("%f, %f, %f\n", pro.Score[i], pro.RandomScore[i], pro.DiffScore[i])
	}
	err := ioutil.WriteFile(pro.Name+".dat", []byte(out), 0644)
	if err != nil {
		fmt.Println("Erro ao salvar o arquivo")
		panic(err)
	}
}

func (pro *Protein) Analysis(peps []pep.Peptide, randpeps []pep.Peptide) {
	// fmt.Print(match(s1, s2))

	pro.Score = pro.calcScore(peps)
	pro.RandomScore = pro.calcScore(randpeps)
	n := len(randpeps)
	for i := 0; i < pro.Length; i++ {
		pro.Score[i] = pro.Score[i] / (float64(n) / 100.0)
		pro.RandomScore[i] = pro.RandomScore[i] / (float64(n) / 100.0)
	}
	pro.DiffScore = make([]float64, pro.Length)
	for i := 0; i < pro.Length; i++ {
		pro.DiffScore[i] = pro.Score[i] - pro.RandomScore[i]
		pro.Total += pro.DiffScore[i]
	}
	pro.save()
	pro.plot()

	// score := pro.calcScore(p)
	// for i := 0; i < len(score)-1; i++ {
	// 	fmt.Printf("%f, ", score[i])
	// }
	// fmt.Printf("%f\n", score[len(score)-1])
}
