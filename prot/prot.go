package prot

import (
	"fmt"
	"io/ioutil"
	"math"
	"sort"

	"code.google.com/p/plotinum/plot"
	"code.google.com/p/plotinum/plotter"
	"code.google.com/p/plotinum/plotutil"
	"github.com/gonum/stat"
	"github.com/jgcarvalho/PhageAnalysis/pep"
)

type Protein struct {
	Name             string
	Seq              string
	Score            []float64
	RandomScore      [][]float64
	RandomScoreMean  []float64
	RandomScoreSD    []float64
	RandomScoreError []float64
	DiffScore        [][]float64
	DiffScoreMean    []float64
	DiffScoreSD      []float64
	DiffScoreError   []float64
	Total            []float64
	TotalMean        float64
	TotalSD          float64
	TotalError       float64
	Length           int
}

type Proteins []Protein

func match(s1, s2 string, exp float64) []float64 {
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
	value := math.Pow((float64(x)-1.0)/(float64(max)-1.0), exp)
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

func (pro *Protein) calcScore(p []pep.Peptide, exp float64) []float64 {
	pro.Length = len(pro.Seq)
	score := make([]float64, pro.Length)
	l := len(p[0].Seq)
	ps := make([]float64, l)
	for i := 0; i < pro.Length-(l-1); i++ {
		for j := 0; j < len(p); j++ {
			ps = match(pro.Seq[i:i+l], p[j].Seq, exp)
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
		ptsdiff[i].Y, ptsscore[i].Y, ptsrandom[i].Y = pro.DiffScoreMean[i], pro.Score[i], pro.RandomScoreMean[i]
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
	out += fmt.Sprintf("NRes, Res, Score, Random Score Mean, Random Score SD, Random Score Error, " +
		"Diff Score Mean, Diff Score SD, Diff Score Error\n")
	for i := 0; i < len(pro.Score); i++ {
		out += fmt.Sprintf("%d, %c, %f, %f, %f, %f, %f, %f, %f\n", i+1, pro.Seq[i], pro.Score[i],
			pro.RandomScoreMean[i], pro.RandomScoreSD[i], pro.RandomScoreError[i],
			pro.DiffScoreMean[i], pro.DiffScoreSD[i], pro.DiffScoreError[i])
	}
	err := ioutil.WriteFile(pro.Name+".dat", []byte(out), 0644)
	if err != nil {
		fmt.Println("Erro ao salvar o arquivo")
		panic(err)
	}
}

func (pro *Protein) Analysis(peps []pep.Peptide, randpeps [][]pep.Peptide, exp float64) {
	// fmt.Print(match(s1, s2))
	nlibrand := len(randpeps)

	fmt.Println("Calculating score to protein: ", pro.Name)
	pro.Score = pro.calcScore(peps, exp)

	pro.RandomScore = make([][]float64, nlibrand)
	pro.DiffScore = make([][]float64, nlibrand)

	pro.RandomScoreMean = make([]float64, pro.Length)
	pro.RandomScoreSD = make([]float64, pro.Length)
	pro.RandomScoreError = make([]float64, pro.Length)
	pro.DiffScoreMean = make([]float64, pro.Length)
	pro.DiffScoreSD = make([]float64, pro.Length)
	pro.DiffScoreError = make([]float64, pro.Length)

	pro.Total = make([]float64, nlibrand)

	for i := 0; i < nlibrand; i++ {
		pro.RandomScore[i] = pro.calcScore(randpeps[i], exp)
	}

	// pro.RandomScore = pro.calcScore(randpeps[0])
	// n := pro.Length
	for i := 0; i < pro.Length; i++ {
		// pro.Score[i] = pro.Score[i] / 1000.0
		pro.Score[i] = pro.Score[i] / float64(len(randpeps[0])/10)
	}
	for j := 0; j < nlibrand; j++ {
		pro.DiffScore[j] = make([]float64, pro.Length)
		for i := 0; i < pro.Length; i++ {
			// pro.RandomScore[j][i] = pro.RandomScore[j][i] / 1000.0
			pro.RandomScore[j][i] = pro.RandomScore[j][i] / float64(len(randpeps[0])/10)
			pro.DiffScore[j][i] = pro.Score[i] - pro.RandomScore[j][i]
			pro.Total[j] += pro.DiffScore[j][i]
		}
		// pro.DiffScore = make([]float64, pro.Length)
		// for i := 0; i < pro.Length; i++ {
		// 	pro.DiffScore[j][i] = pro.Score[i] - pro.RandomScore[j][i]
		// 	pro.Total[j] += pro.DiffScore[j][i]
		// }
		fmt.Printf("Total %d = %.3f\n", j, pro.Total[j])
	}

	pro.TotalMean, pro.TotalSD = stat.MeanStdDev(pro.Total, nil)
	pro.TotalError = stat.StdErr(pro.TotalSD, float64(nlibrand))
	fmt.Println("Done")

	diff := make([]float64, nlibrand)
	randscore := make([]float64, nlibrand)
	for i := 0; i < pro.Length; i++ {
		for j := 0; j < nlibrand; j++ {
			diff[j] = pro.DiffScore[j][i]
			pro.DiffScoreMean[i], pro.DiffScoreSD[i] = stat.MeanStdDev(diff, nil)
			pro.DiffScoreError[i] = stat.StdErr(pro.DiffScoreSD[i], float64(nlibrand))
			randscore[j] = pro.RandomScore[j][i]
			pro.RandomScoreMean[i], pro.RandomScoreSD[i] = stat.MeanStdDev(randscore, nil)
			pro.RandomScoreError[i] = stat.StdErr(pro.RandomScoreSD[i], float64(nlibrand))
		}
	}

	pro.save()
	// pro.plot()

	// score := pro.calcScore(p)
	// for i := 0; i < len(score)-1; i++ {
	// 	fmt.Printf("%f, ", score[i])
	// }
	// fmt.Printf("%f\n", score[len(score)-1])
}

func (p Proteins) Len() int {
	return len(p)
}

func (p Proteins) Less(i, j int) bool {
	return p[i].TotalMean < p[j].TotalMean
}

func (p Proteins) Swap(i, j int) {
	p[i], p[j] = p[j], p[i]
}

func (p Proteins) String() string {
	output := fmt.Sprintf("#Rank\tID\tLength\tTotalMeanScore\tTotalSDScore\tTotalErrorScore\n")
	sort.Sort(sort.Reverse(p))
	for i := 0; i < len(p); i++ {
		if p[i].Length > 1 {
			output += fmt.Sprintf("%d\t%s\t%d\t%f\t%f\t%f\n", i, p[i].Name, p[i].Length, p[i].TotalMean, p[i].TotalSD, p[i].TotalError)
		}
	}
	return output
}

func (p Proteins) SaveRank(fn string) {
	out := p.String()
	err := ioutil.WriteFile(fn, []byte(out), 0644)
	if err != nil {
		fmt.Println("Erro ao salvar o arquivo")
		panic(err)
	}
}
