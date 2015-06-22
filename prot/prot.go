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

// Dixon Q test table (David B. Rorabacher Anal. Chem., 1991, 63 (2), 139-146• DOI: 10.1021/ac00002a010 )
// 3 0.886 0.941 0.970 0.976 0.988 0.994
// 4 0.679 0.765 0.829 0.846 0.889 0.926
// 5 0.557 0.642 0.7 10 0.729 0.780 0.821
// 6 0.482 0.560 0.625 0.644 0.698 0.740
// 7 0.434 0.507 0.568 0.586 0.637 0.680
// 8 0.399 0.468 0.526 0.543 0.590 0.634
// 9 0.370 0.437 0.493 0.510 0.555 0.598
// 10 0.349 0.412 0.466 0.483 0.527 0.568
// 11 0.332 0.392 0.444 0.460 0.502 0.542
// 12 0.318 0.376 0.426 0.441 0.482 0.522
// 13 0.305 0.361 0.410 0.425 0.465 0.503
// 14 0.294 0.349 0.396 0.411 0.450 0.488
// 15 0.285 0.338 0.384 0.399 0.438 0.475
// 16 0.277 0.329 0.374 0.388 0.426 0.463
// 17 0.269 0.320 0.365 0.379 0.416 0.452
// 18 0.263 0.313 0.356 0.370 0.407 0.442
// 19 0.258 0.306 0.349 0.363 0.398 0.433
// 20 0.252 0.300 0.342 0.356 0.391 0.425
// 21 0.247 0.295 0.337 0.350 0.384 0.418
// 22 0.242 0.290 0.331 0.344 0.378 0.411
// 23 0.238 0.285 0.326 0.338 0.372 0.404
// 24 0.234 0.281 0.321 0.333 0.367 0.399
// 25 0.230 0.277 0.3 17 0.329 0.362 0.393
// 26 0.227 0.273 0.312 0.324 0.357 0.388
// 27 0.224 0.269 0.308 0.320 0.353 0.384
// 28 0.220 0.266 0.305 0.316 0.349 0.380
// 29 0.218 0.263 0.301 0.312 0.345 0.376
// 30 0.215 0.260 0.298 0.309 0.341 0.372
var q99 = [28]float64{0.994, 0.926, 0.821, 0.740, 0.680, 0.634, 0.598, 0.568, 0.542, 0.522, 0.503, 0.488, 0.475, 0.463, 0.452, 0.442, 0.433, 0.425, 0.418, 0.411, 0.404, 0.399, 0.393, 0.388, 0.384, 0.380, 0.376, 0.372}
var q98 = [28]float64{0.988, 0.889, 0.780, 0.698, 0.637, 0.590, 0.555, 0.527, 0.502, 0.482, 0.465, 0.450, 0.438, 0.426, 0.416, 0.407, 0.398, 0.391, 0.384, 0.378, 0.372, 0.367, 0.362, 0.357, 0.353, 0.349, 0.345, 0.341}
var q96 = [28]float64{0.976, 0.846, 0.729, 0.644, 0.586, 0.543, 0.510, 0.483, 0.460, 0.441, 0.425, 0.411, 0.399, 0.388, 0.379, 0.370, 0.363, 0.356, 0.350, 0.344, 0.338, 0.333, 0.329, 0.324, 0.320, 0.316, 0.312, 0.309}
var q95 = [28]float64{0.970, 0.829, 0.710, 0.625, 0.568, 0.526, 0.493, 0.466, 0.444, 0.426, 0.410, 0.396, 0.384, 0.374, 0.365, 0.356, 0.349, 0.342, 0.337, 0.331, 0.326, 0.321, 0.317, 0.312, 0.308, 0.305, 0.301, 0.298}
var q90 = [28]float64{0.941, 0.765, 0.642, 0.560, 0.507, 0.468, 0.437, 0.412, 0.392, 0.376, 0.361, 0.349, 0.338, 0.329, 0.320, 0.313, 0.306, 0.300, 0.295, 0.290, 0.285, 0.281, 0.277, 0.273, 0.269, 0.266, 0.263, 0.260}
var q80 = [28]float64{0.886, 0.679, 0.557, 0.482, 0.434, 0.399, 0.370, 0.349, 0.332, 0.318, 0.305, 0.294, 0.285, 0.277, 0.269, 0.263, 0.258, 0.252, 0.247, 0.242, 0.238, 0.234, 0.230, 0.227, 0.224, 0.220, 0.218, 0.215}

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
	Qvalue           []float64
	// Total            []float64
	// TotalMean        float64
	// TotalSD          float64
	// TotalError       float64
	TotalScore     float64
	TotalDiffScore float64
	Length         int
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
	// out := fmt.Sprintf("%s Score Total: %f Length: %d\n\n", pro.Name, pro.Total, pro.Length)
	out := fmt.Sprintf("%s Score Total: %f Length: %d\n\n", pro.Name, pro.TotalScore, pro.Length)
	out += fmt.Sprintf("NRes, Res, Score, Random Score Mean, Random Score SD, Random Score Error, " +
		"Diff Score Mean, Diff Score SD, Diff Score Error, Q-Value\n")
	for i := 0; i < len(pro.Score); i++ {
		out += fmt.Sprintf("%d, %c, %f, %f, %f, %f, %f, %f, %f, %f\n", i+1, pro.Seq[i], pro.Score[i],
			pro.RandomScoreMean[i], pro.RandomScoreSD[i], pro.RandomScoreError[i],
			pro.DiffScoreMean[i], pro.DiffScoreSD[i], pro.DiffScoreError[i], pro.Qvalue[i])
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
	pro.Qvalue = make([]float64, pro.Length)

	// pro.Total = make([]float64, nlibrand)

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
			// pro.Total[j] += pro.DiffScore[j][i]
		}
		// pro.DiffScore = make([]float64, pro.Length)
		// for i := 0; i < pro.Length; i++ {
		// 	pro.DiffScore[j][i] = pro.Score[i] - pro.RandomScore[j][i]
		// 	pro.Total[j] += pro.DiffScore[j][i]
		// }
		// fmt.Printf("Protein: %s Total %d = %.3f\n", pro.Name, j, pro.Total[j])
	}

	// pro.TotalMean, pro.TotalSD = stat.MeanStdDev(pro.Total, nil)
	// pro.TotalError = stat.StdErr(pro.TotalSD, float64(nlibrand))
	fmt.Println("Done")

	diff := make([]float64, nlibrand)
	randscore := make([]float64, nlibrand)
	for i := 0; i < pro.Length; i++ {
		max, min := 0.0, 0.0
		for j := 0; j < nlibrand; j++ {
			diff[j] = pro.DiffScore[j][i]
			pro.DiffScoreMean[i], pro.DiffScoreSD[i] = stat.MeanStdDev(diff, nil)
			pro.DiffScoreError[i] = stat.StdErr(pro.DiffScoreSD[i], float64(nlibrand))
			randscore[j] = pro.RandomScore[j][i]
			pro.RandomScoreMean[i], pro.RandomScoreSD[i] = stat.MeanStdDev(randscore, nil)
			pro.RandomScoreError[i] = stat.StdErr(pro.RandomScoreSD[i], float64(nlibrand))
			sort.Float64s(randscore)
			// fmt.Println("random", randscore)
			min, max = randscore[0], randscore[len(randscore)-1]
		}
		if pro.Score[i] < max {
			pro.Qvalue[i] = 0.0
		} else {
			pro.Qvalue[i] = (pro.Score[i] - max) / (pro.Score[i] - min)
		}
		if pro.Qvalue[i] > q99[nlibrand-3] {
			pro.TotalScore += pro.Score[i]
			pro.TotalDiffScore += pro.DiffScoreMean[i]
		}
		// fmt.Printf("Qvalue min %f, max %f, score %f, value %f\n", min, max, pro.Score[i], pro.Qvalue[i])

	}
	fmt.Printf("Protein: %s, Score: %f, Diff: %f\n", pro.Name, pro.TotalScore, pro.TotalDiffScore)

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
	// return p[i].TotalMean < p[j].TotalMean
	return p[i].TotalScore < p[j].TotalScore
}

func (p Proteins) Swap(i, j int) {
	p[i], p[j] = p[j], p[i]
}

func (p Proteins) String() string {
	// output := fmt.Sprintf("#Rank\tID\tLength\tTotalMeanScore\tTotalSDScore\tTotalErrorScore\n")
	output := fmt.Sprintf("#Rank\tID\tLength\tTotalScore\tTotalDiffScore\n")

	sort.Sort(sort.Reverse(p))
	for i := 0; i < len(p); i++ {
		if p[i].Length > 1 {
			// output += fmt.Sprintf("%d\t%s\t%d\t%f\t%f\t%f\n", i, p[i].Name, p[i].Length, p[i].TotalMean, p[i].Qvalue, p[i].TotalError)
			output += fmt.Sprintf("%d\t%s\t%d\t%f\t%f\n", i, p[i].Name, p[i].Length, p[i].TotalScore, p[i].TotalDiffScore)

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
