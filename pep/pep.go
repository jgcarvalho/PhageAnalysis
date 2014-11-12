package pep

import (
	"fmt"
	"math/rand"
	"sort"
)

type Peptide struct {
	// Seq  *linear.Seq
	Name string
	Seq  string
	Freq int
}

type selection struct {
	peps []Peptide
}

// func calcFreq(p []Peptide) {
//     rank := make(map[string]Peptide)
//     for i, v := range p {
//         val, exist := rank[v.Seq.String()]
//         if exist {
//             val.Freq += 1
//             rank[v.Seq.String()] = val
//         } else {
//             s := v.Seq.String()
//             rank[s] = p[i]
//         }
//     }
//     // for i, v := range
// }

func Rank(p []Peptide) []Peptide {
	rank := make(map[string]Peptide)
	var sel selection
	for i, v := range p {
		// val, exist := rank[v.Seq.String()]
		val, exist := rank[v.Seq]
		if exist {
			val.Freq++
			// rank[v.Seq.String()] = val
			rank[v.Seq] = val
		} else {
			// s := v.Seq.String()
			s := v.Seq
			rank[s] = p[i]
		}
	}
	sel.peps = make([]Peptide, len(rank))
	i := 0
	for k := range rank {
		sel.peps[i] = rank[k]
		i++
	}
	sort.Sort(sort.Reverse(sel))
	for i := range sel.peps {
		// sel.peps[i].Seq.SetName(fmt.Sprintf("pep %d [Freq: %d]", i+1, sel.peps[i].Freq))
		sel.peps[i].Name = fmt.Sprintf("pep %d [Freq: %d]", i+1, sel.peps[i].Freq)
	}
	return sel.peps
}

func aafreq(p []Peptide) (map[rune]float64, int) {
	count := map[rune]int{
		'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
		'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0,
		'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0,
		'S': 0, 'T': 0, 'V': 0, 'Y': 0, 'W': 0,
	}
	freq := map[rune]float64{
		'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
		'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0,
		'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0,
		'S': 0, 'T': 0, 'V': 0, 'Y': 0, 'W': 0,
	}
	n := 0
	nlib := 0
	for _, v := range p {
		nlib += v.Freq
		// for _, c := range v.Seq.Seq.String() {
		for _, c := range v.Seq {
			count[c] += v.Freq
			n += v.Freq
		}
	}
	// fmt.Println(count)
	for k := range freq {
		freq[k] = float64(count[k]) / float64(n)
	}
	return freq, nlib
}

func randomAA(freq map[rune]float64) string {
	// rand.Seed()
	x := rand.Float64()
	// fmt.Println(x)
	aa := ""
	for k := range freq {
		if x > freq[k] {
			// fmt.Println("not", k)
			x -= freq[k]
		} else {
			// fmt.Println("***", k)
			aa = string(k)
			break
		}
	}
	return aa
}

func randomPep(l int, freq map[rune]float64) string {
	p := ""
	for i := 0; i < l; i++ {
		p += randomAA(freq)
	}
	return p
}

func RandomLibrary(peptides []Peptide) []Peptide {
	freq, nlib := aafreq(peptides)
	// fmt.Println("Library size", nlib)
	randPeptides := make([]Peptide, nlib)
	// pl := len(peptides[0].Seq.String())
	pl := len(peptides[0].Seq)
	// fmt.Println("**********", pl)
	// for k := range freq {
	// 	fmt.Printf("%s=%f, ", string(k), freq[k])
	// }

	for i := 0; i < nlib; i++ {
		// fmt.Println(randomPep(12, freq))
		// randPeptides[i].Seq = linear.NewSeq("rand", alphabet.BytesToLetters([]byte(randomPep(pl, freq))), alphabet.Protein)
		randPeptides[i].Name = "rand"
		randPeptides[i].Seq = randomPep(pl, freq)
		randPeptides[i].Freq = 1
	}
	// Teste para avaliar a similaridade entre as sequencias
	// rfreq, _ := aafreq(randPeptides)
	// for k := range freq {
	// 	fmt.Printf("%s=%f, ", string(k), freq[k])
	// 	fmt.Printf("%s=%f\n ", string(k), rfreq[k])
	// }
	return randPeptides
}

func Teste(peptides []Peptide) {
	randomPeps := RandomLibrary(peptides)
	for _, v := range randomPeps {
		fmt.Println(v)
	}

}

func (p Peptide) String() string {
	return fmt.Sprintf(">%s\n%s\n", p.Name, p.Seq)
}

func (sel selection) Len() int {
	return len(sel.peps)
}

func (sel selection) Less(i, j int) bool {
	return sel.peps[i].Freq < sel.peps[j].Freq
}

func (sel selection) Swap(i, j int) {
	sel.peps[i], sel.peps[j] = sel.peps[j], sel.peps[i]
}
