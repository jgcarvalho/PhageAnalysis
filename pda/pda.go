package pda

import (
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/jgcarvalho/PhageAnalysis/pep"
	"github.com/jgcarvalho/PhageAnalysis/prot"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/io/seqio/fastq"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/seq/sequtils"
)

var codon = map[string]string{
	"AAA": "K", "AAT": "N", "AAC": "N", "AAG": "K",
	"ATA": "I", "ATT": "I", "ATC": "I", "ATG": "M",
	"ACA": "T", "ACT": "T", "ACC": "T", "ACG": "T",
	"AGA": "R", "AGT": "S", "AGC": "S", "AGG": "R",
	"TAA": "*", "TAT": "Y", "TAC": "Y", "TAG": "*",
	"TTA": "L", "TTT": "F", "TTC": "F", "TTG": "L",
	"TCA": "S", "TCT": "S", "TCC": "S", "TCG": "S",
	"TGA": "*", "TGT": "C", "TGC": "C", "TGG": "W",
	"CAA": "Q", "CAT": "H", "CAC": "H", "CAG": "Q",
	"CTA": "L", "CTT": "L", "CTC": "L", "CTG": "L",
	"CCA": "P", "CCT": "P", "CCC": "P", "CCG": "P",
	"CGA": "R", "CGT": "R", "CGC": "R", "CGG": "R",
	"GAA": "E", "GAT": "D", "GAC": "D", "GAG": "E",
	"GTA": "V", "GTT": "V", "GTC": "V", "GTG": "V",
	"GCA": "A", "GCT": "A", "GCC": "A", "GCG": "A",
	"GGA": "G", "GGT": "G", "GGC": "G", "GGG": "G",
}

/*
var sw = align.FittedAffine{
	Matrix: [][]int{
		{0, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50},
		{-50, 50, -50, -50, -4, 1, -4, -1, -4, 1, -4, -1, -4, -1, -4, 0},
		{-50, -50, 50, -50, -50, -4, 1, -1, -4, -4, 1, -1, -4, -4, -1, 0},
		{-50, -50, -50, 50, -50, -2, -2, -1, -4, -2, -2, -1, -4, -3, -3, 0},
		{-50, -50, -50, -50, 50, 1, 1, -1, -4, -4, -4, -4, 1, -1, -1, 0},
		{-50, 1, -4, -2, 1, -1, -2, -1, -4, -2, -4, -3, -2, -1, -3, 0},
		{-50, -4, 1, -2, 1, -2, -1, -1, -4, -4, -2, -3, -2, -3, -1, 0},
		{-50, -1, -1, -1, -1, -1, -1, -1, -4, -3, -3, -2, -3, -2, -2, 0},
		{-50, -4, -4, -4, -4, -4, -4, -4, 5, 1, 1, -1, 1, -1, -1, 0},
		{-50, 1, -4, -2, -4, -2, -4, -3, 1, -1, -2, -1, -2, -1, -3, 0},
		{-50, -4, 1, -2, -4, -4, -2, -3, 1, -2, -1, -1, -2, -3, -1, 0},
		{-50, -1, -1, -1, -4, -3, -3, -2, -1, -1, -1, -1, -3, -2, -2, 0},
		{-50, -4, -4, -4, 1, -2, -2, -3, 1, -2, -2, -3, -1, -1, -1, 0},
		{-50, -1, -4, -3, -1, -1, -3, -2, -1, -1, -3, -2, -1, -1, -2, 0},
		{-50, -4, -1, -3, -1, -3, -1, -2, -1, -3, -1, -2, -1, -2, -1, 0},
		{-50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	},
	GapOpen: 0,
}*/

/*func Translate(in alphabet.Slice) (string, bool) {
	s := fmt.Sprintf("%s", in)
	ok := true
	// fmt.Println(s)
	var p string
	for i := 0; i < len(s); i += 3 {
		// test NNK codons
		if s[i+2] != 'T' && s[i+2] != 'G' {
			ok = false
		}
		p += codon[s[i:i+3]]
	}
	// pep := linear.NewSeq("pep", alphabet.BytesToLetters([]byte(p)), alphabet.Protein)
	pep := p
	return pep, ok
}*/

func reverse(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}

func translate(dna string) string {
	var p string
	for i := 0; i < len(dna); i += 3 {
		/*		if dna[i+2] != 'T' && dna[i+2] != 'G' {
				fmt.Println("not NNK")
			}*/
		p += codon[dna[i:i+3]]
	}

	return p
}

func ReadMFastQ(fn string) ([]seq.Sequence, error) {
	fq, err := os.Open(fn)
	if err != nil {
		fmt.Println("Erro ao ler o arquivo", err)
	}
	var s []alphabet.Letter
	var sequence []seq.Sequence
	//FASTQ encoding Sanger, Illumina, Solexa...
	t := linear.NewSeq("", s, alphabet.DNAredundant)
	reader := fastq.NewReader(fq, t)
	scanner := seqio.NewScanner(reader)
	for scanner.Next() {
		s := scanner.Seq()
		sequence = append(sequence, s)
	}
	return sequence, nil
}

/*//Pode substituir parte da função acima
func ReadMDna(fq *os.File) ([]seq.Sequence, error) {
	var s []alphabet.Letter
	var sequence []seq.Sequence
	//FASTQ encoding Sanger, Illumina, Solexa...
	t := linear.NewSeq("", s, alphabet.DNAredundant)
	reader := fastq.NewReader(fq, t)
	scanner := seqio.NewScanner(reader)
	for scanner.Next() {
		s := scanner.Seq()
		sequence = append(sequence, s)
	}
	fmt.Println(sequence)
	return sequence, nil
}*/

func getpep(dna, template string) (string, bool, int) {
	score := 0
	template = strings.Repeat("#", len(dna)) + template + strings.Repeat("#", len(dna))
	pepLen := strings.Count(template, "X")
	firstN := strings.Index(template, "X")
	lastN := strings.LastIndex(template, "X")
	maxscore := len(dna) - pepLen

	bestScore := maxscore
	bestPep := ""

	var valid bool
	var pep []byte

	//fmt.Println("Max Score", maxscore)
	for j := lastN + 1 - len(dna); j < firstN; j++ {
		score = maxscore
		pep = []byte{}
		for i, v := range dna {
			if byte(v) == template[i+j] {
				score -= 1
			}
			if template[i+j] == 'X' {
				pep = append(pep, byte(v))
			}
		}
		if score < bestScore {
			bestScore = score
			bestPep = string(pep)
		}
		//fmt.Println(j, score)
	}

	bestPep = translate(bestPep)
	if bestScore > 1 || strings.Contains(bestPep, "*") {
		valid = false
	} else {
		valid = true
	}

	if pepLen/3 != len(bestPep) {
		fmt.Println("incorrect peptide")
		fmt.Println(dna)
		valid = false
	}
	return bestPep, valid, bestScore
}

func GetPeptides(dna []seq.Sequence, template seq.Sequence, translatedStrand bool, readLen int) (peptides, unreliable []pep.Peptide) {
	// count the number of N to determine the peptide length
	//pepLen := strings.Count(fmt.Sprintf("%s", template.Slice()), "X")

	// fmt.Println("New function")
	// fmt.Println("count X", pepLen)
	//template.RevComp()
	templateSeq := fmt.Sprintf("%s", template.Slice())
	/*	dnaSeq := fmt.Sprintf("%s", dna[0].Slice())
		maxScore := len(dnaSeq) - pepLen
		fmt.Println(templateSeq)
		fmt.Println(dnaSeq)
		fmt.Println(maxScore)
		pep, valid, score := getpep(dnaSeq, templateSeq)
		fmt.Println(pep, valid, score)*/

	var p pep.Peptide
	//var pepSeq string
	var valid bool
	var score int
	var dnaSeq string

	for _, v := range dna {
		if !translatedStrand {
			v.RevComp()
		}
		dnaSeq = fmt.Sprintf("%s", v.Slice())
		if len(dnaSeq) > readLen {
			fmt.Println("skip read: length > ", readLen)
			continue
		}
		if strings.Contains(dnaSeq, "N") {
			fmt.Println("skip read: undefined base N")
			continue
		}
		p.Seq, valid, score = getpep(dnaSeq, templateSeq)
		p.Freq = 1
		fmt.Println(p.Seq, valid, score)
		if valid {
			peptides = append(peptides, p)
		} else {
			unreliable = append(unreliable, p)
		}
	}
	peptides = pep.Rank(peptides)
	return
}

/*func GetPeptidesOld(dna []seq.Sequence, template seq.Sequence) (peptides, unreliable []pep.Peptide) {
	// count the number of N to determine the peptide length
	pepLen := strings.Count(fmt.Sprintf("%s", template.Slice()), "N")

	var p pep.Peptide
	var isok bool
	for _, v := range dna {
		v.RevComp()
		aln, err := sw.Align(template, v)
		if err != nil {
			panic(err)
		}
		fmt.Println(aln)

		if len(aln) == 1 {
			fa := align.Format(template, v, aln, '-')
			if strings.Count(fmt.Sprint(fa[0]), "N") == pepLen {
				for i, v := range fmt.Sprint(fa[0]) {
					if v == 'N' {
						tmp := fa[1].Slice(i, i+pepLen)
						p.Seq, isok = Translate(tmp)
						p.Freq = 1
						// if strings.Contains(p.Seq.String(), "*") || (!isok) {
						if strings.Contains(p.Seq, "*") || (!isok) {
							unreliable = append(unreliable, p)
						} else {
							peptides = append(peptides, p)
						}

						break
					}
				}

			}
		} else {
			fmt.Println("problema alinhament")
			fmt.Println(len(aln))
		}

	}
	peptides = pep.Rank(peptides)
	return peptides, unreliable
}*/

func ReadTranslatedPeptides(fn string, penLen int) ([]pep.Peptide, error) {
	fa, err := os.Open(fn)
	if err != nil {
		fmt.Println("Erro ao ler o arquivo", err)
	}
	var s []alphabet.Letter
	var sequence []seq.Sequence
	//FASTQ encoding Sanger, Illumina, Solexa...
	t := linear.NewSeq("", s, alphabet.Protein)
	reader := fasta.NewReader(fa, t)
	scanner := seqio.NewScanner(reader)
	for scanner.Next() {
		s := scanner.Seq()
		sequence = append(sequence, s)
	}
	peptides := make([]pep.Peptide, len(sequence))
	for i, s := range sequence {
		// fmt.Println(s.Description())
		freq, err := strconv.Atoi(strings.TrimSuffix(strings.TrimPrefix(s.Description(), "[Freq: "), "]"))
		if err != nil {
			fmt.Println("Erro ao ler a frequencia do peptideo", err)
			return nil, err
		}
		peptides[i].Name = s.Name()
		peptides[i].Freq = freq
		peptides[i].Seq = fmt.Sprintf("%s", s.Slice())
	}
	return peptides, err
}

func ReadMProteins(fn string) ([]prot.Protein, error) {
	fa, err := os.Open(fn)
	if err != nil {
		fmt.Println("Erro ao ler o arquivo", err)
	}
	var s []alphabet.Letter
	var sequence []seq.Sequence
	//FASTQ encoding Sanger, Illumina, Solexa...
	t := linear.NewSeq("", s, alphabet.Protein)
	reader := fasta.NewReader(fa, t)
	scanner := seqio.NewScanner(reader)
	for scanner.Next() {
		s := scanner.Seq()
		sequence = append(sequence, s)
	}
	proteins := make([]prot.Protein, len(sequence))
	for i, s := range sequence {
		proteins[i].Name = fmt.Sprintf("%s", s.Name())
		proteins[i].Seq = fmt.Sprintf("%s", s.Slice())
	}
	return proteins, err
}

func CreateTemplate(fp, rp string, pepLen int) *linear.Seq {
	// fpSeq := alphabet.BytesToLetters([]byte("CCTCTCTATGGGCAGTCGGTGATCCTTTCTATTCTCACTCT"))
	fpSeq := alphabet.BytesToLetters([]byte(fp))
	forw := linear.NewSeq("Forward Primer", fpSeq, alphabet.DNAredundant)
	// forw.RevComp()
	// fmt.Println(forw)

	// rpSeq := alphabet.BytesToLetters([]byte("CCGAACCTCCACC"))
	rpSeq := alphabet.BytesToLetters([]byte(rp))
	reverse := linear.NewSeq("Reverse Primer", rpSeq, alphabet.DNAredundant)
	reverse.RevComp()
	// fmt.Println(reverse)

	varSeq := alphabet.BytesToLetters([]byte(strings.Repeat("X", pepLen)))
	template := linear.NewSeq("template", varSeq, alphabet.DNAredundant)
	// fmt.Println(template)

	sequtils.Join(template, forw, seq.Start)
	sequtils.Join(template, reverse, seq.End)
	// sequtils.Join(template, reverse, seq.Start)
	// sequtils.Join(template, forw, seq.End)
	// fmt.Println(template)

	return template
}
