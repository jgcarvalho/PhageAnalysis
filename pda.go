package main

import (
	"fmt"
	"os"
	"strings"

	"github.com/jgcarvalho/PhageAnalysis/pep"

	"code.google.com/p/biogo/align"
	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/io/seqio"
	"code.google.com/p/biogo/io/seqio/fastq"
	"code.google.com/p/biogo/seq"
	"code.google.com/p/biogo/seq/linear"
	"code.google.com/p/biogo/seq/sequtils"
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
}

func translate(in alphabet.Slice) (*linear.Seq, bool) {
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
	pep := linear.NewSeq("pep", alphabet.BytesToLetters([]byte(p)), alphabet.Protein)
	return pep, ok
}

func readMFastQ(fn string) ([]seq.Sequence, error) {
	fasta, err := os.Open(fn)
	if err != nil {
		fmt.Println("Erro ao ler o arquivo", err)
	}
	var s []alphabet.Letter
	var sequence []seq.Sequence
	//FASTQ encoding Sanger, Illumina, Solexa...
	t := linear.NewSeq("", s, alphabet.DNAredundant)
	reader := fastq.NewReader(fasta, t)
	scanner := seqio.NewScanner(reader)
	for scanner.Next() {
		s := scanner.Seq()
		sequence = append(sequence, s)
	}
	return sequence, nil
}

func getPeptides(dna []seq.Sequence, template seq.Sequence) (peptides, unreliable []pep.Peptide) {
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

		if len(aln) == 1 {
			fa := align.Format(template, v, aln, '-')
			if strings.Count(fmt.Sprint(fa[0]), "N") == pepLen {
				for i, v := range fmt.Sprint(fa[0]) {
					if v == 'N' {
						tmp := fa[1].Slice(i, i+pepLen)
						p.Seq, isok = translate(tmp)
						p.Freq = 1
						if strings.Contains(p.Seq.String(), "*") || (!isok) {
							unreliable = append(unreliable, p)
						} else {
							peptides = append(peptides, p)
						}

						break
					}
				}

			}

		}

	}
	peptides = pep.Rank(peptides)
	return peptides, unreliable
}

func main() {
	//tamanho do peptideo
	pepLen := 12 * 3
	var dna []seq.Sequence

	fpSeq := alphabet.BytesToLetters([]byte("CCTCTCTATGGGCAGTCGGTGATCCTTTCTATTCTCACTCT"))
	forw := linear.NewSeq("Forward Primer", fpSeq, alphabet.DNAredundant)
	// fmt.Println(forw)

	rpSeq := alphabet.BytesToLetters([]byte("CCGAACCTCCACC"))
	reverse := linear.NewSeq("Reverse Primer", rpSeq, alphabet.DNAredundant)
	reverse.RevComp()
	// fmt.Println(reverse)

	varSeq := alphabet.BytesToLetters([]byte(strings.Repeat("N", pepLen)))
	template := linear.NewSeq("template", varSeq, alphabet.DNAredundant)
	// fmt.Println(template)

	sequtils.Join(template, forw, seq.Start)
	sequtils.Join(template, reverse, seq.End)
	// fmt.Println(template)

	files := []string{"./sample1.fastq", "./sample2.fastq", "./sample3.fastq", "./sample4.fastq", "./sample5.fastq"}
	for _, f := range files {
		tmp, err := readMFastQ(f)
		if err != nil {
			fmt.Println("ERRO", err)
		}
		dna = append(dna, tmp...)
	}

	peptides, _ := getPeptides(dna, template)
	for _, p := range peptides {
		fmt.Println(p)
	}

	pep.Teste(peptides)
}
