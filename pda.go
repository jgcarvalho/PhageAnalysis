package main

import (
	"fmt"
	"os"
	"strings"

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

func RevComp(s string) string {
	runes := []rune(s)
	for i, v := range runes {
		if v == 'C' {
			runes[i] = 'G'
		} else if v == 'G' {
			runes[i] = 'C'
		} else if v == 'A' {
			runes[i] = 'T'
		} else if v == 'T' {
			runes[i] = 'A'
		}
	}
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}

func Translate(in alphabet.Slice) {
	s := fmt.Sprintf("%s", in)
	fmt.Println(s)
	s = RevComp(s)
	fmt.Println(s)
	for i := len(s); i > 0; i -= 3 {
		fmt.Printf(" %s ", codon[s[i-3:i]])
		// if s[i:i+3] == "ACT" {
		// 	fmt.Printf(codon)
		// }
	}
	fmt.Println()
}

func ReadMFastQ(fn string) ([]seq.Sequence, error) {
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
		//fmt.Println(s)
		sequence = append(sequence, s)
	}
	// for i, seq := reader.Read(){
	//     fmt.Println(i, seq)
	// }
	//fmt.Println("Read -> ", seq.Alphabet())
	return sequence, nil
}

func main() {
	// s := []alphabet.Letter{"ACTG"}
	fpSeq := alphabet.BytesToLetters([]byte("CCTCTCTATGGGCAGTCGGTGATCCTTTCTATTCTCACTCT"))
	forw := linear.NewSeq("Forward Primer", fpSeq, alphabet.DNAredundant)
	forw.RevComp()
	fmt.Println(forw)

	rpSeq := alphabet.BytesToLetters([]byte("CCGAACCTCCACC"))
	reverse := linear.NewSeq("Reverse Primer", rpSeq, alphabet.DNAredundant)
	fmt.Println(reverse)

	varSeq := alphabet.BytesToLetters([]byte("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"))
	variable := linear.NewSeq("Variable", varSeq, alphabet.DNAredundant)
	fmt.Println(variable)

	sequtils.Join(variable, forw, seq.End)
	sequtils.Join(variable, reverse, seq.Start)
	fmt.Println(variable)

	dna, err := ReadMFastQ("/home/zeh/gocode/src/github.com/jgcarvalho/PhageAnalysis/sample2.fastq")
	if err != nil {
		fmt.Println("ERRO", err)
	}

	sw := align.FittedAffine{
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

	for _, v := range dna {
		//fmt.Println(i, v)
		aln, err := sw.Align(variable, v)
		if err != nil {
			fmt.Println(err)
		}
		if len(aln) == 1 {
			fa := align.Format(variable, v, aln, '-')
			if strings.Count(fmt.Sprint(fa[0]), "N") == 36 {
				for i, v := range fmt.Sprint(fa[0]) {
					if v == 'N' {
						tmp := fa[1].Slice(i, i+36)
						fmt.Printf("%s\n", tmp)
						fmt.Printf("%T\n", tmp)
						Translate(tmp)
						break
					}
				}
				fmt.Printf("%s\n", aln)
				fmt.Printf("%s\n%s\n", fa[0], fa[1])
			}
		}
		// if i > 5 {
		// 	break
		// }
	}

}
