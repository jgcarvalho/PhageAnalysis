package main

import (
	"fmt"
	"os"

	"code.google.com/p/biogo/align"
	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/io/seqio"
	"code.google.com/p/biogo/io/seqio/fastq"
	"code.google.com/p/biogo/seq"
	"code.google.com/p/biogo/seq/linear"
	"code.google.com/p/biogo/seq/sequtils"
)

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

	dna, err := ReadMFastQ("/home/jgcarvalho/sync/colaboracoes/carlos/sample2.fastq")
	if err != nil {
		fmt.Println("ERRO", err)
	}

	sw := align.FittedAffine{
		Matrix: [][]int{
			{0, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100},
			{-100, 50, -50, -50, -4, 1, -4, -1, -4, 1, -4, -1, -4, -1, -4, 0},
			{-100, -50, 50, -50, -50, -4, 1, -1, -4, -4, 1, -1, -4, -4, -1, 0},
			{-100, -50, -50, 50, -50, -2, -2, -1, -4, -2, -2, -1, -4, -3, -3, 0},
			{-100, -50, -50, -50, 50, 1, 1, -1, -4, -4, -4, -4, 1, -1, -1, 0},
			{-100, 1, -4, -2, 1, -1, -2, -1, -4, -2, -4, -3, -2, -1, -3, 0},
			{-100, -4, 1, -2, 1, -2, -1, -1, -4, -4, -2, -3, -2, -3, -1, 0},
			{-100, -1, -1, -1, -1, -1, -1, -1, -4, -3, -3, -2, -3, -2, -2, 0},
			{-100, -4, -4, -4, -4, -4, -4, -4, 5, 1, 1, -1, 1, -1, -1, 0},
			{-100, 1, -4, -2, -4, -2, -4, -3, 1, -1, -2, -1, -2, -1, -3, 0},
			{-100, -4, 1, -2, -4, -4, -2, -3, 1, -2, -1, -1, -2, -3, -1, 0},
			{-100, -1, -1, -1, -4, -3, -3, -2, -1, -1, -1, -1, -3, -2, -2, 0},
			{-100, -4, -4, -4, 1, -2, -2, -3, 1, -2, -2, -3, -1, -1, -1, 0},
			{-100, -1, -4, -3, -1, -1, -3, -2, -1, -1, -3, -2, -1, -1, -2, 0},
			{-100, -4, -1, -3, -1, -3, -1, -2, -1, -3, -1, -2, -1, -2, -1, 0},
			{-100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		},
		GapOpen: -10,
	}
	// nw := align.SWAffine{
	// 	Matrix:  matrix.NUC_4,
	// 	GapOpen: 0,
	// }
	for i, v := range dna {
		fmt.Println(i, v)
		aln, err := sw.Align(variable, v)
		if err == nil {
			fmt.Printf("%s\n", aln)
			fa := align.Format(variable, v, aln, '-')
			fmt.Printf("%s\n%s\n", fa[0], fa[1])
		} else {
			fmt.Println(err)
		}
		// if i > 5 {
		// 	break
		// }
	}

}
