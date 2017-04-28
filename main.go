package main

import (
	"flag"
	"fmt"
	"sync"

	"github.com/biogo/biogo/seq"

	"github.com/BurntSushi/toml"
	"github.com/jgcarvalho/PhageAnalysis/pda"
	"github.com/jgcarvalho/PhageAnalysis/pep"
	"github.com/jgcarvalho/PhageAnalysis/prot"
)

type Config struct {
	Title      string
	Translate  bool   `toml:"translate"`
	Mapping    bool   `toml:"mapping"`
	Nproc      int    `toml:"nproc"`
	PepLen     int    `toml:"peptide-length"`
	ProtFasta  string `toml:"protein-fasta"`
	PepFasta   string `toml:"peptide-fasta"`
	ForwPrimer string `toml:"forward-primer"`
	RevPrimer  string `toml:"reverse-primer"`

	Nrandom    int
	ExpScore   float64 `toml:"exp-score"`
	QuadWindow int     `toml:"quadrant-window"`
}

func main() {
	showCredits()
	var config Config
	var cf string
	// cf = "/home/jgcarvalho/gocode/src/github.com/jgcarvalho/PhageAnalysis/config.toml"
	flag.StringVar(&cf, "c", "", "Configuration file (TOML format)")
	flag.Parse()

	c, err := toml.DecodeFile(cf, &config)
	if err != nil {
		fmt.Println("Arquivo de configuração inválido!", err)
		return
	}
	if len(c.Undecoded()) > 0 {
		fmt.Printf("Chaves desconhecidas no arquivo de configuração: %q\n", c.Undecoded())
		return
	}

	//runtime.GOMAXPROCS(config.Nproc)

	var peptides, unreliable []pep.Peptide
	// Le as sequencias de dna ou aa dos peptideos. Se for dna é necessario traduzir
	if config.Translate {
		//tamanho do peptideo
		pepNBases := (config.PepLen) * 3
		var dna []seq.Sequence
		// fp := config.ForwPrimer
		// rp := config.RevPrimer
		fmt.Println("Creating template to translate peptides")
		template := pda.CreateTemplate(config.ForwPrimer, config.RevPrimer, pepNBases)
		fmt.Println(template)

		// files := []string{fdna}
		fmt.Println("Reading peptides sequencing data")
		dna, err = pda.ReadMFastQ(config.PepFasta)
		if err != nil {
			fmt.Println("ERRO", err)
		}

		// Peptideos obtidos do sequenciamento. Tem que respeitar o padrao NNK e não possuir
		// stop codons na sequencia
		fmt.Println("Translating peptides and computing frequency")
		peptides, unreliable = pda.GetPeptides(dna, template)
		pep.Peptides(peptides).SaveFasta("peptides.fasta")
		pep.Peptides(unreliable).SaveFasta("peptides_unreliable.fasta")

	} else {
		peptides, err = pda.ReadTranslatedPeptides(config.PepFasta, config.PepLen)
	}

	fmt.Println("peptides", peptides)
	for _, p := range peptides {
		fmt.Println(p)
	}

	if config.Mapping {
		// Peptideos Randomicos com o mesmo numero da biblioteca e com a mesma proporcao de
		// aminoacidos
		fmt.Println("Generating random peptides")
		// var RandLibraries [][]pep.Peptides
		RandLibraries := make([][]pep.Peptide, config.Nrandom)
		for i := 0; i < config.Nrandom; i++ {
			RandLibraries[i] = pep.RandomLibrary(peptides)
		}

		// Le as proteinas
		fmt.Println("Reading proteins")
		proteins, err := pda.ReadMProteins(config.ProtFasta)
		if err != nil {
			fmt.Println("ERRO", err)
		}

		var wg1 sync.WaitGroup
		for j := 0; j < config.Nproc; j++ {
			wg1.Add(1)
			go func(p *[]prot.Protein, pep *[]pep.Peptide, rlibpep *[][]pep.Peptide, exp float64, j int) {
				defer wg1.Done()
				for i := (j * len(proteins)) / config.Nproc; i < ((j+1)*len(proteins))/config.Nproc; i++ {
					(*p)[i].Analysis(*pep, *rlibpep, exp)
				}
			}(&proteins, &peptides, &RandLibraries, config.ExpScore, j)
		}

		wg1.Wait()
		prot.Proteins(proteins).SaveRank("rank_proteins")
	}

}

func showCredits() {
	credits :=
		`#############################################################################
# PhageFinger version=0.1                                                   #
# Brazilian Biosciences National Laboratory                                 #
# Author: José Geraldo de Carvalho Pereira                                  #
#############################################################################`

	fmt.Printf("%s\n\n", credits)
}
