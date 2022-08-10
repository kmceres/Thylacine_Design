# Pangenome analysis
**Purpose:** Prokaryotes and some DNA viruses have variable gene content and when designing a species level molecular test it is important to identify conserved genes. The concept of a pangenome, shown below, is the idea that individuals of a species share a common set of "core" genes, but may also have differing gene content in the "accessory" genome. The combined core and accessory genes make up the pangenome of a species. 

<p align="center" width="100%">
	<img src="images/pangenome.png"
     alt="Pangenome concept"
     width="500"/>
</p>

## Pangenome analysis pipeline 
### Requirements:
* **Input data**:
	* fasta sequences for complete and near complete genomes 
	* High quality reference genome in genbank format is required for annotation
* **Prokka** genome annotation:
	* Not strictly required, but does reduce error in gene identification if all genomes are annotated using the same software 
	* Installation information can be found [here](https://github.com/tseemann/prokka)
* **Panaroo**: 
	* Panaroo is a software package by Tonkin-Hill et al. described [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02090-4) that was designed for prokaryotic pangenome analysis 
	* Installation and documentation can be found [here](https://gtonkinhill.github.io/panaroo/#/gettingstarted/quickstart)

* **Post processing**:
	* R scripts are available for downstream analysis 
	* Instructions for downloading R (required) and Rstudio (not required) can be found [here](https://rstudio-education.github.io/hopr/starting.html)
	* Instructions for downloading R packages can be found [here](https://r-coder.com/install-r-packages/#:~:text=Alternatively%2C%20you%20can%20install%20R,mirror%20and%20install%20the%20package.). Generally, a quick google search can help you figure out how to download a specific package if install.package("package_name") doesn't work. 
	* R packages used for pangenome analysis:
		* General purpose packages:
			* dplyr
			* tidyr
			* stringr
			* ggplot2
			* optparse
		* Sequence analysis:
			* seqinr
			* adegenet
			* hierfstat
			* rhierbaps
			* R453Plus1Toolbox


### Annotating genomes
African swine fever genomes downloaded from NCBI Virus can be found in data/fastas. The reference genome can be found in data/ref.

* Step 1: create a sequence list, for example:
```
ls *.fasta | sed "s/\.fasta//g" > seq_list.txt
```
* Step 2: Annotate genomes using prokka. 
```
chmod 755 run_prokka.sh

./run_prokka.sh seq_list.txt /path/to/fastas/ /path/to/refs
```

* For example:
```
./scripts/run_prokka.sh testl.txt ~/Thylacine_Design/pangenome_analysis/data/fastas ~/Thylacine_Design/pangenome_analysis/data/ref/GCF_000858485.1.gbff 
```

### Characterize pangenome
* Step 1: Navigate to directory containing annotation files
```
cd gffs
```
* Step 2: Run panaroo
```
panaroo -i *.gff --clean-mode strict -t 8 --out_dir panaroo_out
```
* Step 3: Use panaroo's built-in mafft aligner to generate core gene alignments
```
panaroo-msa -o . -a core --aligner mafft -t 8
```

### Identify potential target sequences
The R script process_pangenome_output.R takes the gene presence/absence matrix and the gene alignments created by panaroo and identifies potential target regions and primers. 

**Required inputs:** 
* Path to panaroo output directory 

**Optional inputs:** 
* -t target region size, default 300 bp
* -p primer size, default 20 bp
* -o output file name, defualt potential_primers.csv

**Outputs:**
* csv of potentail primer sequences, target sequence, and information about each primer sequence:
	* number of positions with SNPs
	* average heterozygosity
	* GC content
	* sequence complexity (DUST method; 100 = least complex, 0 = most complex)
* plot of locations of potential primers 

**Example usage:**

```
Rscript ./scripts/process_pangenome_output.R -f ~/CXL/ASFV/Complete_near_complete/fastas/fastas/prokka_out/panaroo_out/ -t 300 -p 20
```


### Check for off target matches
The R script compare_off_targets.R checks if potential primers will bind to off target genomes. 

**Required inputs:**
* -t path to directory containing off target whole genome sequences

**Optional inputs:**
* -f csv output from process_pangenome_output.R
* -s maximum number of substitutions allowed between primer and off target genome 

**Example usage:**
```
 Rscript ./scripts/compare_off_targets.R -f potential_primers.csv -t "~/Thylacine_Design/pangenome_analysis/data/off_target_genomes/" 
```


