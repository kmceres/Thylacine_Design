# Pangenome analysis
**Purpose:** Prokaryotes and some DNA viruses have variable gene content and when designing a species level molecular test it is important to identify conserved genes. The concept of a pangenome, shown below, is the idea that individuals of a species share a common set of "core" genes, but may also have differing gene content in the "accessory" genome. The combined core and accessory genes make up the pangenome of a species. 

<p align="center" width="100%">
	<img src="pangenome.png"
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
			* stringr
			* ggplot2
		* Sequence analysis:
			* seqinr
			* adegenet
			* hierfstat

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
	For example:
	```
	./scripts/run_prokka.sh testl.txt ~/Thylacine_Design/pangenome_analysis/data/fastas ~/Thylacine_Design/pangenome_analysis/data/ref/GCF_000858485.1.gbff 
	```
* Step 3: Run panaroo on annotated genomes
```
	cd gffs
	panaroo -i *.gff --clean-mode strict -t 8 --out_dir panaroo_out
```





