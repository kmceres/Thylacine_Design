#!/bin/bash
# in get_variable_snps conda env
input=$1
path=$2
ref=$3
while IFS= read -r line
do
	mkdir prokka_out
	mkdir gffs
	prokka $path/${line}.fasta --proteins $ref --outdir prokka_out/${line}_prokka
	cp prokka_out/${line}_prokka/*.gff gffs/${line}.gff

done <"$input"