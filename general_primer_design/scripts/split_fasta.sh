input=$1
out=fastas

seqkit split2 --by-size 1 $input -O $out

find $out -name "*.fasta" \
    | while read f; do \
        mv $f $out/$(seqkit seq --name --only-id $f).fasta; \
    done