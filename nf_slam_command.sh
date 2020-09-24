#!/bin/bash
# Run as follows: ./nf_slam_command.sh [run names]

nextflow run nf-core/slamseq 	--input RTS0_nf_input.tsv \
								-name $1 \
							 	-profile docker,test \
								--fasta GRCm38.primary_assembly.genome.fa \
								--bed 3UTRs.bed \
								--mapping 3UTRs.bed \
								--trim5 2 \
								--polyA 8 \
								--multimappers false \
							    --quantseq true \
								--endtoend true \
								--base_quality 20 \
								--pvalue 0.05 \
								--readlength 300 \
								--skip_deseq2 true \
								--email lecka48@liu.se \
								--email_on_fail lecka48@liu.se
								
								