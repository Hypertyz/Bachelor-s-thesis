#!/bin/bash
# was producing an error

nextflow run nf-core/slamseq 	--input RTS0_nf_input.tsv \
								-name test_run_lk \
							 	-profile test \
								--fasta GRCm38.primary_assembly.genome.fa \
								--bed 3UTRs.bed \
								--mapping 3UTRs.bed \
						  		--igenomes_ignore=true \
								--trim5 0 \
								--multimappers=true \
							    --quantseq=true \
								--endtoend=false \
								--base_quality 20 \
								--read_length \
								--pvalue 0.05 \
								--skip_trimming=true \
								--skip_deseq2=false \
								--email lecka48@liu.se \
								--email_on_fail lecka48@liu.se
								
								