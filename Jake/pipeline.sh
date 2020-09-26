# Trim read1
./Tools/bbmap/bbduk.sh -Xmx24g in=./RTS0/Data_RTS0/first_trimtagger-pooled-library-2S_S1_L001_R1_001.fastq.gz out=./first_trim/trimmed_R1.fastq.gz ftl=15

# Trim read2
./Tools/bbmap/bbduk.sh -Xmx24g in=./RTS0/Data_RTS0/first_trimtagger-pooled-library-2S_S1_L001_R2_001.fastq.gz out=./first_trim/trimmed_R2.fastq.gz minlen=6 ftl=2 ftr=7 mbq=10

# Repair
./Tools/bbmap/repair.sh -Xmx24g in=./first_trim/trimmed_R1.fastq.gz in2=./first_trim/trimmed_R2.fastq.gz out=repaired/R1_repaired.fastq.gz out2=repaired/R2_repaired.fastq.gz

# Demultiplex
tagdust -t 16 -fe 1 -arch arch1.txt -o demux/RTS0 repaired/R1_repaired.fastq.gz repaired/R2_repaired.fastq.gz

# Nextflow
nextflow run nf-core/slamseq --input RTS0_nf_input.tsv --genome GRCm38 --max_memory '32.GB' --quantseq true
