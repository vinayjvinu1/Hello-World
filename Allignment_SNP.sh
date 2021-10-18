bwa index FRG1.fasta
bwa mem -t 4 FRG1.fasta 1T-Old_R1.fastq.gz 1T-Old_R2.fastq.gz > 3T-New.sam

samtools view -S -b 1T-Old.sam > 1T-Old.bam

samtools sort -@ 4 1T-Old.bam -o 1T-Old.sorted.bam


samtools faidx FRG1.fasta

samtools mpileup -g -f FRG1.fasta 1T-Old.sorted.bam > 1T-Old.raw.bcf

bcftools call -O b -vc 1T-Old.raw.bcf > 1T-Old.var.bcf

bcftools view 1T-Old.var.bcf | vcfutils.pl varFilter -> 1T-Old.var-final.vcf
