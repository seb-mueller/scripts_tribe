# file=/data/public_data/tomato/Solanum_pennellii_complete/Solanum_pennellii_genome/Solanum_pennellii_genome.fasta 
# genome=/data/public_data/tomato/S_lycopersicum_chromosomes.3.00.fa
# lyc=/data/public_data/tomato/S_lycopersicum/S_lycopersicum_chromosomes.3.00

# cd /data/public_data/tomato/merged_genomes
merged=/data/public_data/tomato/merged_genomes/genomes_SL30_penn_organelles_merged
PATH=$PATH:/applications/bedtools/bedtools2/bin/:/applications/UCSC-tools/:/applications/bowtie/bowtie-1.2.2


for file in *.fq.gz; do
  base=${file%_pooled_trimmed.fq.gz}
  echo $base
  bowtie --wrapper basic-0 -v 0 -k 1 -m 1 --best -p 20 -q -S $lyc $file ${base}_SolLyc_zeroMM_bowtie_unique.sam > unique.log
  samtools view -bS ${base}_SolLyc_zeroMM_bowtie_unique.sam > ${base}_SolLyc_zeroMM_bowtie_unique.unsorted.bam
  samtools sort ${base}_SolLyc_zeroMM_bowtie_unique.unsorted.bam > ${base}_SolLyc_zeroMM_bowtie_unique.bam
  samtools index ${base}_SolLyc_zeroMM_bowtie_unique.bam
done

for file in *.fq.gz; do
  base=${file%_pooled_trimmed.fq.gz}
  echo $base
  bowtie --wrapper basic-0 -v 0 -k 1 -m 1 --best -p 4 -q -S $merged $file ${base}_merged_zeroMM_bowtie_unique.sam 2>> merged_unique.log
  samtools view -bS ${base}_merged_zeroMM_bowtie_unique.sam > ${base}_merged_zeroMM_bowtie_unique.unsorted.bam
  samtools sort ${base}_merged_zeroMM_bowtie_unique.unsorted.bam > ${base}_merged_zeroMM_bowtie_unique.bam
  samtools index ${base}_merged_zeroMM_bowtie_unique.bam
done
