# https://www-ncbi-nlm-nih-gov.ezp.lib.cam.ac.uk/Traces/study/?acc=SRP127908&o=acc_s%3Aa&s=SRR7637630,SRR7637631,SRR7637632,SRR6436049,SRR6436050,SRR6436051

lyc=/data/public_data/tomato/S_lycopersicum/S_lycopersicum_chromosomes.3.00
export PATH=$PATH:/applications/sratoolkit/sratoolkit.2.8.1-3-ubuntu64/bin
export PATH=$PATH:/applications/bedtools/bedtools2/bin/:/applications/UCSC-tools/:/applications/bowtie/bowtie-1.2.2

SRRs=(SRR6436051 SRR6436050 SRR6436048 SRR6436055)
SRRs=(SRR7637630 SRR7637631 SRR7637632 SRR6436049)

for SRR in $SRRs; do
  print $SRR
  fastq-dump -A $SRR -O sra_fastq
done

mv SRR6436051.fastq.gz wt-1.fastq.gz 
mv SRR6436050.fastq.gz wt-2.fastq.gz 
mv SRR6436048.fastq.gz dclab-1.fastq.gz 
mv SRR6436055.fastq.gz dcllab-2.fastq.gz

rename 's/SRR6436049/wt-3/g' *
rename 's/SRR7637630/wt-1_tmv/g' *
rename 's/SRR7637631/wt-2_tmv/g' *
rename 's/SRR7637632/wt-3_tmv/g' *

source activate srna_mapping
for file in *.fastq.gz; do

conda activate srnas
adaptor=AGATCGGAAGAGCAC
file=wt-1.fastq.gz

for file in *.fastq.gz; do
cutadapt --cores=8 -m 15 -M 40  -e 0.1 -q 20 -O 1 -a $adaptor  -o ${file%%.fastq.gz}_trimmed.fastq $file > ${file%%fastq}_trimmed.report
done

gzip *.fastq


for file in *trimmed.fq.gz; do
  base=${file%_trimmed.fastq.gz}
  echo $base
  bowtie --wrapper basic-0 -v 0 -k 1 --best -p 10 -q -S $lyc $file ${base}_SolLyc_zeroMM_bowtie_multi.sam > multi.log
  samtools view -bS ${base}_SolLyc_zeroMM_bowtie_multi.sam > ${base}_SolLyc_zeroMM_bowtie_multi.unsorted.bam
  samtools sort ${base}_SolLyc_zeroMM_bowtie_multi.unsorted.bam > ${base}_SolLyc_zeroMM_bowtie_multi.bam
  samtools index ${base}_SolLyc_zeroMM_bowtie_multi.bam
done

for file in *.fq.gz; do
  base=${file%_pooled_trimmed.fq.gz}
  echo $base
  bowtie --wrapper basic-0 -v 0 -k 1 -m 1 --best -p 4 -q -S $merged $file ${base}_merged_zeroMM_bowtie_unique.sam 2>> merged_unique.log
  samtools view -bS ${base}_merged_zeroMM_bowtie_unique.sam > ${base}_merged_zeroMM_bowtie_unique.unsorted.bam
  samtools sort ${base}_merged_zeroMM_bowtie_unique.unsorted.bam > ${base}_merged_zeroMM_bowtie_unique.bam
  samtools index ${base}_merged_zeroMM_bowtie_unique.bam
done

PATH=$PATH:/applications/bedtools/bedtools2/bin/:/applications/UCSC-tools/
chromsize=/data/public_data/tomato/S_lycopersicum/S_lycopersicum_chromosomes.3.00_chrlens.txt
chromsize=/data/public_data/tomato/merged_genomes/genomes_SL30_penn_organelles_merged.chrsize

# file=22_KP4042-2_SolLyc_zeroMM_bowtie.sorted.bam
for file in *multi.bam; do
  echo $file
  # for l in $(seq 21 24); do
  for l in 21 22 24; do
    echo $l;
    length=${l};
    echo $length;
    samtools view -h $file | awk -v l=$length -F"\t" ' { if (/^@/ || $6 == l"M") {print} } ' | samtools view -bS - | genomeCoverageBed -bg -ibam stdin -g $chromsize > ${file%.bam}$l.bedgraph
    # samtools view -h $file | awk -v l=$length -F"\t" ' { if (/^@/ || $6 == l) {print} } ' >  tmp
    bedGraphToBigWig ${file%.bam}$l.bedgraph $chromsize ${file%.bam}$l.bw
  done
done
