# RNA-Seq dcl2
dir_raw=/projects/TRIBE/RNA_seq/dcl2_rnaseq/raw
dir_mapped=/projects/TRIBE/RNA_seq/dcl2_rnaseq/mapped
rawdir_STAR=/home/sl750/STAR-2.5.3a/bin/Linux_x86_64/
# dir_genome=/home/sl750/genomes/Heinz/STAR_genome_plus_organelles
dir_genome=/data/public_data/tomato/S_lycopersicum/star_index

export PATH=$PATH:/applications/sratoolkit/sratoolkit.2.8.1-3-ubuntu64/bin
export PATH=$PATH:/applications/bedtools/bedtools2/bin/:/applications/UCSC-tools/:/applications/bowtie/bowtie-1.2.2

SRRs=(SRR6866906 SRR6866908 SRR6866909 SRR6872415 SRR6872416 SRR6872417)
#second attempt since first was only virus infected dcl2
SRRs=(SRR6872533 SRR6872532 SRR6872535)
SRR6868335 SRR6868336 SRR6868333)
SRRs=(SRR6868335)

for SRR in $SRRs; do
  print $SRR
  fastq-dump --split-files -A $SRR -O sra_fastq
done

rename 's/SRR68669/WT_/g' *
rename 's/SRR68724/dcl2b_/g' *
rename 's/SRR68683/dcl2b_/g' *
# -rwxrwxr-x 1 sm934 epigen 3.4G Aug  3 10:50 WT_06.fastq.gz
# -rwxrwxr-x 1 sm934 epigen 3.8G Aug  3 12:13 WT_08.fastq.gz
# -rwxrwxr-x 1 sm934 epigen 3.6G Aug  3 14:09 WT_09.fastq.gz
# -rwxrwxr-x 1 sm934 epigen 3.9G Aug  3 16:23 dcl2b_15.fastq.gz
# -rwxrwxr-x 1 sm934 epigen 3.6G Aug  3 18:54 dcl2b_16.fastq.gz
# -rwxrwxr-x 1 sm934 epigen 4.1G Aug  3 23:15 dcl2b_17.fastq.gz

now=$(date +"%Y_%m_%d")
date >> ${dir_mapped}/logfile.STAR.all.$now.txt

ulimit -n 10000
for file in $(ls ${dir_raw}/*_1.fastq); do
	base=$(basename ${file%_1.fastq})
	pair1=$base"_1"
	pair2=$base"_2"
	outpair=$base"_"
	
  ls -hl  ${dir_raw}/${pair1}.fastq ${dir_raw}/${pair2}.fastq  >> ${dir_mapped}/logfile.STAR.all.$now.txt

	command="STAR --runThreadN 40 --genomeDir ${dir_genome} --readFilesIn ${dir_raw}/$pair1.fastq ${dir_raw}/$pair2.fastq --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignIntronMin 20 --alignIntronMax 10000 --bamRemoveDuplicatesType UniqueIdentical --outFilterMismatchNmax 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${dir_mapped}/${base} >> ${dir_mapped}/logfile.STAR.all.$now.txt 2>&1"

	echo "$base STAR: mapping to tomato genome plus organelles"  >> ${dir_mapped}/logfile.STAR.all.$now.txt
	echo $command  >> ${dir_mapped}/logfile.STAR.all.$now.txt
	eval $command  >> ${dir_mapped}/logfile.STAR.all.$now.txt
	echo "\n *********** \n" >> /${dir_mapped}/logfile.STAR.all.$now.txt
done

gzip *.fastq

parallel 'samtools index {}' ::: *out.bam
parallel 'bamCoverage -b {} -o {.}.rpkm.bw --normalizeUsing RPKM --binSize 20' ::: *out.bam
rename 's/Aligned.sortedByCoord.out//g' *
