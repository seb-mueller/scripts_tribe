# base= /projects/TRIBE/srnas/tomato_external_libs/
base=/projects/TRIBE/srnas/trimmed_set3_4
dir_data=/projects/TRIBE/srnas
dir_genome=/tribe/sebastian/genomes/tomato
dir_git=~/git/resources_tomato/
# mygenome="${dir_genome}/S_lycopersicum_chromosomes.2.50.fa"
sly_mirs="${base}/mature_miRNAs_r22_UtoT_Sly.fa_clean_revcomp.fa"
# tissues=(bud flower leaf open_flower pollen root root_tip seedling unripe_fruit)
tissues=(F4s Parents)

SRRs=(SRR8763692 SRR8763691 SRR8763697 SRR8763696 SRR8763698 SRR8763699 SRR8763700 SRR8763701 SRR8763745)


SRRs=(SRR7975909 SRR7975905 SRR7975907 SRR7975929 SRR7975931 SRR7975927 SRR8607402 SRR8607401 SRR8607403 SRR9656585 SRR9656582 SRR6315778 SRR6315775 SRR6315774 SRR6315776 SRR6837350 SRR6866905 SRR6866907 SRR6866904 SRR6872537 SRR6872534 SRR6872536)

mv sra_fastq data

for SRR in $SRRs; do
  print $SRR
  fastq-dump $SRR -O data
done

# adaptors= ~/code/scripts/adapter/Small_RNA_Adapter_2.fa

export adapters=/home/sm934/code/scripts/adapter/adapter_list_8bp.txt
fastqc -t 4 -a $adapters $file

git clone git@github.com:seb-mueller/snakemake_sRNAseq.git
ln -s snakemake_sRNAseq/Snakefile .

# copy and edit both!
cp snakemake_sRNAseq/config.yaml .
cp snakemake_sRNAseq/samples.csv .

fastqc -t 4 -a $adapters root_tip_SRR8763691.fastq.gz
# trimming reads
snakemake --list
conda env list
conda activate --stack srnas
snakemake --cores 32 --keep-going -rp --show-failed-logs fq2fa
# parallel 'fastq_to_fasta -i <(zcat {}) -o {.}.fa' ::: *.gz

base=/projects/TRIBE/srnas/banana
dir_data=/projects/TRIBE/srnas
# dir_genome=/tribe/sebastian/genomes/tomato
# dir_git=~/git/resources_tomato/
mygenome="${base}/musa_acuminata_v2_pseudochromosome.fna"
# sly_mirs="${base}/mature_miRNAs_r22_UtoT_Sly.fa_clean_revcomp.fa"

cd /projects/TRIBE/srnas/tomato_external_libs/trimmed/fasta
conda activate --stack phasis
alias phasmerge="python3 ~/software/PHASIS/phasmerge"
alias phasdetect="python3 ~/software/PHASIS/phasdetect"
alias phastrigs="python3 ~/software/PHASIS/phastrigs"

# copy phasis.set into trimmed/fasta and edit manually
phasdetect
phasmerge -mode merge -dir phased_22 -debug F -pval 1
phasmerge -mode merge -dir phased_21 -debug F -pval 1

phastrigs -mode auto -dir summary_22 -mir $sly_mirs
phastrigs -mode auto -dir summary_21 -mir $sly_mirs
