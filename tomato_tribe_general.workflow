cat /data/public_data/tomato/habrochaites/CBYS010000001-CBYS010042990.fasta > Sol_Hab_plus_penn.fa
echo -e "\n" >> Sol_Hab_plus_penn.fa

cat /data/public_data/tomato/Solanum_pennellii_complete/Solanum_pennellii_genome/Solanum_pennellii_genome.fasta >> Sol_Hab_plus_penn.fa

/applications/ncbi-blast-2.2.29+/bin/makeblastdb -in Sol_Hab_plus_penn.fa -out Shabpenn -dbtype 'nucl'
#hen1_vs_Shab.out

for file in $(ls *.fasta); do
	/applications/ncbi-blast-2.2.29+/bin/blastn -db Shabpenn -query $file -out ${file%.fasta}_vs_Shabpenn.txt
done

#Hen1
samtools faidx Sol_Hab_plus_penn.fa "ENA|HG975441|HG975441.1":43134000-43149000 | revseq -sequence stdin -outseq stdout > hen1_penn.fa

samtools faidx Sol_Hab_plus_penn.fa "ENA|CBYS010010298|CBYS010010298.1":24000-41000 > hen1_hab.fa

cat hen1_hab.fa hen1_penn.fa hen1_XM_004233836.fasta HEN1_S_Lycopersicum_gDNA.fa | muscle -clw -out hen1_muscle.aln

#met1
samtools faidx Sol_Hab_plus_penn.fa "ENA|HG975450|HG975450.1":20814476-20821185 | revseq -sequence stdin -outseq stdout > met1_penn.fa

samtools faidx Sol_Hab_plus_penn.fa "ENA|CBYS010021371|CBYS010021371.1":9500-16102 > met1_hab.fa 

seqret -sbegin1 9000 -sequence MET1_S_Lycopersicum_gDNA.fa -outseq MET1_S_Lycopersicum_gDNA_trimmed.fa

cat MET1_S_Lycopersicum_gDNA_trimmed.fa hen1_XM_004233836.fasta met1_hab.fa met1_penn.fa | muscle -clw -out met1_muscle.aln

#DRM2
samtools faidx Sol_Hab_plus_penn.fa "ENA|CBYS010011811|CBYS010011811.1":1-15000 > drm2_hab.fa
samtools faidx Sol_Hab_plus_penn.fa "ENA|HG975441|HG975441.1":37217193-37228598 > drm2_penn.fa

cat DRM2_NM_001246974.fasta DRM5_DRM2homologue_S_Lycopersicum_gDNA.fa drm2_hab.fa drm2_penn.fa | muscle -clw -out drm2_muscle.aln
Slyref_SPennqryCTr.vcf
#Pol4
samtools faidx Sol_Hab_plus_penn.fa "ENA|HG975447|HG975447.1":68071558-68085174 | revseq -sequence stdin -outseq stdout > pol4_penn.fa

samtools faidx Sol_Hab_plus_penn.fa "ENA|CBYS010011132|CBYS010011132.1":48000-61000 | revseq -sequence stdin -outseq stdout > pol4_hab.fa

cat Pol4_XM_004245609.fasta pol4_hab.fa pol4_penn.fa POLIV_S_Lycopersicum_gDNA.fa | muscle -clw -out pol4_muscle.aln

##DRM2 like



######
###compute repeat track
####computing repeat-stretch track
library(Biostrings)
library(rtracklayer)
genometomato <- readDNAStringSet("/data/public_data/tomato/assembly_build_2.50/S_lycopersicum_chromosomes.2.50.fa")
genometomatomod <- genometomato
genometomatomod <- chartr("N","T",genometomatomod)


genometomatomodrc <- reverseComplement(genometomatomod)
genometomatomodc <- complement(genometomatomod)
genometomatomodr <- reverse(genometomatomod)

pb <- txtProgressBar(1,length(names(genometomato))); count <- 1
clist <- list()
for (chr in names(genometomato)) {
	starts <- (1:(length(genometomatomod[[chr]])-21))
	seqsdict <- PDict(Views(genometomatomod[[chr]],start=starts,end=starts+20))
	clist[[chr]] <- rowSums(vcountPDict(seqsdict,genometomatomod))
	setTxtProgressBar(pb, count <- count + 1)
}
close(pb)
export.bw(RleList(clist),"Sly25repeats20mers.bw")

#reverse complement
clist <- list()
for (chr in names(genometomato)) {
	starts <- (1:(length(genometomatomod[[chr]])-21))
	seqsdict <- PDict(Views(genometomatomod[[chr]],start=starts,end=starts+20))
	clist[[chr]] <- rowSums(vcountPDict(seqsdict,genometomatomodrc))
}
export.bw(RleList(clist),"Sly25repeats20mers_revcomp.bw")

## computing homologie with Arabidopsis
genomeat <- readDNAStringSet("/data/public_data/arabidopsis/TAIR_9/arabidopsis_genome")
genomeatrc <- reverseComplement(genomeat)
# pb <- txtProgressBar(1,length(names(genometomato))); count <- 1
# clist <- list()
# for (chr in names(genometomato)) {
# 	starts <- (1:(length(genometomatomod[[chr]])-21))
# 	seqsdict <- PDict(Views(genometomatomod[[chr]],start=starts,end=starts+20))
# 	clist[[chr]] <- countPDict(seqsdict,genomeat[[1]])
# 	for (scaf in 2:length(genomeat)) {
# 		genomeatsub <- genomeat[[scaf]]
# 		clist[[chr]] <- clist[[chr]] + countPDict(seqsdict,genomeatsub)
# 	}
# 		setTxtProgressBar(pb, count <- count + 1)
# }

pb <- txtProgressBar(1,length(names(genometomato))); count <- 1
clist <- list()
for (chr in names(genometomato)) {
	starts <- (1:(length(genometomatomod[[chr]])-21))
	seqsdict <- PDict(Views(genometomatomod[[chr]],start=starts,end=starts+20))
	clist[[chr]] <- rowSums(vcountPDict(seqsdict,genomeat))
	clist[[chr]] <- clist[[chr]] + rowSums(vcountPDict(seqsdict,genomeatrc))
	setTxtProgressBar(pb, count <- count + 1)
}
close(pb)
export.bw(RleList(clist),"Sly25vsAth20mers.bw")

#cd /home/bioinf/sm934/analysis/tomato/custom_tracks
##penelli
genomepenn <- readDNAStringSet("/data/public_data/tomato/Solanum_pennellii_complete/Solanum_pennellii_genome/Solanum_pennellii_genome.fasta")
genomepennrc <- reverseComplement(genomepenn)

pb <- txtProgressBar(1,length(names(genometomato))); count <- 1
clist <- list()
for (chr in names(genometomato)) {
	starts <- (1:(length(genometomatomod[[chr]])-21))
	seqsdict <- PDict(Views(genometomatomod[[chr]],start=starts,end=starts+20))
	clist[[chr]] <- rowSums(vcountPDict(seqsdict,genomepenn))
	clist[[chr]] <- clist[[chr]] + rowSums(vcountPDict(seqsdict,genomepennrc))
	setTxtProgressBar(pb, count <- count + 1)
}
close(pb)
export.bw(RleList(clist),"Sly25vsSpenn20mers.bw")

##habrochaites
genome <- readsly $shab
$mummer/show-snps -Clr Slyref_SHabqrDNAStringSet("/data/public_data/tomato/habrochaites/CBYS010000001-CBYS010042990.fasta")

genometmp <- do.call(paste0,genome)
genome <- DNAStringSet(genometmp)

genomerc <- reverseComplement(genome)

pb <- txtProgressBar(1,length(names(genometomato))); count <- 1
clist <- list()
for (chr in names(genometomato)) {
	starts <- (1:(length(genometomatomod[[chr]])-21))
	seqsdict <- PDict(Views(genometomatomod[[chr]],start=starts,end=starts+20))
	clist[[chr]] <- rowSums(vcountPDict(seqsdict,genome))
	clist[[chr]] <- clist[[chr]] + rowSums(vcountPDict(seqsdict,genomerc))
	setTxtProgressBar(pb, count <- count + 1)
}
close(pb)
export.bw(RleList(clist),"Sly25vsShab20mers.bw")


##potato
genomepot <- readDNAStringSet("/data/public_data/potato/v403/PGSC_DM_v4.03_pseudomolecules.fasta")
genomepotrc <- reverseComplement(genomepot)

pb <- txtProgressBar(1,length(names(genometomato))); count <- 1
clist <- list()
for (chr in names(genometomato)) {
	starts <- (1:(length(genometomatomod[[chr]])-21))
	seqsdict <- PDict(Views(genometomatomod[[chr]],start=starts,end=starts+20))
	clist[[chr]] <- rowSums(vcountPDict(seqsdict,genomepot))
	clist[[chr]] <- clist[[chr]] + rowSums(vcountPDict(seqsdict,genomepotrc))
	setTxtProgressBar(pb, count <- count + 1)
}
close(pb)
export.bw(RleList(clist),"Sly25vsPotato20mers.bw")


#####
##Finding SNPs between M82 and penellii
# make CPPFLAGS="-O3 -DSIXTYFOURBITS"
#
sly=/data/public_data/tomato/assembly_build_2.50/S_lycopersicum_chromosomes.2.50.fa
sly30=/data/public_data/tomato/S_lycopersicum_chromosomes.3.00.fa
spenn=/data/public_data/tomato/Solanum_pennellii_complete/Solanum_pennellii_genome/Solanum_pennellii_genome.fasta
shab=/data/public_data/tomato/habrochaites/CBYS010000001-CBYS010042990.fasta
spimp=/data/public_data/tomato/pimpinellifolium/Spimpinellifolium_genome.contigs.fasta
sm82=/data/public_data/tomato/Solanum_lycopersicum_M82_genome/M82_genome/M82_genome.fasta
sI=/home/syngenta_ftp/files/nrgene/I_Tomato/S_lycopersicum_I_pseudo_chromosomes_plus_organelles.fa
sC=/home/syngenta_ftp/files/nrgene/C_Tomato/S_habrochaites_C_pseudo_chromosomes_plus_organelles.fa


mummer=/applications/MUMmer/MUMmer3.23/
#identify repeats
$mummer/nucmer --maxmatch --nosimplify --prefix=Sly_repeats $sly $sly
## was killed after 1 week of running..

#Lyc vs Penn
$mummer/nucmer --prefix=Slyref_qry $sly $spenn
$mummer/show-snps -Clr Slyref_qry.delta > Slyref_qry.snps
#-T 	Switch to tab-delimited format
#-I			Do not output indels
#-r			Sort SNPs by reference position
$mummer/show-snps -CTr Slyref_qry.delta > Slyref_SPennqryCTr.snps


#/applications/MUMmer/MUMmer3.23/show-coords -rcl Slyref_qry.delta > Slyref_qry.coords
#/applications/MUMmer/MUMmer3.23/show-aligns Slyref_qry.delta Sly Spenn > ref_qry.aligns
#mummerplot Slyref_qry.delta $sly $spenn -Q qry.fasta --filter --layout

##Lyc vs habrochaites
$mummer/nucmer --prefix=Slyref_qry $sly $shab
$mummer/show-snps -Clr Slyref_SHabqry.delta > Slyref_SHabqry.snps
#-T 	Switch to tab-delimited format
$mummer/show-snps -CTr Slyref_qry.delta > Slyref_SHabqryCTr.snps

##ussing maxmatch to loosen uniqueness criteria.
$mummer/nucmer --maxmatch --prefix=Slyref_maxmatch_qry $sly $shab
##after 1 week of high mem computing: ERROR: mummer and/or mgaps returned non-zero

/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_HabENA_Ly250 $shab $sly"
cat Solanum_HabENA_Ly250.delta | tr \| _ > tmp;cp tmp Solanum_HabENA_Ly250.delta

/scripts/csmit -m 300G -c 1  "$mummer/nucmer --mincluster 100 -maxgap 300 --prefix=Slyref_qryc100g300 $sly $shab"

#$mummer/show-snps -Clr Slyref_maxmatch_SHabqry.delta > Slyref_maxmatch_SHabqry.snps
#-T 	Switch to tab-delimited format
#-I 	Do not report indels
$mummer/show-snps -CTr Slyref_maxmatch_qry.delta > Slyref_maxmatch_SHabqryCTr.snps


$mummer/delta-filter -r -q Slyref_qry.delta > Slyref_qry.filter.delta
$mummer/delta-filter -u 10 -r -q Slyref_qry.delta > Slyref_qry.filter.u10.delta

#coords
base=Slyref_qry.filter.u10.delta
$mummer/show-coords -rcl $base > $base.coords
cat $base.coords | tail -n +6 | perl -lane '{print "$F[17]\tNucmer_delta_filter\thomolog\t$F[0]\t$F[1]\t\.\t\.\t\.\tID=@F"}' > $base.coords.gff3



$mummer/show-coords -rcl Slyref_qry.filter.delta > Slyref_qry.filter.delta.coords
$mummer/show-snps -CTr Slyref_qry.delta > Slyref_SHabqryCTr.filter.snps
$mummer/show-snps -Tr Slyref_qry.delta > Slyref_SHabqryCTr.filter.repeats.snps
#june15
$mummer/show-snps -CTr Slyref_qry.filter.delta > Slyref_SHabqryCTr.filter.filter.snps
#converting coords in gff:



##Lyc vs pimpinellifolium
$mummer/nucmer --prefix=Slyref_qry $sly $spimp
#$mummer/show-snps -Clr Slyref_SPimpqry.delta > Slyref_SPimpqry.snps
#-T 	Switch to tab-delimited format
$mummer/show-snps -CTr Slyref_qry.delta > Slyref_SPimpqryCTr.snps

##Lyc vs M82
$mummer/nucmer --prefix=Slyref_qry $sly $sm82
#$mummer/show-snps -Clr Slyref_SPimpqry.delta > Slyref_M82qry.snps
#-T 	Switch to tab-delimited format
$mummer/show-snps -CTr Slyref_qry.delta > Slyref_M82qryCTr.snps

#on node7
cd /data/public_data/tomato/assembly_build_3.00/inhouse_annotation
##Lyc vs Penn v2
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Ly250_Penn $sly $spenn"
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Penn_Ly250 $spenn $sly"

##Heinz30 vs Penn v2
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Heinz30_Penn $sly30 $spenn"
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Penn_Heinz30 $spenn $sly30"

##Heinz30 vs M82 
/scripts/csmit -b -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Heinz30_M82 $sly30 $sm82"
/scripts/csmit -b -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_M82_Heinz30 $sm82 $sly30"

##Heinz30 vs C
/scripts/csmit -b -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Heinz30_C $sly30 $sC"
/scripts/csmit -b -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_C_Heinz30 $sC $sly30"

#Heinz30 vs Heinz25,Shab,I
/scripts/csmit -b -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Heinz30_Heinz25 $sly30 $sly"
/scripts/csmit -b -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Heinz30_SHab $sly30 $shab"
/scripts/csmit -b -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Heinz30_I $sly30 $sI"
##C vs I
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_C_I $sC $sI"
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_I_C $sI $sC"


##Heinz vs I
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Ly250_I $sly $sI"
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_I_Ly250 $sI $sly"

##Heinz vs C
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Ly250_C $sly $sC"
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_C_Ly250 $sC $sly"

##M82 vs Penn
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_M82_Penn $sm82 $spenn"
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Penn_M82 $spenn $sm82"


##Heinz vs M82
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_Ly250_M82 $sly $sm82"
/scripts/csmit -m 300G -c 1  "$mummer/nucmer --prefix=Solanum_M82_Ly250 $sm82 $sly"



#delta filter
for file in $(ls *[^r].delta); do
	echo ${file%.delta}.filter.delta
	if [ ! -e "${file%.delta}.filter.delta" ]
	then $mummer/delta-filter -r -q $file > ${file%.delta}.filter.delta &
	fi
done

#creating coord, snps and its counterparts gff and vcf:
for file in $(ls Solanum*[^r].filter.delta); do
# for file in $(ls Solanum*M82*.filter.delta); do
	echo $file
	$mummer/show-coords -rcl $file > $file.coords
	cat $file.coords | tail -n +6 | perl -lane '{print "$F[17]\tNucmer_delta_filter\thomolog\t$F[0]\t$F[1]\t\.\t\.\t\.\tID=@F"}' > $file.coords.gff3
	$mummer/show-snps -CTr $file > $file.snps
	cat header.vcf > ${file}.vcf
	cat ${file}.snps | tail -n +5 | perl -lane 'if ($F[1]=~/[ACGT]/ & $F[2]=~/[ACGT]/) {print "$F[8]\t$F[0]\t\.\t$F[1]\t$F[2]\t10\tPASS\tAF=1"}' >> ${file}.vcf
	/applications/IGVTools/v2_3_53/igvtools index ${file}.vcf
done


snpbase=Slyref_SHabqryCTr.filter
snpbase=Slyref_SHabqryCTr.filter.filter
cat header.vcf > ${snpbase}.vcf
cat ${snpbase}.snps | tail -n +5 | perl -lane 'if ($F[1]=~/[ACGT]/ & $F[2]=~/[ACGT]/) {print "$F[8]\t$F[0]\t\.\t$F[1]\t$F[2]\t10\tPASS\tAF=1"}' >> ${snpbase}.vcf
/home/sm934/data/software/IGVTools/igvtools index ${snpbase}.vcf


snpbase=Slyref_SHabqryCTr.filter.repeats
##has 2 more collums in snps file and tends to not be sorted!!
cat header.vcf > ${snpbase}.vcf
cat ${snpbase}.snps | tail -n +5 | perl -lane 'if ($F[1]=~/[ACGT]/ & $F[2]=~/[ACGT]/) {print "$F[10]\t$F[0]\t\.\t$F[1]\t$F[2]\t10\tPASS\tAF=1"}' | sort -rb -k1n,1 -k2n,2 >> ${snpbase}.vcf

#cat Slyref_SHabqryCTr.vcf | perl -lane 'print $F[0]' | uniq > Slyref_SHabqryCTr_chromosomes.txt
#/applications/samtools-0.1.19/bcftools/bcftools view -Sb -D Slyref_SHabqryCTr_chromosomes.txt Slyref_SHabqryCTr.vcf > Slyref_SHabqryCTr.bcf

file="Slyref_SPennqryCTr"
cat header.vcf > $file.vcf
cat $file.snps | tail -n +5 | perl -lane 'if ($F[1]=~/[ACGT]/ & $F[2]=~/[ACGT]/) {print "$F[8]\t$F[0]\t\.\t$F[1]\t$F[2]\t10\tPASS\tAF=1"}' >> $file.vcf

file="Slyref_SPimpqryCTr"
cat header.vcf > $file.vcf
cat $file.snps | tail -n +5 | perl -lane 'if ($F[1]=~/[ACGT]/ & $F[2]=~/[ACGT]/) {print "$F[8]\t$F[0]\t\.\t$F[1]\t$F[2]\t10\tPASS\tAF=1"}' >> $file.vcf

file="Slyref_M82qryCTr"
cat header.vcf > $file.vcf
cat $file.snps | tail -n +5 | perl -lane 'if ($F[1]=~/[ACGT]/ & $F[2]=~/[ACGT]/) {print "$F[8]\t$F[0]\t\.\t$F[1]\t$F[2]\t10\tPASS\tAF=1"}' >> $file.vcf


##injecting snps in genome
#using picard tools to create dict file for genome
java -jar /home/sm934/hydro/CreateSequenceDictionary.jar R=S_lycopersicum_chromosomes.2.50.fa O= S_lycopersicum_chromosomes.2.50.dict

##creater empty vcf as dummy
head -n 18 Slyref_SHabqryCTr.vcf > dummy.vcf

##use GATK to inject SNPs from vcf
java -jar /home/sm934/data/software/GenomeAnalysisTK.jar -R /home/sm934/igv/genomes/S_lycopersicum_chromosomes.2.50.fa -o S_lycopersicum_chromosomes_SNPfromSHamasked.2.50.fa --snpmask Slyref_SHabqryCTr.vcf -T FastaAlternateReferenceMaker --variant dummy.vcf

##same for pennellii
java -jar /applications/GenomeAnalysisTK/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -R /data/public_data/tomato/assembly_build_2.50/S_lycopersicum_chromosomes.2.50.fa -o S_lycopersicum_chromosomes250_SNPfromSpenninjected.fa --snpmask Slyref_SPennqryCTr.vcf

##chromsome names have been anoyingly been changed, have to rename them:
##
cat S_lycopersicum_chromosomes_SNPfromSHamasked.2.50.fa | perl -plane 's/>(\d+)/$1-1/e' | perl -plane 's/(\d{2})/>SL2.50ch$1/' | perl -plane 's/^(\d{1})/>SL2.50ch0$1/' > tmp.fa
mv tmp.fa S_lycopersicum_chromosomes_SNPfromSHamasked.2.50.fa

##again for repeat SNPs
java -jar /home/sm934/data/software/GenomeAnalysisTK.jar -R /home/sm934/igv/genomes/S_lycopersicum_chromosomes.2.50.fa -o S_lycopersicum_chromosomes_SNPfromSHamasked.repeats.2.50.fa --snpmask Slyref_SHabqryCTr.filter.repeats.vcf -T FastaAlternateReferenceMaker --variant dummy.vcf

##chromsome names have been anoyingly been changed, have to rename them:
##
cat S_lycopersicum_chromosomes_SNPfromSHamasked.2.50.fa | perl -plane 's/>(\d+)/$1-1/e' | perl -plane 's/(\d{2})/>SL2.50ch$1/' | perl -plane 's/^(\d{1})/>SL2.50ch0$1/' > tmp.fa
mv tmp.fa S_lycopersicum_chromosomes_SNPfromSHamasked.2.50.fa



##compare results for sanity check
samtools faidx S_lycopersicum_chromosomes_SNPfromSHamasked.2.50.fa "SL2.50ch00:358-370"
samtools faidx $sly "SL2.50ch00:358-370"

##add organelles
#/data/public_data/tomato/assembly_build_2.50/S_lycopersicum_chromosomes.2.50.SNPinjected

cat /home/sm934/analysis/tomato/mummer/S_lycopersicum_chromosomes_SNPfromSHamasked.2.50.fa /data/public_data/tomato/mitochondrion/mitochondria_Solyc_merged.fa /data/public_data/tomato/chloroplast/solanum_lycopersicum_chloroplast_AM087200_simpleheader.fasta > Sly250ShabSnps.plus_organelles.fa

perl ~/bin/fasta2chrLengths.pl Sly250ShabSnps.plus_organelles.fa > tomato250.plus_organelles.chrom.sizes.txt


##Computing Nucleosome Occupancy
library(NuPoP)
library(rtracklayer)
chrsizesSly25 <- read.table("/data/public_data/tomato/assembly_build_2.50/tomato250.chrom.sizes.txt")
perl ~/bin/splitfasta.pl /data/public_data/tomato/assembly_build_2.50/S_lycopersicum_chromosomes.2.50.fa
predNuPoP("/home/sm934/analysis/zhengming/NucPos/SL2.50ch00",species=0,model=4)
results=readNuPoP("/home/bioinf/sm934/analysis/zhengming/NucPos/SL2.50ch00_Prediction1.txt",startPos=1,endPos=chrsizesSly25[chrsizesSly25$V1=="SL2.50ch00",2])
export.bw(RleList(SL2.50ch00=results$Affinity), "NucAffinitySly250model4.bw")
export.bw(RleList(SL2.50ch00=results$Occup), "NucOccupSly250model4.bw")

##mrsFast
#-n max mappings per read
mates1=Plantd_Sly140718_I126_FCC4RWRACXX_L4_index14_1.fq
mates2=Plantd_Sly140718_I126_FCC4RWRACXX_L4_index14_2.fq
mrsfast=/applications/mrsFAST/mrsfast3.3.3/mrsfast
genome=/data/public_data/tomato/assembly_build_2.50/S_lycopersicum_chromosomes.2.50.fa
$mrsfast --search $genome --pe --seq1 $mates1 --seq2 $mates2 --threads 10 -n 100 --min 100 --max 500 -e 3 --max-discordant-cutoff 50 


#####investigating amiguity mapping:#############

sm934@node4:~/test$ bowtie2 -f -x genome --score-min L,-0,-0 --n-ceil L,3,0.150 --np 0 -U seq_1mmN.fa
1 reads; of these:
  1 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    1 (100.00%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
100.00% overall alignment rate
@HD     VN:1.0  SO:unsorted
@SQ     SN:c1   LN:43
@PG     ID:bowtie2      PN:bowtie2      VN:2.2.3        CL:"/applications/bowtie2/bowtie2-2.2.3/bowtie2-align-s --wrapper basic-0 -f -x genome --score-min L,-0,-0 --n-ceil L,3,0.150 --np 0 -U seq_1mmN.fa"
cel-miR-35-3p   0       c1      4       42      22M     *       0       0       TNACCGGGTGGAAACTAGCAGT  IIIIIIIIIIIIIIIIIIIIII  AS:i:0  XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:1C20       YT:Z:UU
sm934@node4:~/test$ bowtie2 -f -x genome --score-min L,-0,-0 --n-ceil L,3,0.150 --np 0 -U seq_1mm.fa
1 reads; of these:
  1 (100.00%) were unpaired; of these:
    1 (100.00%) aligned 0 times
    0 (0.00%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
0.00% overall alignment rate
@HD     VN:1.0  SO:unsorted
@SQ     SN:c1   LN:43
@PG     ID:bowtie2      PN:bowtie2      VN:2.2.3        CL:"/applications/bowtie2/bowtie2-2.2.3/bowtie2-align-s --wrapper basic-0 -f -x genome --score-min L,-0,-0 --n-ceil L,3,0.150 --np 0 -U seq_1mm.fa"
cel-miR-35-3p   4       *       0       0       *       *       0       0       TGACCGGGTGGAAACTAGCAGT  IIIIIIIIIIIIIIIIIIIIII  YT:Z:UU
sm934@node4:~/test$ bowtie2 -f -x genome --score-min L,-0,-0 --n-ceil L,3,0.150 --np 0 -U seq.fa
1 reads; of these:
  1 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    1 (100.00%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
100.00% overall alignment rate
@HD     VN:1.0  SO:unsorted
@SQ     SN:c1   LN:43
@PG     ID:bowtie2      PN:bowtie2      VN:2.2.3        CL:"/applications/bowtie2/bowtie2-2.2.3/bowtie2-align-s --wrapper basic-0 -f -x genome --score-min L,-0,-0 --n-ceil L,3,0.150 --np 0 -U seq.fa"
cel-miR-35-3p   0       c1      4       42      22M     *       0       0       TCACCGGGTGGAAACTAGCAGT  IIIIIIIIIIIIIIIIIIIIII  AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:22 YT:Z:UU
sm934@node4:~/test$ bowtie2 -f -x genome --score-min L,-0,-0 --n-ceil L,3,0.150 --np 0 -U seq.fa

###calling hight quality from Input chipseq data
using samtools pileup ..
#http://samtools.sourceforge.net/mpileup.shtml
#cd /home/sm934/analysis/syngenta/rawdata/SNPcalling_from_shortreads
export PATH=$PATH:/applications/samtools/samtools-1.1/:/applications/bcftools/bcftools-1.1

bamfile="../INPUT-alpha_PlantA_GenotypeSha_Lane1_Cycles16_140628_I168_FCC4U5HACXX_L4_Index7.RmDup.trimmed.k100stringentASonSly250.bam"
bamfile="../INPUT-alpha_Plante_GenotypeSly_Lane2_Cycles14_140704_I883_FCC4RNVACXX_L3_Index5.RmDup.trimmed.k100stringentASonSly250.bam"
genome="/data/public_data/tomato/assembly_build_2.50/S_lycopersicum_chromosomes.2.50.fa"
/applications/samtools/samtools-1.1/samtools mpileup -g -f $genome $bamfile > my-raw.bcf
#bcftools view -bvcg my-raw.bcf > my-var2.bcf
bcftools call -vmA -O v my-raw.bcf > my-var_vmA.vcf 
bcftools filter --include 'TYPE="snp" && QUAL>=5' -O v my-var_vmA.bcf > my-var_vmA_filter.vcf


###mapping real data to snp injected sly genome
#using chipseq data for testing (plant d), only the first pair no simplify results

/applications/bowtie2/bowtie2-2.2.3/bowtie2-build S_lycopersicum_chromosomes_SNPfromSHamasked.2.50.fa Sly250ShabSnps


data=/home/sm934/analysis/syngenta/rawdata/INPUT-alpha_Plantd_GenotypeSly_Lane2_Cycles16_140704_I883_FCC4RNVACXX_L3_Index16_1.RmDup.trimmed.fq.gz
cd /home/sm934/analysis/tomato/mummer/mappingtest
bowtiepath=/applications/bowtie2/bowtie2-2.2.3/
index=/data/public_data/tomato/assembly_build_2.50/S_lycopersicum_chromosomes.2.50
indexsnps=/home/sm934/analysis/tomato/mummer/Sly250ShabSnps
sampath=/applications/samtools-0.1.19/


echo "k10L02stdgenomenp1.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${index} -U $data -S k10L02stdgenome.sam -p 14 >> log.txt 2>&1
echo "k10L02stdgenomenp6.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${index} -U $data -S k10L02stdgenome_np6.sam -p 14 --n-ceil L,0,0.150 --np 6  >> log.txt 2>&1

echo "k10L02SNPgenomenp1.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${indexsnps} -U $data -S k10L02SNPgenome_np1.sam -p 14 --n-ceil L,0,0.150 --np 1 >> log.txt 2>&1
echo "k10L02SNPgenomenp6.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${indexsnps} -U $data -S k10L02SNPgenome_np6.sam -p 14 --n-ceil L,0,0.150 --np 6 >> log.txt 2>&1
echo "k10L02SNPgenomenp0.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${indexsnps} -U $data -S k10L02SNPgenome_np0.sam -p 14 --n-ceil L,0,0.150 --np 0  >> log.txt 2>&1
echo "k10L02SNPgenomenp0_10Ns.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${indexsnps} -U $data -S k10L02SNPgenome_np0_10Ns.sam -p 14 --n-ceil L,10,0 --np 0  >> log.txt 2>&1

echo "k10L02SNPgenomenp0_10NsL0.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${indexsnps} -U $data -S k10L02SNPgenome_np0_10NsL0.sam -p 14 --n-ceil L,0,0 --np 0  >> log.txt 2>&1



data=/home/sm934/analysis/syngenta/rawdata/INPUT-alpha_Plantd_GenotypeSly_Lane2_Cycles16_140704_I883_FCC4RNVACXX_L3_Index16_1.RmDup.trimmed.fq.gz
data2=/home/sm934/analysis/syngenta/rawdata/INPUT-alpha_Plantd_GenotypeSly_Lane2_Cycles16_140704_I883_FCC4RNVACXX_L3_Index16_2.RmDup.trimmed.fq.gz


echo "k10L02stdgenomenp0.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${index} -U $data -S k10L02stdgenome_np0.sam -p 14 --n-ceil L,0,0.150 --np 0  >> log.txt 2>&1
echo "k10L02stdgenomenp6.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${index} -U $data -S k10L02stdgenome_np6.sam -p 14 --n-ceil L,0,0.150 --np 6  >> log.txt 2>&1

echo "k10L02SNPgenomenp0.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${indexsnps} -U $data -S k10L02SNPgenome_np0.sam -p 14 --n-ceil L,0,0.150 --np 0  >> log.txt 2>&1
echo "k10L02SNPgenomenp6.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${indexsnps} -U $data -S k10L02SNPgenome_np6.sam -p 14 --n-ceil L,0,0.150 --np 6 >> log.txt 2>&1

#Shab
data=/home/sm934/analysis/syngenta/rawdata/INPUT-alpha_PlantA_GenotypeSha_Lane1_Cycles16_140628_I168_FCC4U5HACXX_L4_Index7_1.fq.gz
echo "k10L02stdSHAnp0.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${index} -U $data -S k10L02stdSHA_np0.sam -p 14 --n-ceil L,0,0.150 --np 0  >> log.txt 2>&1
echo "k10L02stdSHAnp6.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${index} -U $data -S k10L02stdSHA_np6.sam -p 14 --n-ceil L,0,0.150 --np 6  >> log.txt 2>&1

echo "k10L02SNPSHAnp0.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${indexsnps} -U $data -S k10L02SNPSHA_np0.sam -p 14 --n-ceil L,0,0.150 --np 0  >> log.txt 2>&1
echo "k10L02SNPSHAnp6.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 10 --score-min L,-0.2,-0.2 -x ${indexsnps} -U $data -S k10L02SNPSHA_np6.sam -p 14 --n-ceil L,0,0.150 --np 6 >> log.txt 2>&1

for file in $(ls *.sam); do
	base=${file%.sam}
 	$sampath/samtools view -bS $base.sam > $base"_unsorted".bam
 	rm $base.sam
 	$sampath/samtools sort $base"_unsorted".bam $base
 	rm $base"_unsorted".bam
 	$sampath/samtools index $base.bam
done

###
--
--n-ceil <func> max # of Ns.
--np [1] Sets penalty
--score-min <func> L,-0.2,-0.2=-20.2 ~= 3mismatches (1 == score of -6)

###another run for PE using K27 since plant e and A have similar # of reads of rm dup.
data=/home/sm934/analysis/syngenta/rawdata/K27_Plante_GenotypeSly_Lane2_Cycles16_140704_I883_FCC4RNVACXX_L3_Index2_1.RmDup.trimmed.fq.gz
data2=/home/sm934/analysis/syngenta/rawdata/K27_Plante_GenotypeSly_Lane2_Cycles16_140704_I883_FCC4RNVACXX_L3_Index2_2.RmDup.trimmed.fq.gz

echo "k10L02stdgenomenp0.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 100 --score-min L,-0.4,-0.4 -x ${index} -1 $data -2 $data2 -S k10L02stdgenome_np0.sam -p 14 --n-ceil L,0,0.150 --np 0.1  >> log.txt 2>&1
echo "k10L02stdgenomenp6.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 100 --score-min L,-0.4,-0.4 -x ${index} -1 $data -2 $data2 -S k10L02stdgenome_np6.sam -p 14 --n-ceil L,0,0.150 --np 6  >> log.txt 2>&1

echo "k10L02SNPgenomenp0.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 100 --score-min L,-0.4,-0.4 -x ${indexsnps} -1 $data -2 $data2 -S k10L02SNPgenome_np0.sam -p 14 --n-ceil L,0,0.150 --np 0.1  >> log.txt 2>&1
echo "k10L02SNPgenomenp6.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 100 --score-min L,-0.4,-0.4 -x ${indexsnps} -1 $data -2 $data2 -S k10L02SNPgenome_np6.sam -p 14 --n-ceil L,0,0.150 --np 6 >> log.txt 2>&1

#Shab
data=/home/sm934/analysis/syngenta/rawdata/K27_PlantA_GenotypeSha_Lane1_Cycles16_140628_I168_FCC4U5HACXX_L4_Index5_1.RmDup.trimmed.fq.gz
data2=/home/sm934/analysis/syngenta/rawdata/K27_PlantA_GenotypeSha_Lane1_Cycles16_140628_I168_FCC4U5HACXX_L4_Index5_2.RmDup.trimmed.fq.gz


echo "k10L02stdSHAnp0.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 100 --score-min L,-0.4,-0.4 -x ${index} -1 $data -2 $data2 -S k10L02stdSHA_np0.sam -p 14 --n-ceil L,0,0.150 --np 0.1  >> log.txt 2>&1
echo "k10L02stdSHAnp6.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 100 --score-min L,-0.4,-0.4 -x ${index} -1 $data -2 $data2 -S k10L02stdSHA_np6.sam -p 14 --n-ceil L,0,0.150 --np 6  >> log.txt 2>&1

echo "k10L02SNPSHAnp0.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 100 --score-min L,-0.4,-0.4 -x ${indexsnps} -1 $data -2 $data2 -S k10L02SNPSHA_np0.sam -p 14 --n-ceil L,0,0.150 --np 0.1  >> log.txt 2>&1
echo "k10L02SNPSHAnp6.sam:" >> log.txt
${bowtiepath}/bowtie2 -q -k 100 --score-min L,-0.4,-0.4 -x ${indexsnps} -1 $data -2 $data2 -S k10L02SNPSHA_np6.sam -p 14 --n-ceil L,0,0.150 --np 6 >> log.txt 2>&1

##dez14: SNP aware mapping on chipseq data with bowtie2 since tophat failed:
#/home/sm934/analysis/syngenta/rawdata_rnaseq/snptest
#PlantA

 index=/data/public_data/tomato/assembly_build_2.50/S_lycopersicum_chromosomes.2.50
 indexsnps=/home/sm934/analysis/tomato/mummer/Sly250ShabSnps
 indexhab=/data/public_data/tomato/habrochaites/SHab
 bowtiepath=/applications/bowtie2/bowtie2-2.2.3/
 sampath=/applications/samtools-0.1.19/
 btpath=/applications/bedtools-2.17.0/bin/
 chromsizes=/data/public_data/tomato/assembly_build_2.50/tomato250.chrom.sizes.txt
 threads=40
 k=10

file=PlantA_SHa140718_I126_FCC4RWRACXX_L4_index13_1.fq.gz
base=${file%_1.fq.gz}
pair1=$base"_1"
pair2=$base"_2"
 
$bowtiepath/bowtie2 -q -k $k --score-min L,-0.2,-0.2 -x $indexsnps -1 ../$pair1.RmDup.trimmed.fq.gz -2 ../$pair2.RmDup.trimmed.fq.gz -S $base"_unsorted".sam -p $threads  >> log.txt 2>&1
#50%

$bowtiepath/bowtie2 -q -k $k --score-min L,-0.2,-0.2 -x $indexsnps -1 ../$pair1.RmDup.trimmed.fq.gz -2 ../$pair2.RmDup.trimmed.fq.gz -S $base"_unsorted_snpaware".sam -p $threads --np 0  >> log.txt 2>&1
#51%

$bowtiepath/bowtie2 -q -k $k --score-min L,-0.2,-0.2 -x $indexsnps -1 ../$pair1.RmDup.trimmed.fq.gz -2 ../$pair2.RmDup.trimmed.fq.gz -S $base"_unsorted_np6".sam -p $threads --np 6  >> log.txt 2>&1
#44%

##normal genome (without snps)
$bowtiepath/bowtie2 -q -k $k --score-min L,-0.2,-0.2 -x $index -1 ../$pair1.RmDup.trimmed.fq.gz -2 ../$pair2.RmDup.trimmed.fq.gz -S $base"_unsorted_np6".sam -p $threads --np 6  >> log.txt 2>&1
#45%

echo "with seed mm -N 1" >> log.txt
$bowtiepath/bowtie2 -q -k $k --score-min L,-0.2,-0.2 -x $indexsnps -1 ../$pair1.RmDup.trimmed.fq.gz -2 ../$pair2.RmDup.trimmed.fq.gz -S $base"_unsorted_snpaware".sam -p $threads --np 0 -N 1  >> log.txt 2>&1
#-N 1 increase of mapping by ~2% -> 51%

echo "map against hab genome" >> log.txt
$bowtiepath/bowtie2 -q -k $k --score-min L,-0.2,-0.2 -x $indexhab -1 ../$pair1.RmDup.trimmed.fq.gz -2 ../$pair2.RmDup.trimmed.fq.gz -S $base"_unsorted_snpaware".sam -p $threads --np 0 -N 1  >> log.txt 2>&1
##60%

echo "map against Sly genome with N 1 np0" >> log.txt
$bowtiepath/bowtie2 -q -k $k --score-min L,-0.2,-0.2 -x $index -1 ../$pair1.RmDup.trimmed.fq.gz -2 ../$pair2.RmDup.trimmed.fq.gz -S $base"_unsorted_snpaware".sam -p $threads --np 0 -N 1  >> log.txt 2>&1
##45%
##conclusion for stringent map comparison of Shab sample: 45% Sly; 51% SlySNPs; 60%Shab

## is difference as big if less stringent bowtie? e.g. default
echo "map against Sly genome with N 1 np0, relaxed" >> log.txt
$bowtiepath/bowtie2 -q -k $k -x $index -1 ../$pair1.RmDup.trimmed.fq.gz -2 ../$pair2.RmDup.trimmed.fq.gz -S $base"_unsorted_snpaware".sam -p $threads --np 0 -N 1  >> log.txt 2>&1
##82.3

echo "map against SlySNP genome with N 1 np0, relaxed" >> log.txt
$bowtiepath/bowtie2 -q -k $k -x $indexsnps -1 ../$pair1.RmDup.trimmed.fq.gz -2 ../$pair2.RmDup.trimmed.fq.gz -S $base"_unsorted_snpaware".sam -p $threads --np 0 -N 1  >> log.txt 2>&1
#83.5%

echo "map against Hab genome with N 1 np0, relaxed" >> log.txt
$bowtiepath/bowtie2 -q -k $k -x $indexhab -1 ../$pair1.RmDup.trimmed.fq.gz -2 ../$pair2.RmDup.trimmed.fq.gz -S $base"_unsorted_snpaware".sam -p $threads --np 0 -N 1  >> log.txt 2>&1
#87.3

##trying tophat
tophatpath=/applications/tophat/tophat-2.0.12.Linux_x86_64/
bowtiepath=/applications/bowtie2/bowtie2-2.2.3/
sampath=/applications/samtools-0.1.19/
index=/data/public_data/tomato/assembly_build_2.50/S_lycopersicum_chromosomes.2.50
indexsnps=/home/sm934/analysis/tomato/mummer/Sly250ShabSnps
chromsizes=/data/public_data/tomato/assembly_build_2.50/tomato250.chrom.sizes.txt
threads=20
k=1
mismatches=20

PATH=$PATH:$sampath:$bowtiepath


	${tophatpath}tophat -p ${threads} -o tophat --mate-inner-dist 20 --mate-std-dev 100 --read-mismatches ${mismatches} --read-edit-dist ${mismatches} --min-intron-length 40 --max-intron-length 10000 --max-multihits ${k} ${index} ../$pair1.RmDup.trimmed.fq.gz ../$pair2.RmDup.trimmed.fq.gz
	
	${tophatpath)tophat -p ${threads) -o tophatsnp  --mate-inner-dist 20 --mate-std-dev 100 --read-mismatches ${mismatches) --read-edit-dist ${mismatches) --min-intron-length 40 --max-intron-length 10000 --max-multihits ${k) ${indexsnps) ../$pair1.RmDup.trimmed.fq.gz ../$pair2.RmDup.trimmed.fq.gz

	${tophatpath)tophat -p ${threads) -o tophatsnp2 --b2-np 0 --b2-N 1 --b2-L 10 --mate-inner-dist 20 --mate-std-dev 100 --read-mismatches ${mismatches) --read-edit-dist ${mismatches) --min-intron-length 40 --max-intron-length 10000 --max-multihits ${k) ${indexsnps) ../$pair1.RmDup.trimmed.fq.gz ../$pair2.RmDup.trimmed.fq.gz
	
	tail -n 1000000 PlantA_SHa140718_I126_FCC4RWRACXX_L4_index13_1.RmDup.trimmed.fq > pair1.fq
	tail -n 1000000 PlantA_SHa140718_I126_FCC4RWRACXX_L4_index13_2.RmDup.trimmed.fq > pair2.fq
	
	${tophatpath}tophat --b2-np 0 -p ${threads} -o tophatsnp --mate-inner-dist 20 --mate-std-dev 100 --read-mismatches ${mismatches} --read-edit-dist ${mismatches} --min-intron-length 40 --max-intron-length 10000 --max-multihits ${k} ${indexsnps} pair1.fq pair2.fq
	
	${tophatpath}tophat --b2-np 0 -p ${threads} -o tophat --mate-inner-dist 20 --mate-std-dev 100 --read-mismatches ${mismatches} --read-edit-dist ${mismatches} --min-intron-length 40 --max-intron-length 10000 --max-multihits ${k} ${index} pair1.fq pair2.fq
	
	${tophatpath}tophat --b2-np 0 --b2-score-min L,-0.4,-0.4 -p ${threads} -o tophatsnpsminscore --mate-inner-dist 20 --mate-std-dev 100 --read-mismatches ${mismatches} --read-edit-dist ${mismatches} --min-intron-length 40 --max-intron-length 10000 --max-multihits ${k} ${indexsnps} pair1.fq pair2.fq
	
##conclusion
## -b2-score-min L,-0.4,-0.4 needs to be set since tophat seem to ignore default and goes by number of mismatches. Since Ns are counted as mismatches (even thought not penalized in alignment score AS), one needs to set mismatches hight (e.g. 20) and using AS as cutoff by setting score-min. 

##renaming NCBI annotation chromosomes
head -n 9 ref_SL2.50_top_level.gff3 > header.gff3
less ref_SL2.50_top_level.gff3 | grep "NC_015444.2" | perl -pe 's/NC_015444.2/SL2.50ch07/g' >> ref_SL2.50_top_level.rnamed.chr7.gff3 

sed -e 's/^/s%/' -e 's/=/%/' -e 's/$/%g/' mappingfile | sed -f - ref_SL2.50_gnomon_top_level.gff3 > ref_SL2.50_gnomon_top_level.renamed.gff3
