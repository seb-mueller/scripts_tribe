# cd /projects/TRIBE/scripts/seqs
# cd /projects/TRIBE/EPRV_giri/
cd /projects/TRIBE/EPRV_annotation/EPRV_giri
seqs=eprv_P4042_+DESL.fa
seqs=ORF_ntseq_fromGIRI
seqs=LycEPRV_I_plus_subsections.fa

ref=/data/public_data/tomato/S_lycopersicum/S_lycopersicum_chromosomes.3.00.fa
refpenn=/data/public_data/tomato/Solanum_pennellii_complete/Solanum_pennellii_genome/Solanum_pennellii_genome.fasta

# minimap2 -t 10 -x asm20 -N 50 -a $ref sub.fa -o alignment2.sam
# minimap2 -t 10 -x asm20 -N 50 -a $ref $seqs -o ${seqs}.sam
# minimap2 -t 10 -x splice -N 50 -a $ref $seqs -o ${seqs}.splice.sam

# BLAT

blat $ref sub.fa output.psl
file=LycEPRV_I_plus_subsections.fa

# blat $ref $seqs ${seqs}.psl
# blat $refpenn $seqs ${seqs}_penn.psl

blat -minScore=25 -minIdentity=80 -noHead -maxIntron=1000 $ref $seqs ${seqs}_relaxed.nohead.psl

blat -minScore=25 -minIdentity=80 -noHead -maxIntron=1000 $refpenn $seqs ${seqs}_penn_relaxed.nohead.psl
