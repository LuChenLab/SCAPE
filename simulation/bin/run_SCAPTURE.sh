source ~/.zshrc

conda activate SCAPTURE_env
# bamfile='/home/zhouran/data/proj/2021-0918-apa_evaluation/03.simulation_data_with_polyA/align/simu.Aligned.sortedByCoord.out.bam'
bamfile=`realpath $1`
workdir=$2

cd $2

# export PATH=/home/zhouran/data/proj/2021-0918-apa_evaluation/02.evaluation/SCAPTURE/code/SCAPTURE-main/:$PATH

/home/zhouran/data/proj/2021-0918-apa_evaluation/02.evaluation/SCAPTURE/code/SCAPTURE-main/scapture \
-m annotation -o SCAPTURE_annotation -g /mnt/raid/Ref/MusMus/release101/mm10/fasta/genome.fa \
--gtf /mnt/raid/Ref/MusMus/release101/mm10/genes/genes.gtf \
--cs /home/zhouran/data/proj/2021-0918-apa_evaluation/02.evaluation/SCAPTURE/mm10.chromsize \
--extend 2000 &> annotation.log

scapture -m PAScall \
-a SCAPTURE_annotation \
-g /mnt/raid/Ref/MusMus/release101/mm10/fasta/genome.fa -b $bamfile  \
-l 100 -o simu -p 16 --species mouse \
--polyaDB /home/zhouran/data/proj/2021-0918-apa_evaluation/02.evaluation/SCAPTURE/mm10.SupTab_KnownPASs_fourDBs.nochr.txt &> PBMC3k.PAScall.log


perl -alne '$,="\t";print @F[0..11] if $F[12] > 0 | $F[13] eq "positive";' simu.exonic.peaks.evaluated.bed simu.intronic.peaks.evaluated.bed > simu.PASquant.bed
PAS="simu.PASquant.bed"
PREFIX="simu"

bedtools intersect -a $PAS -b $PAS -s -split -f 0.3 -wo \
| perl -alne '$n1=3;$n2=15; @s1=split(/\|/,$F[$n1]);@s2=split(/\|/,$F[$n2]);  next if $s1[0] eq $s2[0]; if( $s1[0] =~ /^A\W\d+\.\d+$/ & $s2[0] !~ /^A\W\d+\.\d+$/ ){print $F[$n1];};if( $s1[0] !~ /^A\W\d+\.\d+/ & $s2[0] =~ /^A\W\d+\.\d+/ ){print $F[$n2];};  if( $s1[0] !~ /^A\W\d+\.\d+$/ & $s2[0] !~ /^A\W\d+\.\d+$/ ){  if( $s1[0] =~ /\-$s2[0]/ | $s1[0] =~ /$s2[0]\-/ ){print $F[$n1];};if( $s2[0] =~ /\-$s1[0]/ | $s2[0] =~ /$s1[0]\-/ ){print $F[$n2];}; }; ' \
| perl -alne 'if($#ARGV==0){$g{$F[0]}=1;}else{next if $g{$F[$n1]};print;}' - $PAS \
| sort -k3,3 -k1,1 -k2,2n \
| perl -alne '$,="\t"; @s=split(/\|/,$F[3]); $gene{$s[0]}++; print @F[0..2],"$s[0]-$gene{$s[0]}",@F[4..11],split(/\|/,$F[3]);' \
| sort -k1,1 -k2,2n > $PREFIX".KeepPAS.metadata"

cut -f 1-12 $PREFIX".KeepPAS.metadata" > $PREFIX".KeepPAS.bed"
cat $PREFIX".KeepPAS.bed" | perl -alne '$F[6]=$F[2];$F[7]=$F[2];$,="\t";print @F[0..11];' \
| bedToGenePred /dev/stdin /dev/stdout \
| genePredToGtf file /dev/stdin /dev/stdout \
| perl -alne '$s=$_;$s=~s/\/dev\/stdin/SCAPTURE/;print $s;' > $PREFIX".KeepPAS.gtf"

featureCounts -M -O --largestOverlap \
-F GTF -t exon -g gene_id -s 1 -T 4 -R BAM \
-a $PREFIX".KeepPAS.gtf" -o $PREFIX".KeepCell.assigned" $bamfile

# final result: simu.KeepCell.assigned
