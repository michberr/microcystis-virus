virome=Sample_49639
cd $virome
pwd
echo -e "[`date`]\t Running bowtie2 on virome $virome"

for genome in lmm01 MaMV; do    
    R1=$(ls *_R1_*.fastq.gz)
    R2=$(ls *_R2_*.fastq.gz)
    echo $R1
    echo $R2
#    bowtie2 --no-unal --non-deterministic -p 20 -X 2000 -x ../${genome}/${genome} -1 $R1 -2 $R2 | samtools view -Sb - > ${virome}-${genome}-align.bam
done

echo -e "[`date`]\t Finished with bowtie alignment to microcystis phages"

