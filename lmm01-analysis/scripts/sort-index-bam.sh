phage=$lmm-01

# Sort and index each bamfile
cd $phage 
for file in $(ls *align.bam)
do
  samtools sort $file -o ${file}_sorted.bam
  samtools index ${file}_sorted.bam
done 

# Run bam_to_coverage
python bam_file_to_coverage.py -b . -o .

# Concatenate all files
echo "name,length,length_covered,coverage,number_bp,number_reads" > All-${phage}-coverage.csv
cat *.csv | grep gi >> All-${phage}-coverage.csv
