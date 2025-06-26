module load samtools

for i in fastq/*.fq.gz; do
        i=$(basename $i)
        i=${i%.fq.gz}
        echo "file $i"
        bwa-mem2 mem -c 1000000 -L 100 -t 56 -a ../bwa-index-mm39 fastq/$i.fq.gz | samtools view -bS - > bam/$i.bam
done