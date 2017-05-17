#!/bin/bash
#$ -l mem=3G,time=1:00:  -cwd -j y -N bowtie_build
#$ -o bowtie_build.out
#$ -e bowtie_build.err
echo "starting bowtie"
for org in "$@"
do
	bowtie2-build ${org}.fasta ${org}
done
