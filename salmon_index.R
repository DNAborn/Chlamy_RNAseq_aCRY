
# get Chlamydomonas genomic data
# https://data.jgi.doe.gov/refine-download/phytozome?q=Chlamydomonas+reinhardtii+CC-4532+v6.1
# (Create Phytozome Account)

## Files
# Genome (DNA):            CreinhardtiiCC_4532_707_v6.0.hardmasked.fa.gz
# Transcriptome (RNA):     CreinhardtiiCC_4532_707_v6.1.transcript.fa.gz
# Gene Info:               CreinhardtiiCC_4532_707_v6.1.gene.gff3.gz
# Master Annotation: 	     CreinhardtiiCC_4532_707_v6.1.master_annotation_table.tsv

genomicdir=""
genome=""
transcripts=""
genes=""
annotation=""

# create list of chromosomes
cd $genomicdir
grep "^>" <(gunzip -c $genome) | cut -d " " -f 1 > decoys_chlamy_v6.1.txt

# to modify, copy to linux partition
mkdir ~/decoys
cp decoys_chlamy_v6.1.txt ~/decoys

# remove ">"
sed -i.bak -e 's/>//g' ~/decoys/decoys_chlamy_v6.1.txt

# Optional: check txt file
vim ~/decoys/decoys_chlamy_v6.1.txt
# to exit: -press 'ESC', type ":q!", press 'ENTER'

# copy back
cp ~/decoys/decoys_chlamy_v6.1.txt $genomicdir

# combine cdna ncdna and dna in one file
cd $genomicdir

cat $genome $transcripts > gentrome_chlamy_v6.1.fa.gz

# make index with salmon
# is salmon activated?
conda activate salmon
cd $gdir
indexname=chlamy_v6.1_index
echo $indexname

salmon index -t gentrome_chlamy_v6.1.fa.gz -d decoys_chlamy_v6.1.txt -p 10 -i $indexname
