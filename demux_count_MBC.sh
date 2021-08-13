#!/bin/bash
# A shell script to demultiplex fastq files.
#this is written for the analysis of Maryland age samples. They were sequenced with Miseq MICRO 300 cycle (300 bp SR), 
#thus the script needs to be changed to accomodate single end reads 

datadir=/Users/asbjorh/PhD/ANK3_HTS-amplicons/Maryland/data/201124_M07166.Project_Hughes-amplicon3-2020-11-10/
cutadapt=~/.local/bin/cutadapt
#PATH=~/.local/bin/:$PATH
#export PATH

barcodes=/Users/asbjorh/PhD/ANK3_HTS-amplicons/MiSeq/demultx/test_cutadapt/N-index_sequences_ANK3-GUSB.fa
sequences=/Users/asbjorh/PhD/ANK3_HTS-amplicons/ANK3_isoform_exons.fa

## remove the sequences (+ downstream) of Illumina adapter trimming sequences. Remove 3' bases after low quality encountered.
#adapter sequences can be in any direction (3'/5' of PCR product)
for file in $(basename $datadir/*.fastq.gz |  sed 's/_L001_R1_001.fastq.gz//g' | sort -u)
do 
                cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
                -o "trimmed."$file".fastq.gz" \
                -q 10 \
                $datadir/$file"_L001_R1_001.fastq.gz"  




#Demultiplexing with barcodes. These are 4N-index, anchored to 5'. Allow for 1 errors in demultiplexing. 
#Otherwise R13 is very overrepresentated. This will lose the sequences with 1-2 N instead of 4, 
#but it is not a big nor biased problem

                cutadapt -g file:$barcodes  \
                -o $file"."{name}".fastq.gz" \
                -e 0.15 \
                -y ' round1_demultx_barcode:{name}'   --action=none \
                "trimmed."$file".fastq.gz"  



#error allowed 0.15 equals 1 error in 8-nt barcode. 



                for name in $(grep ">" $barcodes | sed 's/>//g') #get the barcode names from the sequence fasta file  
                do

#Cutadapt again to search for the isoforms of interest. Modify read name with hit.
                        cutadapt -a file:$sequences   \
                        -o $file"."$name".seq-present.fastq" \
                        -e 0.10 --revcomp \
                        -y ' present_exon:{name}'   --action=none \
                        $file"."$name".fastq.gz" 

                        #insert a double check for short unique sequences from Middle, little, Alt0. 
                        #This rids of the false positives for Little, while removing very few TP

                         grep -B 1 -A 2 -e TAAAGCCAA -e TAGATGCA -e GAAGGTGA -e GGTGATCG $file"."$name".seq-present.fastq" \
                         | grep -v "\--" >  $file"."$name".seq-present.check.fastq"

#count occurences of the different isoforms
                        for seq in $(grep ">" $sequences | sed 's/>//') "present_exon:no_adapter"
                        do
                                awk -v seq=$seq '$0~seq {n++} END {printf(FILENAME" " seq " %d\n", n)}' $file"."$name".seq-present.check.fastq" >> seq_counts.txt
                        done
                done

done


#Prepare result file for R analysis

#sed 's/present_exon:no_adapter/present-exon_no-adapter/g' seq_counts.txt > seq_counts.corrected.txt
#this is moot after double checking for presence of specific sequences. Just remove it with inv grep instead
1-ML-ANK3-pool-C1_S2.R1.seq-present.check.fastq ENSE00000997958-ENSE00001786716_Middle-Little 23
1-ML-ANK3-pool-C1_S2.R1.seq-present.check.fastq ENSE00001786716_Little 32
1-ML-ANK3-pool-C1_S2.R1.seq-present.check.fastq ENSE00000997958_Middle 492
##Awk to make table of counts
cut -d "." -f 1,2,5 seq_counts.txt |grep -v "exon:no" |\ #remove the (very few) counts of no sequence
 sed 's/ML-ANK3-pool-//g' | sed 's/C1_//g' | sed 's/.fastq//g' | \
   awk 'BEGIN { FS = "." ; OFS = "\t"} ; { print $1,$2}' | \
    awk 'BEGIN { FS = "_| " ; OFS = "\t"} ; { print $1,$2,$3,$4,$5 }' \
    > exon_counts_full.tsv


##Prepare tab file with read counts per fastq file
#fastq files are 4 lines per read. 
#Divide line number by 4, use awk to print the result next to the file name (tab separated), in stead of beneath,

for i in $(ls trimmed.R*); do echo $i; expr $(gunzip -c $i | wc -l) / 4 | awk -v i=$i '{print $0,NR%2?"\t" i :RS}' >> reads_all_fastq.tsv ; done
for i in $(ls *all-R*.R1*); do echo $i; expr $(gunzip -c $i | wc -l) / 4 | awk -v i=$i '{print $0,NR%2?"\t" i :RS}' >> reads_all_fastq.tsv ; done
for i in $(ls *all-R*.R2*); do echo $i; expr $(gunzip -c $i | wc -l) / 4 | awk -v i=$i '{print $0,NR%2?"\t" i :RS}' >> reads_all_fastq.tsv ; done
for i in $(ls *.seq-present.fastq); do echo $i; expr $(cat $i | wc -l) / 4 | awk -v i=$i '{print $0,NR%2?"\t" i :RS}' >> reads_all_fastq.tsv ; done
for i in $(ls *.seq-present.check.fastq); do echo $i; expr $(cat $i | wc -l) / 4 | awk -v i=$i '{print $0,NR%2?"\t" i :RS}' >> reads_all_fastq.tsv ; done
  
datadir=/Users/asbjorh/PhD/ANK3_HTS-amplicons/MiSeq/data/200722_M02980.Project_Hughes-amplicon2-2020-07-09/
for i in $(ls $datadir/*fastq.gz); do echo $i; expr $(gunzip -c $i | wc -l) / 4 | awk -v i=$i '{print $0,NR%2?"\t" i :RS}' >> reads_all_fastq.tsv ; done
