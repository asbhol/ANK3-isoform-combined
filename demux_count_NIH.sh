#!/bin/bash
# A shell script to demultiplex fastq files, combine R1 and R2 and look for sequences of interest.
#this runs on the full set, previous versions ran on a subset of reads

datadir=/Users/asbjorh/PhD/ANK3_HTS-amplicons/MiSeq/data/200722_M02980.Project_Hughes-amplicon2-2020-07-09/
cutadapt=~/.local/bin/cutadapt
#PATH=~/.local/bin/:$PATH
#export PATH

barcodes=/Users/asbjorh/PhD/ANK3_HTS-amplicons/MiSeq/demultx/test_cutadapt/N-index_sequences_ANK3-GUSB.fa
sequences=/Users/asbjorh/PhD/ANK3_HTS-amplicons/ANK3_isoform_exons.fa

## remove the sequences (+ downstream) of Illumina adapter trimming sequences. Remove 3' bases after low quality encountered.

for file in $(basename $datadir/*.fastq.gz |  sed 's/_L001_R._001.fastq.gz//g' | sort -u)
do 
                cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
                -o "trimmed.R1."$file".fastq.gz" -p "trimmed.R2."$file".fastq.gz" \
                -q 10 \
                $datadir/$file"_L001_R1_001.fastq.gz" $datadir/$file"_L001_R2_001.fastq.gz" 




#Demultiplexing with barcodes. These are 4N-index, anchored to 5'. Allow for 1 errors in demultiplexing. 
#Otherwise R13 is very overrepresentated. This will lose the sequences with 1-2 N instead of 4, 
#but it is not a big nor biased problem

                cutadapt -g file:$barcodes  \
                -o $file"-round1-"{name}".R1.fastq.gz" \
                -p $file"-round1-"{name}".R2.fastq.gz" \
                -e 0.15 \
                -y ' round1_demultx_barcode:{name}'   --action=none \
                "trimmed.R1."$file".fastq.gz"  "trimmed.R2."$file".fastq.gz"



#error allowed 0.15 equals 1 error in 8-nt barcode. 


#Then all the read pairs in which no barcode could be found will end up in round1-unknown.R1.fastq.gz and round1-unknown.R2.fastq.gz.
#This will also include the pairs in which the barcode was not actually in R1, but in R2. To demultiplex these reads as well, run Cutadapt a second time with
#those “unknown” files as input, but also reverse the roles of R1 and R2


                cutadapt -g file:$barcodes    \
                -o $file"-round2-"{name}".R1.fastq.gz" \
                -p $file"-round2-"{name}".R2.fastq.gz" \
                -e 0.15 \
                -y ' round2_demultx_barcode:{name}'   --action=none \
                $file-round1-unknown.R2.fastq.gz $file-round1-unknown.R1.fastq.gz

#Flash to combine R1 and R2, after round1 and round2 is combined with cat

                for name in $(grep ">" $barcodes | sed 's/>//g')
                do
                        cat $file"-round1-"$name".R1.fastq.gz" $file"-round2-"$name".R1.fastq.gz" > $file"-all-"$name".R1.fastq.gz"
                        cat $file"-round1-"$name".R2.fastq.gz" $file"-round2-"$name".R2.fastq.gz" > $file"-all-"$name".R2.fastq.gz"

                        flash2 --max-overlap 300 -z -o $file"-all-"$name $file-all-$name.R1.fastq.gz $file-all-$name.R2.fastq.gz

#Cutadapt again to search for the isoforms of interest. Modify read name with hit.
                        cutadapt -a file:$sequences   \
                        -o $file"-all-"$name".extendedFrags.seq-present.fastq" \
                        -e 0.10 --revcomp \
                        -y ' present_exon:{name}'   --action=none \
                        $file"-all-"$name".extendedFrags.fastq.gz"

                        #insert a double check for short unique sequences from Middle, little, Alt0. 
                        #This rids of the false positives for Little, while removing very few TP

                         grep -B 1 -A 2 -e TAAAGCCAA -e TAGATGCA -e GAAGGTGA -e GGTGATCG $file"-all-"$name".extendedFrags.seq-present.fastq" \
                         | grep -v "\--" >  $file"-all-"$name".extendedFrags.seq-present.check.fastq"

#count occurences of the different isoforms
                        for seq in $(grep ">" $sequences | sed 's/>//') "present_exon:no_adapter"
                        do
                                awk -v seq=$seq '$0~seq {n++} END {printf(FILENAME" " seq " %d\n", n)}' $file"-all-"$name".extendedFrags.seq-present.check.fastq" >> seq_counts.txt
                        done
                done

done


#Prepare result file for R analysis

#sed 's/present_exon:no_adapter/present-exon_no-adapter/g' seq_counts.txt > seq_counts.corrected.txt
#this is moot after double checking for presence of specific sequences. Just remove it with inv grep instead

##Awk to make table of counts
cut -d . -f 1,5 seq_counts.txt | sed 's/AH-ANK3-pool-//g' \
| cut -d "-" -f 1,2,4,5,6 | sed 's/.fastq//g' \
| awk 'BEGIN { FS = "_| " ; OFS = "\t"} ; { print $1,$2,$3,$4,$5 }' > exon_counts_full.tsv


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
