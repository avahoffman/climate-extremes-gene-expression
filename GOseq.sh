#!/bin/bash

#SBATCH --time=10000:00:00
#SBATCH --job-name=goseq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=ERRORLOGjan7
#SBATCH --output=OUTLOGjan7

#  GOseq.sh
#  
#
#  Created by Ava Hoffman on 12/13/16.
#

#http://trinotate.github.io/ for ALL downloads. There are many..

export PATH=$PATH:~/scratch_dir/ncbi-blast-2.3.0+
export PATH=$PATH:~/scratch_dir/ncbi-blast-2.3.0+/bin
export PATH=$PATH:~/scratch_dir/Trinotate-3.0.1
export PATH=$PATH:~/scratch_dir/TransDecoder-2.1.0
export PATH=$PATH:~/scratch_dir/hmmer/binaries
#export PATH=$PATH:~/scratch_dir/blastdbs
export PATH=$PATH:~/scratch_dir/signalp-4.1
#export PATH=$PATH:~/scratch_dir/chlorop-1.1
export PATH=$PATH:~/scratch_dir/sqlite
#export PATH=$PATH:~/scratch_dir/tmhmm-2.0c/bin
#export PATH=$PATH:~/scratch_dir/trinityrnaseq-2.1.1/


#~/scratch_dir/ncbi-blast-2.3.0+/bin/makeblastdb -in uniprot_sprot.pep -dbtype prot
#~/scratch_dir/hmmer/src/hmmpress Pfam-A.hmm

#generate pep file, "likely proteins" - takes about a minute
#~/scratch_dir/TransDecoder-2.1.0/TransDecoder.LongOrfs -t ~/scratch_dir/microarray/final_maize_genes_to_blast.fasta
#~/scratch_dir/TransDecoder-2.1.0/TransDecoder.Predict -t ~/scratch_dir/microarray/final_maize_genes_to_blast.fasta

#search "Trinity" transcripts - this part takes 9 minutes
#blastx -query ~/scratch_dir/microarray/final_maize_genes_to_blast.fasta -db ~/scratch_dir/microarray/uniprot_sprot.pep -num_threads 40 -max_target_seqs 1 -outfmt 6 > blastx.MICROARRAY.trinotate.outfmt6

#search Transdecoder-predicted proteins - almost 15 mins
#blastp -query ~/scratch_dir/microarray/final_maize_genes_to_blast.fasta.transdecoder.pep -db ~/scratch_dir/microarray/uniprot_sprot.pep -num_threads 40 -max_target_seqs 1 -outfmt 6 > blastp.MICROARRAY.trinotate.outfmt6

#HMMER identifies protein domains - almost 15 mins
#~/scratch_dir/hmmer/binaries/hmmscan --cpu 16 --domtblout TrinotatePFAMmicroarray.out ~/scratch_dir/microarray/Pfam-A.hmm ~/scratch_dir/microarray/final_maize_genes_to_blast.fasta.transdecoder.pep > pfam.log

#SignalP predicts signal peptide cleavage sites --- not doing
#~/scratch_dir/signalp-4.1/signalp -f short -n signalp.microarray.out ~/scratch_dir/microarray/final_maize_genes_to_blast.fasta.transdecoder.pep

#predict transmembrane helices in proteins ----- not doing
#~/scratch_dir/tmhmm-2.0c/bin/tmhmm --short < ~/scratch_dir/andro/Trinotate/Trinity_andro_out.Trinity.02122016.fasta.transdecoder.pep > tmhmm.andro.out

#predict chloroplast stuff _--- not doing
#~/scratch_dir/chlorop-1.1/chlorop -F ~/scratch_dir/andro/Trinotate/Trinity_andro_out.Trinity.02122016.fasta.transdecoder.pep > chlorop.andro.out

#run Trinotate!!!
#~/scratch_dir/trinityrnaseq-2.1.1/util/support_scripts/get_Trinity_gene_to_trans_map.pl ~/scratch_dir/microarray/final_maize_genes_to_blast.fasta > microarray.gene_trans_map
#~/scratch_dir/Trinotate-3.0.1/Trinotate ~/scratch_dir/Trinotate.sqlite init --gene_trans_map ~/scratch_dir/microarray/microarray.gene_trans_map --transcript_fasta ~/scratch_dir/microarray/final_maize_genes_to_blast.fasta --transdecoder_pep ~/scratch_dir/microarray/final_maize_genes_to_blast.fasta.transdecoder.pep
#~/scratch_dir/Trinotate-3.0.1/Trinotate ~/scratch_dir/Trinotate.sqlite LOAD_swissprot_blastp ~/scratch_dir/microarray/blastp.MICROARRAY.trinotate.outfmt6
#~/scratch_dir/Trinotate-3.0.1/Trinotate ~/scratch_dir/Trinotate.sqlite LOAD_swissprot_blastx ~/scratch_dir/microarray/blastx.MICROARRAY.trinotate.outfmt6
#~/scratch_dir/Trinotate-3.0.1/Trinotate ~/scratch_dir/Trinotate.sqlite LOAD_pfam ~/scratch_dir/microarray/TrinotatePFAMmicroarray.out
#~/scratch_dir/Trinotate-3.0.1/Trinotate ~/scratch_dir/Trinotate.sqlite report > trinotate_annotation_report_MICROARRAY.xls

#extract GO assignments per gene
#~/scratch_dir/Trinotate-3.0.1/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls ~/scratch_dir/microarray/trinotate_annotation_report_MICROARRAY.xls -T --include_ancestral_terms > go_annotations.txt

#run assignment using GOseq. First, highly expressed on average across all samples

~/scratch_dir/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling factor_labs.txt --GO_assignments go_annotations.txt --lengths gene-lengths-sppcontrast.txt
