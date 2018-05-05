#!/bin/sh

#  blast_script.sh
#  
#
#  Created by Ava Hoffman on 11/1/16.
#
# note to user: run as 'perl blast_script.sh' in terminal

# decided not to do evalue -5 because it doesnt add that much
# length is actually included in outfmt 6 and 10 by default

~/Documents/CSU/Research/msmi_arrays/CEE_microarrays/ncbi-blast-2.5.0+/bin/blastn -query ~/Documents/CSU/Research/msmi_arrays/CEE_microarrays/final_assembled.fasta -db ~/Documents/CSU/Research/msmi_arrays/CEE_microarrays/output_andro_database -out blastn.maizeVSandro.eval-10.withalignlength.outfmt10 -evalue 1e-10 -max_target_seqs 1 -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore length"

#specifiy length criteria

~/Documents/CSU/Research/msmi_arrays/CEE_microarrays/ncbi-blast-2.5.0+/bin/blastn -query ~/Documents/CSU/Research/msmi_arrays/CEE_microarrays/final_assembled.fasta -db ~/Documents/CSU/Research/msmi_arrays/CEE_microarrays/output_sorgh_database -out blastn.maizeVSsorgh.eval-10.withalignlength.outfmt10 -evalue 1e-10 -max_target_seqs 1 -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore length"