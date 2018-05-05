# with open("/Users/avahoffman/Documents/CSU/Research/msmi_arrays/CEE_microarrays/final_assembled_maize.fasta") as readfile, open("/Users/avahoffman/Documents/CSU/Research/msmi_arrays/CEE_microarrays/maize_gene_sequences.fasta",'w') as outfile:
# 	linecount=1
# 	for line in readfile:
# 		linecount+=1
# 		if linecount%2==0:
# 			outfile.write(line)

with open("/Users/avahoffman/Documents/CSU/Research/msmi_arrays/CEE_microarrays/8-annotations/maize_genes_to_annotate.txt") as readfile, open("/Users/avahoffman/Documents/CSU/Research/msmi_arrays/CEE_microarrays/8-annotations/final_maize_genes_to_blast.fasta",'w') as outfile:
    genenumber=1
    for line in readfile:
        #values=[]
        newlist=line.split('\t')
        for index, value in enumerate(newlist):
            if index == 0:
                print(value)
                outfile.write(">TRINITY_")
                outfile.write(value)
                outfile.write("_c0_g"+str(genenumber)+"_i1")
                outfile.write("\n")
                genenumber+=1
            if index == 1:
                start = 0
                for nucleotide in value:
                    start+=1
                    outfile.write(nucleotide)
                    if start > 60:
                        outfile.write("\n")
                        start=0
