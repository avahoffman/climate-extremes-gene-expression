'''
splitting up files
'''
# import os
# directory="/Users/avahoffman/Documents/CSU/Research/msmi_arrays/CEE_microarrays"
# for filename in os.listdir(directory):
# 	if filename.endswith(".gpr"): 
# 		print(os.path.join(directory, filename))
# 		with open(directory+"/"+filename) as open_file, open(directory+"/"+filename+"_qualtrimmed.txt","w") as new_file:
# 			linecount=0
# 			for line in open_file:
# 				linecount+=1
# 				newlist=line.split('\t')
# 				if linecount > 31:
# 					new_file.write(line)
# 		continue
# 	else:
# 		continue
'''
reformat
'''
# import os
# directory="/Users/avahoffman/Documents/CSU/Research/msmi_arrays/CEE_microarrays"
# for filename in os.listdir(directory):
# 	if filename.endswith(".fasta"): 
# 		print(os.path.join(directory, filename))
# 		with open(directory+"/"+filename) as open_file, open(directory+"/"+filename+"_fastafixed.fasta","w") as new_file:
# 			genecount=0
# 			for line in open_file:
# 				if '>' in line:
# 					genecount+=1
# 					new_file.write("\n"+line)
# 				else:
# 					new_file.write(line.strip())
# 			print(genecount)
# 		continue
# 	else:
# 		continue
'''
concatenate files
'''
import os
directory="/Users/avahoffman/Documents/CSU/Research/msmi_arrays/CEE_microarrays"
with open(directory+"/"+"final_assembled.fasta","w") as new_file:
	genecount=0
	for filename in os.listdir(directory):
		if filename.endswith(".fasta"): 
			print(os.path.join(directory, filename))
			with open(directory+"/"+filename) as open_file:
				for line in open_file:
					if '>' in line:
						genecount+=1
						new_file.write(line)
					else:
						new_file.write(line)
	print(genecount)
