def csv_read(file1,file2,file3,file4,file5,file6,file7,file8,outfile):
	with open(outfile,'w') as outfile:
		linecount=0
		totalcount=0
		allcount=0
		with open(file1) as file1:
			for line in file1:
				newlist=line.split('\t')
				if linecount > 0 and newlist[7] != '#N/A':
					outfile.write(line)
					totalcount+=1
				linecount+=1
				allcount+=1
		outfile.write('\n')
		linecount=0
		with open(file2) as file2:
			for line in file2:
				newlist=line.split('\t')
				if linecount > 0 and newlist[7] != '#N/A':
					outfile.write(line)
					totalcount+=1
				linecount+=1
				allcount+=1
		outfile.write('\n')
		linecount=0
		with open(file3) as file3:
			for line in file3:
				newlist=line.split('\t')
				if linecount > 0 and newlist[7] != '#N/A':
					outfile.write(line)
					totalcount+=1
				linecount+=1
				allcount+=1
		outfile.write('\n')
		linecount=0
		with open(file4) as file4:
			for line in file4:
				newlist=line.split('\t')
				if linecount > 0 and newlist[7] != '#N/A':
					outfile.write(line)
					totalcount+=1
				linecount+=1
				allcount+=1
		outfile.write('\n')
		linecount=0
		with open(file5) as file5:
			for line in file5:
				newlist=line.split('\t')
				if linecount > 0 and newlist[7] != '#N/A':
					outfile.write(line)
					totalcount+=1
				linecount+=1
				allcount+=1
		outfile.write('\n')
		linecount=0
		with open(file6) as file6:
			for line in file6:
				newlist=line.split('\t')
				if linecount > 0 and newlist[7] != '#N/A':
					outfile.write(line)
					totalcount+=1
				linecount+=1
				allcount+=1
		outfile.write('\n')
		linecount=0
		with open(file7) as file7:
			for line in file7:
				newlist=line.split('\t')
				if linecount > 0 and newlist[7] != '#N/A':
					outfile.write(line)
					totalcount+=1
				linecount+=1
				allcount+=1
		outfile.write('\n')
		linecount=0
		with open(file8) as file8:
			for line in file8:
				newlist=line.split('\t')
				if linecount > 0 and newlist[7] != '#N/A':
					outfile.write(line)
					totalcount+=1
				linecount+=1
				allcount+=1
		print("total features appearing in Trinity: "+str(totalcount))
		print("proportion features used : "+str(totalcount / allcount))
				

def unique_genes(inputfile):
	with open(inputfile) as inputfile:
		feature_dictionary={}
		newlist=[]
		for line in inputfile:
			newlist=line.split('\t')
			try:
				key=newlist[6]
				feature_dictionary[key]=1
			except:
				continue
		print(len(feature_dictionary))
		return feature_dictionary
		
def compare_libs(lib1,lib2):
	common=0
	for i in lib1:
		for j in lib2:
			if i == j:
				common+=1
	print("number of common genes overall is: "+str(common))		
			
if __name__ == '__main__' :
# 	csv_read('ange_MSMI_13464770_matchtrimmed.txt','ange_MSMI_13464771_matchtrimmed.txt','ange_MSMI_13464772_matchtrimmed.txt','ange_MSMI_13464773_matchtrimmed.txt','ange_MSMI_13464774_matchtrimmed.txt','ange_MSMI_13469182_matchtrimmed.txt','ange_MSMI_13469183_matchtrimmed.txt','ange_MSMI_13469185_matchtrimmed.txt','ange_combined.txt')
# 	csv_read('sonu_MSMI_13464433_matchtrimmed.txt','sonu_MSMI_13464434_matchtrimmed.txt','sonu_MSMI_13469184_matchtrimmed.txt','sonu_MSMI_13469186_matchtrimmed.txt','sonu_MSMI_13512650_matchtrimmed.txt','sonu_MSMI_13512651_matchtrimmed.txt','sonu_MSMI_13512652_matchtrimmed.txt','sonu_MSMI_13512653_matchtrimmed.txt','sonu_combined.txt')
	libandro=unique_genes('ange_combined.txt')
	libsonu=unique_genes('sonu_combined.txt')
	compare_libs(libandro,libsonu)
		