# from Bio import Entrez
# 
# def file_read(filename):
# 	with open(filename) as file_handle:
# 		newlist=[]
# 		for line in file_handle:
# 			newlist.append(line.strip())
# 		return newlist
# 
# if __name__ == '__main__' :
# # 	idfile=file_read('/Users/avahoffman/Documents/CSU/Research/msmi_arrays/CEE_microarrays/Batch_entrez_trial.txt')
# # 	print(idfile)
# 	Entrez.email = "avamariehoffman@gmail.com"     # Always tell NCBI who you are
# 	handle = Entrez.efetch(db="EST", id="CD650876", rettype="fasta", retmode="text")
# 	print(handle.read())
	
import urllib2, urllib
url="https://www.ncbi.nlm.nih.gov/nucest/"
params={
    'name':'CD650876'
    }
data=urllib.urlencode(params)
headers = {"Accept" : "*/*"}
req = urllib2.Request(url, data, headers)
print(urllib2.urlopen(req).read())