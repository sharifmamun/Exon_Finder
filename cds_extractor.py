# Got this from Mike on June 7, 2012  
# from
# http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/
# Formatting from:
# http://mail.python.org/pipermail/tutor/2003-July/024352.html

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq

Entrez.email = "mamun@cs.umanitoba.ca" ## replace
def index_genbank_features(gb_record, feature_type) :
    index = -1
    string = feature_type
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==string:
                return index
    return index

title = 'AJ973401'
handle2 = Entrez.efetch ( db='nucleotide',id= title,rettype='fasta')
seq_fasta = SeqIO.read(handle2,format='fasta')
handle1 = Entrez.efetch ( db='nucleotide',id= title,rettype='gb',retmode="text")
seq_gb = SeqIO.read(handle1,format='gb')
x = index_genbank_features(seq_gb,'CDS')
location = seq_gb.features[x].location
print 'CDS .. '
f = open('/home/mamun/Desktop/workfile.gbk', 'w')
f.write(str(seq_fasta[location.start:location.end]))
print repr(seq_fasta[location.start:location.end])
print seq_fasta[location.start.position:location.end.position]
#record_dict = SeqIO.index("workfile.gbk", "fasta")


#y = index_genbank_features(seq_gb,'source')
#print 'chromosome #',
#print seq_gb.features[y].qualifiers['chromosome']
handle1.close()
handle2.close()

print "\n"
print "CDS in Full Length:\n"
result = ''
f=open("temp2.fasta",'w')
for i in range(int(location.start.position),int(location.end.position)):
	result += str(seq_fasta[i])
#	print '\x08%s' %seq_fasta[i],
f.write(result)
print result
print "\n"
print "*************************************************************"
print "\n"
#print record_dict[location.i:location.(i+1)]

handle1.close()
handle2.close()
