# Last modified on July 17, 2012
# on Bobolink
# handling the Tblastx output upto muscle alignment

###################################################################
#---------------------All imported stuffs here---------------------
import re
import os
import string
import time
import shutil
import shlex, subprocess
import curses
import sys, subprocess
from StringIO import StringIO
from Bio import AlignIO
from curses import ascii
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Emboss.Applications import NeedleCommandline

# importing sample program
# make sure, following two scripts are in the same directory
import sample
import checker

# for muscle alignment
from Bio.Align.Applications import MuscleCommandline

#------------------------ End of Imports---------------------------
####################################################################


# grep works similar as the UNIX/LINUX grep command
# receives two argument
# returns the position

def grep(pattern,fileObj):
        r=[]
# If we need the line number, then we can use the following 6 lines!
#       linenumber=0
#       for line in fileObj:
#               linenumber +=1
#               if re.search(pattern,line):
#                    r.append((linenumber,line))
#       return r
        for line in fileObj:
                if re.search(pattern,line):
                        r.append(line)
        return r


# While working with Entrez Database
# To connect, we always need email address
# Please, modify with your own
Entrez.email = "mamun@cs.umanitoba.ca" ## replace

def extractor(infile):
        """
        1. This Function is checking all the directories in Chromosomes Directory
        2. Then checking the Span/Range of intron/exon in that .est2genome file
        3. Extracting the Start and End point of the region
        ***************************************************************************
        4. Then based on the file name, this function will search a chromosome file
        5. Finally, will extract the region from the .est2genome file
        """

        # **************************************
        # This part handles the directory issuue
        # **************************************
        path='/home/mamun/Desktop/'+infile
        #listing = os.listdir(path)
        print "current file is: " + infile
        chrome_temp=infile.split('.')
        chrome=chrome_temp[0].split('_')
        print chrome[3]
        f = open(path, 'r')
        read_data = f.readlines()
        temp=grep("Span", read_data)

        #***********************************
        #this part extracts 4th & 5th column
        #***********************************
        allrows=[]
        #       for item in temp:
        #               allrows.append = temp[item]
        line=[i.split()for i in temp]
        string = str(line)
        line=string[1:-1].split(' ')
        #       print 
        print temp
        begin=int(line[3][1:-2])
        end=int(line[4][1:-2])
        #       print line[3][1:-2]
        #       print line[4][1:-2]
        print begin
        print end
        print "\n"
        print "****************************"
        print "****************************"
        # End of directory traversing

        path_temp='/home/mamun/Desktop/Thesis_Work/Chromosomes/'+str(chrome[3])+'.fa'
        genome=SeqIO.read(open(path_temp),'fasta') #you MUST tell SeqIO what format is being read
        #result=extractseq(str(genome), 3546620, 3547507)
        #following lines are for testing purpose
        #print result
        #if(infile=="nm_lg16_001163463.est2genome"):
        print genome.seq[begin-1:end-1]
	return genome.seq[begin-1:end-1]

#----------------End of extractor-------------------#


# This function receives
# GenBank sequence and feature
# & returns the index of that feature

def index_genbank_features(gb_record, feature_type) :
    index = -1
    string = feature_type
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==string:
                return index
    return index

#----------------End of index_genbank_features--------------------#



# Funcion created on June 14, 2012
# This function 
# receives different parts to create est2genome file name
# using concatenation

def est_genome(cam_no,lg_no, genome1, est):
	out= cam_no+"_"+lg_no+".est2genome"
	print "testing........."
	f=open("temp.fasta",'w')
	f.write(est)
#	f1=open("temp.fasta",'r')
#	print str(genome1)
	f.close()
	#args= ["est2genome", "temp.fasta", genome, "-outfile", out]

# if the length is too long, then ignore it. No need to run est2genome anymore!
	if(len(est)<3000):
		output=os.system("est2genome temp.fasta %s %s" %(genome1,out)) # Success!
	else:
		output = "GARBAGE"
	print output
#	f1=open(out,'w')
#	f1.write(p)
	return output
#---------------End of est_genome--------------------------------------------------------#


# This function works with muscle handler
# works with cds, chromosome, and finally with EST Contigs
# Right now, begin and end are useless! We will remove them later on

def muscle_handler(cds, chromosome, est_contig,begin,end):
	handle2 = open("/home/mamun/Desktop/opuntia.fasta", "w")
	global value
	global value1
	global chrome
	chrome=str(chromosome)
	
	print "+++++++++++++++++++++++++++++++++++++++\n"
	print chrome
	print "+++++++++++++++++++++++++++++++++++++++\n"
# Step 1: 
# Open opuntia.fasta as temporary file
# Write CDS and Chromosome one after another to that file

	handle2.write(">CDS")
	handle2.write("\n")
	handle2.write(cds)
	handle2.write("\n")
	handle2.write(">Chromosome")
	handle2.write("\n")
	handle2.write(chrome)
	handle2.write("\n")
	handle2.close()
# Close the temporary file

#calling muscle_aligner() from checker.py
	checker.muscle_aligner()

# calling cds_handler 1st time
# it returns CDS values after alignment

	message=cds_handler()	
	value1=""
	value=""
	#TESTING validity
#	sample.func()

# if beginning or end of the CDS contains any gap
#	then call checker.x_and_y - here, 0 argument means I am handling CDS
#				  - for more details see checker.x_and_y
# 	remove the beginning and end gaps from the CDS sequence
# Then write the new CDS sequence, and Chromosome again to the temporary file
#				  - here, temporary file is opuntia.fasta
 
	if(message[0]=="-" or message[:-1]=="-"):
		value1=checker.x_and_y(0)
		handle2 = open("/home/mamun/Desktop/opuntia.fasta", "w")
	        #chromosome=str(value)
        	handle2.write(">CDS")
	        handle2.write("\n")
        	handle2.write(cds)
	        handle2.write("\n")
        	handle2.write(">Chromosome")
		handle2.write("\n")
		print "Inside CDS and Chromosome"
		print value1
	        handle2.write(value1)
        	handle2.close()

# else there will be no change in the chromosome

	else:
		print "CDS and Chromosome part is done!\n"
		value1=chrome

#calling check() function to insert est contig into fasta file
	value=checker.check(est_contig)

#aligning with all three (CDS, Chromosome, and EST Contig)
##	checker.muscle_aligner()
	print "\nHandler Number: 2\n"
	#time.sleep(4)
#	message1=cds_handler()

# argument 1 means, we are handling EST Contig in checker.x_and_y
#	test=checker.x_and_y(1)
	print "*********************************\n"
	print "TESTing.........................\n"
#	print test
        #TESTING validity
#       sample.func()
#        if(message1[0]=="-" or message1[:-1]=="-"):
        #value=checker.x_and_y(1)
	print value
	os.remove("/home/mamun/Desktop/opuntia.fasta")
#	else:
#		print "Life goes on!"
	#time.sleep(3)
	handle2 = open("/home/mamun/Desktop/opuntia.fasta", "w")
	handle2.write(">CDS")
        handle2.write("\n")
	handle2.write(cds)
	handle2.write("\n")
        handle2.write(">Chromosome")
	handle2.write("\n")
	if (value1==""):
		handle2.write(str(chrome))
	else:
		handle2.write(str(value1))
	handle2.write("\n")
	handle2.write(">"+est_contig)
	handle2.write("\n")
	print "Inside EST part"
	print "Chrome......................"
	print value1
	print "EST.........................."
	print value
	#Getting first occurane position of ATG
	if value.find('ATG'): 
		position=value.find('ATG')
		temp_value=value[position:]
	       	handle2.write(str(temp_value))
		value=str(temp_value)
	else:
		handle2.write(str(value))
	handle2.close()
	checker.muscle_aligner()
	print "COMPLETE ALIGNMENT is DONE"

	message1=cds_handler()
	if(message1[0]=="-" or message1[:-1]=="-"):
        	value=checker.x_and_y(1)
		print value
		print "Getting perfect EST\n"
		#time.sleep(5)
	 
		handle2 = open("/home/mamun/Desktop/opuntia.fasta", "w")
	        handle2.write(">CDS")
	        handle2.write("\n")
        	handle2.write(cds)
	        handle2.write("\n")
        	handle2.write(">Chromosome")
	        handle2.write("\n")
		if (value1==""):
                	handle2.write(str(chrome))
        	else:
                	handle2.write(str(value1))
		handle2.write("\n")
	        handle2.write(">"+est_contig)
	        handle2.write("\n")
		handle2.write(str(value))
		handle2.close()
		checker.muscle_aligner()
	else:
		print "We are really done!\n"
		#time.sleep(3)
  
	## Working with Amino Acid Alignment
 
# End of iteration and alignment

# Copying temp files to seperate Directory
	shutil.copy2('opuntia.fasta','/home/mamun/Desktop/Output/'+est_contig+'.fasta')
	shutil.copy2('opuntia.aln','/home/mamun/Desktop/Output/'+est_contig+'.aln')
	
	handle2 = open("/home/mamun/Desktop/opuntia.fasta", "w")
	coding_dna1 = Seq(cds, IUPAC.ambiguous_dna)	
	handle2.write(">CDS")
        handle2.write("\n")
	handle2.write(str(coding_dna1.translate()))
        handle2.write("\n")
	handle2.write("\n")
        handle2.write(">Chromosome")
        handle2.write("\n")
        if (value1==""):
		coding_dna2 = Seq(chrome, IUPAC.ambiguous_dna)
                handle2.write(str(coding_dna2.translate()))
		handle2.write("\n")
		handle2.write("\n")
        else:
		coding_dna2 = Seq(value1, IUPAC.ambiguous_dna)
                handle2.write(str(coding_dna2.translate()))
	        handle2.write("\n")
		handle2.write("\n")
        handle2.write(">"+est_contig)
        handle2.write("\n")
	coding_dna3 = Seq(value, IUPAC.ambiguous_dna)
	if (coding_dna3!=""):
        	handle2.write(str(coding_dna3.translate()))
	else:	
		handle2.write("")
        handle2.close()
        checker.muscle_aligner()
  
	# Copyting current file to system!  

	shutil.copy2('opuntia.fasta','/home/mamun/Desktop/Output/'+est_contig+'_acid.fasta')
        shutil.copy2('opuntia.aln','/home/mamun/Desktop/Output/'+est_contig+'_acid.aln')
#--------------------------------muscle_handler---------------------------------#


# This function opens the temporary alignment file
# returns the CDS part

def cds_handler():
	handle1 = open("/home/mamun/Desktop/opuntia.aln", "r")
        parsedfasta = SeqIO.parse(handle1, "fasta")
        for i in parsedfasta:
                print i.id
                print i.seq
                if ("CDS" in str(i.id)):
                        message=str(i.seq)
                        break
  	handle1.close()
	print "////////////////////////////"
	#print message
	return message

#----------------------End of cds_handler()-------------------#



def main():
	handle = open("/home/mamun/Dropbox/result_Campo1.txt", "r")
#	parsedfasta = SeqIO.parse(handle, "fasta")
	temp=[]
	#blast_records = NCBIXML.parse(result_handle)
	#blast_record = blast_records.next()
	#save_file = open("/home/mamun/Dropbox/result_handler.txt", 'w')
	read_data = handle.readlines()
	print "test1"
	temp=grep("Campoletis_sonorensis", read_data)
	j=0
	print temp
	print "test2"
	handle1 = open("/home/mamun/Dropbox/result_Campo1.txt", "r")
	parsedfasta = SeqIO.parse(handle1, "fasta")
	new_file = open('/home/mamun/Desktop/notfound.txt', 'a')
	chromosome_file1 = '/home/mamun/Desktop/Thesis_Work/Chromosomes/'

#	Handling result file
#	getting >gi|......|ref|ID.?|............
#	extracting ID from the fasta format file
	for i in parsedfasta:
            #save_file.write('>%s\n' % (alignment.title,))
		print i.id   
		print "=================================="
		temp_string = str(i.id)
		format_db=string.split(temp_string,'|')
		print format_db[3]
		temp_val=str(format_db[3])
		temp_db=string.split(temp_val,'.')
		print temp_db[0]
		
		handle2=""
		handle1=""	
		title = str(temp_db[0])
		for x in range(10):
			try:
				handle2 = Entrez.efetch ( db='nucleotide',id= title,rettype='fasta')
				seq_fasta = SeqIO.read(handle2,format='fasta')
				handle1 = Entrez.efetch ( db='nucleotide',id= title,rettype='gb',retmode="text")
				seq_gb = SeqIO.read(handle1,format='gb')
			except Exception, e:
				print "Inside Exception\n"
				print "errors : ", e
				time.sleep(1)
			finally:
				print 'In finally block for cleanup'
				break
		
		if (handle1=="") or (handle2==""):
			break
		x = index_genbank_features(seq_gb,'CDS')
		location = seq_gb.features[x].location
		print 'CDS .. '

		# For testing purpose
		f = open('/home/mamun/Desktop/workfile.gbk', 'w')
		f.write(str(seq_fasta[location.start:location.end]))
		# test Ends

		print repr(seq_fasta[location.start:location.end])
		print seq_fasta[location.start.position:location.end.position]
		#record_dict = SeqIO.index("workfile.gbk", "fasta")

		print "\n"

	        # Working with Contig Number
                print temp[j]
                temp_contig= str(temp[j])
                contig_db=string.split(temp_contig,'*')
                print contig_db[0]

		# Working with features
		# Finding out, which chromosome to work with

		y = index_genbank_features(seq_gb,'source')
		print 'chromosome #',
		print seq_gb.features[y]
		s=str(seq_gb.features[y])
		temp_str=""
		#print s.find('chromosome,')
		#print  seq_gb.features[y].qualifiers.get('chromosome')

		# if feature contains 'chromosome' or 'linkage_group'
		#	then work with that as a chromosome number
		# else
		#	not sure!

		if(seq_gb.features[y].qualifiers.get('chromosome')):
			if (str(seq_gb.features[y].qualifiers.get('chromosome')).lower()[2:-2]=='unknown'):
				print "Unknown: Not working"
			else:
				temp_str=str(seq_gb.features[y].qualifiers['chromosome']).lower()[2:-2]
				print temp_str
		elif(seq_gb.features[y].qualifiers.get('linkage_group')):
			temp_str=str(seq_gb.features[y].qualifiers['linkage_group']).lower()[2:-2]
			print temp_str
		else:
			new_file.write(str(contig_db[0]))
			new_file.write("\n")
			new_file.write(str(temp_db[0]))
			new_file.write("\n")
		print "CDS in Full Length:\n"

		# Working with Complete CDS
		result = ''
		
		# Concatenating all into one :)
		for i in range(int(location.start.position),int(location.end.position)):
		        result += str(seq_fasta[i])
		#       print '\x08%s' %seq_fasta[i],
		# DONE, result variable contains our CDS
		print result

		est2gen_output=""
		# Working est2genome
		if (temp_str!=""):
			chromosome_final= chromosome_file1+temp_str+".fa"
			est2gen_output = est_genome(str(contig_db[0]),temp_str,chromosome_final, result)
		else:
			print "*********************************************\n"
			print "Keys do not include chromosome/linkage group!\n"
			print "*********************************************\n"
			pass

		print est2gen_output
		
		# For testing purpose
		print "\n"
		print "*************************************************************"
		print "\n"
		# Test Ends

		#print record_dict[location.i:location.(i+1)]
#Calling extractor
		param=contig_db[0]+"_"+temp_str+".est2genome"

#Edited on July 19, 2012
		# To check both:
		# 1. File exists or not
		# 2. if exists, then is it empty or not!		
		est_result=""
		if os.path.exists("/home/mamun/Desktop/"+param) and os.stat("/home/mamun/Desktop/"+param)[6]!=0:
			est_result=extractor(param)
			if len(est_result)<5000:
				muscle_handler(str(result), str(est_result), str(contig_db[0]),101,101)
			else:
				print "Chromosome length is too LONGER!"
		else:
			print param + " does not exist in the directory"
			print "For MUSCLE alignment, we will need both CDS and Chromosome!"
#Calling muscle function for the alignment
#		muscle_handler(str(result), str(est_result), str(contig_db[0]))

		handle1.close()
		handle2.close()
		
		# For testing purpose
                j+=1
#		time.sleep(3)
#		print temp
		# Test ends

#	handle.close()
	print "test3"
	print j
#	new_file.close()

#------------------------------------------------------End of main---------------------------------------------------------#


main()
