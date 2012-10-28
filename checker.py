# developed on Bobolink
# on July 21, 2012
# checks the cds, chromosome, and est contigs
# before and during muscle alignment

#--------------------Import Stuffs-------------------#

from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO

import os
from Bio import SeqIO
from Bio.Seq import Seq
import curses
from curses import ascii

#-------------------End of all Imports---------------#

# The pupose of this function is
# To append EST Contig to the temporary opuntia.fasta file

def check(est_contig):
	handle1 = open("/home/mamun/Desktop/opuntia.fasta", "a")
	handle2 = open("/home/mamun/Dropbox/Campoletis_sonorensis.fasta", "r")
	temp_id = ""
	temp_seq = ""			
	parsedfasta = SeqIO.parse(handle2, "fasta")
	for i in parsedfasta:
	#       print i.id
               	temp_id=i.id
	        if (est_contig in str(i.id)):
        #       print i.seq
                	temp_seq = i.seq
	                break

	print "inside check FUNCTION\n"
	print temp_seq
	handle1.write("\n")
	handle1.write(">"+temp_id)
	handle1.write("\n")
	handle1.write(str(temp_seq))

	handle1.close()
	handle2.close()
	return str(temp_seq)

#-----------------------------------End of check function------------------------#
			


# Extracts the length of begin and end gaps
# Then removes them from the original Chroomosome/EST Contig
# 	if signal==0,  then works with Chromosome
#	if signal==1, then works with EST Contig
# & makes everything set for muscle alignment

def x_and_y(signal):
	handle1 = open("opuntia.aln", "r")
	global message2
	message2=""

#	str_val = str(value)
	parsedfasta = SeqIO.parse(handle1, "fasta")
	for i in parsedfasta:
	#       print i.id
	#       print i.seq
	        if ("CDS" in str(i.id)):
		        message=str(i.seq)
		        break

	temp = 0 
        # Taking the reverse complement
        # So that we can easily pick the gaps from end
        print "MESSAGE"
        message1=message[::-1]
     #  print message
        print "MESSAGE1"
     #  print message1

        # Counding the number of end gaps!
        print len(message)
        f = False
        for j in message:
           if f == False:
               if curses.ascii.isalnum(j):
                   f = True
                   temp+=1
               #    temp = temp + j
           else:
        #      temp = temp + j
               temp+=1
        #message = temp

        temp1 = 0
        temps=""
        f1 = False
        for j in message1:
           if f1 == False:
               if curses.ascii.isalnum(j):
                   f1 = True
                   temp1+=1
                   temps = temps + j
           else:
               temps = temps + j
               temp1+=1

        print len(temps)
        print "Beginning"
        print len(message1)
        start = (len(message) - temp)
        print start
        print "Ending"
        end = (len(message1) - temp1)
        print end
        handle1.close()

        message2=""
        message3=""
	fin_val=""
        print "handling fasta file"
        handle2 = open("opuntia.fasta", "r")
        parsedfasta = SeqIO.parse(handle2, "fasta")
        for i in parsedfasta:
                print i.id
		if ((signal==0) and ("Chromosome" in str(i.id))):
                        message2=str(i.seq)
                        print "Is it working?"
                #       print message2
                        fin_val=message2[start:-end]
                elif ((signal==1) and ("Campoletis" in str(i.id))):
                        message2=str(i.seq)
                        print "Is it working?"
                #       print message2
                        fin_val=message2[start:-end]
#               elif ((signal==2) and (str_val in str(i.id))):
#                        message2=str(i.seq)
 #                       print "Is it working?"
  #              #       print message2
   #                     fin_val=message2[start:-end]
                else:
                        print "LOOP ENDS\n"
# For testing purpose

        print "+++++++==============+++++++++++++++++++++++"
        print fin_val
#calling RECURSIVELY
        print "********************************\n"
        print "Recursive Function Working\n"
        print "********************************\n"
#       muscle_handler(cds, chromosome, message2[start:-end],0,0)
# Test Ends

        handle2.close()
#	return str(message2[start:-end])
	return fin_val

#---------------------------End of x_and_y---------------------------------#



def muscle_aligner():
        # from http://etal.myweb.uga.edu/biopywork.pdf
        # outfile=open("/home/mamun/Desktop/opuntia.aln",'w')
        muscle_cline=MuscleCommandline(input="/home/mamun/opuntia1.fasta", out="/home/mamun/opuntia1.aln")        
        print muscle_cline
        align=os.system(str(muscle_cline))

# For testing purpose
        #align = AlignIO.read(child.stdout, "fasta")
        #SeqIO.write()
#       print align
        #handling muscle alignment
#       print "ALIGNMENT\n"
# Test Ends
