from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import os.path
import time
from os.path import commonprefix
import sys
import shutil
import os
import glob

import checker

#subgroups=[] 
#all_results=[] 
#results_to_cds=[]
#results_to_DNA=[]

cds = ""
dna = ""

#cds = Seq("ATGTATGCTTGTACTAGATTTATTGTACCTGTCGCCAAATCGACCTTGGTCTCTGGAACTAAGAATTATATTCGACCATTGAGCAGTGCTGTGGTAAGTCACAGCCAGTCTATACAGCAGAATCAGATTCAGAACCCATTCAGTTCTCCATTGGGACAATCTTCAATTATACGTAATTTCCAAACATCTACAATTAGTCGTGATATTGACTCTGCTGCAAAATTCATTGGTGCTGGTGCTGCAACTGTAGGTGTTGCTGGATCGGGAGCTGGGATTGGTTCAGTGTTTGGTTCTTTAATCGTAGGCTATGCCAGGAATCCTTCTCTTAAGCAACAGTTGTTTTCATATGCCATTTTGGGCTTTGCCTTGTCTGAAGCCATGGGTCTTTTTTGCCTTATGATGGCTTTCTTACTTTTGTTCGCGTTCTAA", IUPAC.ambiguous_dna)

#dna=Seq("ATGTATGCTTGTACTAGATTTATTGTACCTGTCGCCAAATCGACCGTAAGTAATTTTCTCTTTATTTCGAATATTTTCAAATTTTTATATTATTATATAATCTTTGAGTAATTTGTTTTATTTTCAGTTGGTCTCTGGAACTAAGAATTATATTCGACCATTGAGCAGTGCTGTGGTAAGTCACAGCCAGTCTATACAGCAGAATCAGATTCAGAACCCATTCAGTGTAAGTAACATTTAATATTTTTCGCATTGCATAACCTCACCCCATTACGAGGTTACATAAGCAAGTTTATTTTCATAATTCAAAGATTTTTTTTAGTAAGAAATAAAAAATTTTTTTAAATAAATCAATTTGTAATAAACATGTATGATCTTATATTTCCAGACTTTTCATCTTTTGAATTGATTTTAAATTAACTAATCATTTTATAGTTTAGGATTATGTCATGAAAAATATTAATTTAAATTGAAAAAAGATTTGATTTAAATGGAATAGTTGAAATAAATAATTTAAAATATAAGATGAATATAAAATTCTAAAATAATTGTCTGTACAAAATTTATTTTTTATGTCAAAATATTAAAACAAAAATAATAAAGATTATATAAGTTGTATTATTTGATGAAATAATTTACTGTATTTAATATATAACTGAAATCATTATTAAATTTGATAATTTTTTATTAAATTTTTTTTTAAATTTAAATATTTTTTATATTTCGTAAGATTTTTAAATATGAAAAATTTTTTCTTTAATATTACTAATTTTTTTCTAATTATTTGAATAATTTTTATATAGTCTCCATTGGGACAATCTTCAATTATACGTAATTTCCAAACATCTACAATTAGTCGTGATATTGACTCTGCTGCAAAATTCATTGGTGCTGGTGCTGCAACTGTAGGTGTTGCTGGATCGGGTATGTATATGATTTTTATTTTTGATAACTATAAAAATTATAATATTTTGATTCCTATTCCTATATTCATTTTTCTTGTTAAAAAAAATATAAAATTTACATTCCAGGGATTGATATGATAAATATTTAATTTTTGGATCCCTAATTGGATATACCCAAAATCCATCCTTAAAGCAACAAATTTCCAATTATGCTATCCTTGGCTTTGTATTATCAGAAGCAATGGACTTTTTTTATTTGATGATAGCATTCTTTATATTATTTGCCTTCTGATCCTAATTCTAAAAAGAATCATCCTAATACCATCATAGTAGTAAAGTTATATGATAAATTATATAAATTTAAATAATTTTTTTTTTAATTATAATAAATAAGATTATACATGTTATGATTTCAGGAGCTGGGATTGGTTCAGTGTTTGGTTCTTTAATCGTAGGCTATGCCAGGAATCCTTCTCTTAAGCAACAGTTGTTTTCATATGCCATTTTGGGCTTTGCCTTGTCTGAAGCCATGGGTCTTTTTTGCCTTATGATGGCTTTCTTACTTTTGTTCGCGTTCTA", IUPAC.ambiguous_dna)

def LongestCommonSubstring(S1, S2): 
    M = [[0]*(1+len(S2)) for i in xrange(1+len(S1))] 
    longest, x_longest = 0, 0 
    for x in xrange(1,1+len(S1)): 
        for y in xrange(1,1+len(S2)): 
            if S1[x-1] == S2[y-1]: 
                M[x][y] = M[x-1][y-1]+1 
                if M[x][y] > longest: 
                    longest = M[x][y] 
                    x_longest = x 
                if longest >= 3: 
                    subgroup=S1[x_longest-longest:x_longest] 
                    subgroups=set([subgroup]) 
                    print(subgroups) 
            else: 
                    M[x][y] = 0 

    return S1[x_longest-longest:x_longest] 


# Defining two global variables to work with both CDS and DNA
tracker_cds = 0
tracker_dna = 0
chrome = ""
CDS = ""
CHROME = ""


def work(val1, val2):
#Working with CDS 

#Working with Apis Chromosome

	print cds
##	temp_rev = dna[::-1].translate()
	temp_cds = val1.translate()
	temp_dna = val2.translate()

	global tracker_cds
	global tracker_dna
#	global dna
	global chrome
	global CDS
	global CHROME
	print temp_cds
	#time.sleep(10)
	i=0
	result=""
	temp=""
	temp_X=0
	max_match=-1
	result=""
	print "tracker_cds:", tracker_cds
	print "tracker_dna:", tracker_dna
	#time.sleep(2)

	ret_val = ""
	
	while (i<3):
		temp_dna = (val2[i:]).translate()
		print temp_dna
		
		# Getting the hypothesis value
		# For example, the minimum length of Exon is 6
		print temp_cds[1:5]
	
		max_match = temp_dna.find(temp_cds[1:5])
		
		print max_match
		if (max_match>=0):
			print "Showing DNA>>"
			print temp_dna[max_match:], "No of shift", i
                        if (i==1):
                                #CDS += "--"
                                CHROME += "~~"
                        elif (i==2):
                                #CDS += "-"
                                CHROME += "~"

			try:
				result=commonprefix([temp_cds[1:], temp_dna[max_match:]])
#				all_results.append(str(result))
#				all_results.append(i)
#				results_to_DNA.append(str(val2[i:(max_match-1+len(result))*3+i]))
				if max_match==0:
					chrome += "---"

				chrome += str(val2[i:(max_match+len(result))*3+i])
				CHROME += str(val2[0:(max_match+len(result))*3+i])

				#chrome += str("*")
				print "CHROME IS WORKING..........."				
				#if len(chrome)%3==1:
				#	chrome += "~~"
				#elif len(chrome)%3 == 2:
				#	chrome += "~"
				print chrome
				#time.sleep(3)
				#results_to_DNA.append(str(dna[i+tracker_dna:(max_match+len(result))*3]))
#				results_to_DNA.append(len(val2[i:(max_match-1+len(result))*3+i]))
				# Appending to subgroups
				# Here, subgroups are the translated Chromosome strands
				#results_to_cds.append(str(val1[0:len(result)*3]))
#				subgroups.append(str(temp_dna[0:max_match-1+len(result)]))
			except:
				print "IN EXCEPTION"
#				all_results.append(str(temp_dna[max_match:]))
#				all_results.append(i)
#				print all_results
				#subgroups.append(str(dna[0:]))
				print "-------------------"
#				print subgroups
				print "\nConverting into DNA>>>\n"
#				results_to_DNA.append(str(val2[i:]))

				chrome += str(val2[i:])
				print str(val2[i:])
				CHROME += str(val2[0:])
#				print results_to_DNA
				#print results_to_cds
				break
	
			finally:			
				print "---------------------"
				print "Common Prefix:", result
				print "*********************"
				print temp_cds[0]+result
				print (len(result)+1)*3
				tracker_cds = tracker_cds + ((len(result) + 1)*3)
				tracker_dna = tracker_dna + ((max_match+len(result))*3) + i
				#if (i==1):
				#	CDS += "~~~"
				CDS += str(((max_match-1)*3)*"~")
				if (i>0 and max_match > 0):
					CDS += "~~~"
				CDS += str(val1[0:(len(result)+1)*3])
					

			print tracker_cds
			print tracker_dna
			print "Test........"
			print chrome
			#time.sleep(2)
			if (tracker_cds<len(cds)-9):
				work(cds[tracker_cds:], dna[tracker_dna:])
			else:
				CDS += str(val1[tracker_cds:])
			#else:
			#	if tracker_cds<len(cds):
			#		cds += 
#		print "Outside of While Loop"
#		print result
#		print i	
		i+=1
	
	print "We are done!"
	#sys.exit()
	print tracker_cds
	print len(cds)

def main():
	#print cds.translate(cds=True)
	#chrome=work(cds, dna)
	global cds
	global dna	
	global chrome
	global CDS
	global CHROME
		
	global tracker_cds
	global tracker_dna

	path = "/home/mamun/Desktop/Output_test/"
	for infile in glob.glob( os.path.join(path, '*.fasta') ):
	
		handle1 = open(infile, "r")
		parsedfasta = SeqIO.parse(handle1, "fasta")

		for i in parsedfasta:
		        if (i.id=="CDS"):
                		cds = Seq(str(i.seq),IUPAC.ambiguous_dna)
				print "CDS is", cds
		        elif (i.id=="Chromosome"):
                		dna = Seq(str(i.seq), IUPAC.ambiguous_dna)
				print "DNA is", dna
		        else:
                		campoletis = Seq(str(i.seq), IUPAC.ambiguous_dna)
				print campoletis

		## making empty list!!
		## for the global variables!!

		tracker_cds = 0
		tracker_dna = 0
#		chrome = ""
		pos = 0

		work(cds, dna)
		#if(chrome[-7:]!=(cds[-7:]) and (len(chrome)<len(cds))):
		#	try:
		#		pos = chrome[-7:].find(cds)
		#		chrome += cds[pos:]
		#	except:
		#		pass

		print "New DNA", dna
		print "New Chromosome", chrome
		
		print "********************"
		print "current file is: " + infile
        	print "** Working with ALIGNMENT file **"
	        handle2 = open("opuntia1.fasta", "w")
        	handle2.write(">CDS")
	        handle2.write("\n")
        	handle2.write(str(cds))
	        handle2.write("\n")
        	handle2.write(">Chromosome")
	        handle2.write("\n")
        	handle2.write(str(chrome))
	        #handle2.write("\n")
        	#handle2.write(">EST")
	        #handle2.write("\n")
        	#handle2.write(str(campoletis))
	        handle2.close()

		checker.muscle_aligner()

		## Checking whether there is a stop codon or not
		## at the very beginning of the EST Contig!
		## My hypothesis is to check within the first 20 seq
		temp = Seq(str(campoletis), IUPAC.ambiguous_dna)
		i=0
		while(i<3):
			temp1 = temp[i:].translate()
			value = temp1[0:20].find("*")
			if (value==-1):
				break
			else:
				i+=1
		print CDS
		print chrome
                handle4 = open("manual.aln", "w")
                handle4.write(">CDS")
                handle4.write("\n")
                handle4.write(str(CDS))
                handle4.write("\n")
                handle4.write(">Chromosome")
                handle4.write("\n")
                handle4.write(str(CHROME))
       		handle4.close()

		handle3 = open("opuntia2.aln", "w")
		handle3.write(">EST")
                handle3.write("\n")
                handle3.write(str(campoletis[i:]))
		handle3.close()
		print "About to end!!"
		## Editing for profile-profile alignment
		os.system("muscle -profile -in1 manual.aln -in2 opuntia2.aln -out combinedAlignment.fasta -maxmb 15000")
		#muscle_cline=MuscleCommandline(input="opuntia1.fasta", out="opuntia1.afa")
        	#print muscle_cline
	        #align=os.system(str(muscle_cline))
		
        	#checker.muscle_aligner()
		shutil.copy2('combinedAlignment.fasta', infile+'.aln')

		#time.sleep(4)
		chrome = ""
	print "------------End of main---------------"
	#print cds.translate()
	
main()

