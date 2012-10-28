#!/usr/bin/env python
# Developed By: S.M. Al Mamun
# Time period: May 1 - May 7
# This program Blasts over Internet using NCBIWWW, NCBIXML [from biopython]
# Input: EST sequences
# Output: 
#	Step 1: relevant sequences from 'nr' database
#	Step 2: extract the top ID from the extracted sequences	

from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from StringIO import StringIO
from Bio import Entrez
from Bio.Seq import Seq 
from math import fabs
import re
import os
import sys 
# copied from http://stackoverflow.com/questions/5000138/need-help-with-python-exception-handling
import time
import urllib
import  urllib2
#from urllib2.request import Request, urlopen
#from urllib.error import URLError, HTTPError
import httplib
RETRY_DELAY = 5

# working as grep
# finding a 7 or more digit value following by a word
# is being called for
# Entrez.efetch
def get_id_from_title ( title ):
  first_try = re.search('\w+\d{7,}',title )
  if first_try:
#   print str(first_try.group(0))
    return str(first_try.group(0))
  else: 
    return None

# Entrez email ID
# Must for accessing Entrez Database, change if you need!
Entrez.email = 'mamun@cs.umanitoba.ca'


#handling the input
if (len(sys.argv) != 3): 
  print 'usage: python get_EST.py <input_file.fas> yes/no'
  sys.exit()

name = sys.argv[1]
user_input = sys.argv[2]
if (user_input == "Yes" or user_input == "yes"):
	predict=True
else:
	predict=False

print 'processing file ' + name

# Extracting each ESTs from the main contig file 
# and processing them 
# maintaining a criteria of having Apis mellifera in them
# modification during: May 8 - May 13

fasta_sequences = SeqIO.parse(open(name),'fasta')
f = open('/home/mamun/Dropbox/output.txt', 'w')
result = open('/home/mamun/Dropbox/result_Campo1.txt', 'w')

#with open(result_file, "w") as f:
for seq in fasta_sequences:
        #if seq.id in wanted:
        #SeqIO.write([seq], f, "fasta")
        print seq.id
        print seq.seq
        print len(seq)

# Modified on May 8, 2012 on Bobolink
# processing EST Blast 
#	fasta_string = open(name).read()
	# the following for loop is for checking 10 times if 
	# for any reason the connection is lost or reset
	# due to network activity
	# http://stackoverflow.com/questions/5000138/need-help-with-python-exception-handling
	fasta_string = str(seq.seq)
	for x in range(10):
		try:
			result_handle = NCBIWWW.qblast("tblastx", "nt", fasta_string,hitlist_size=500,format_type='XML')			 
			print "working with NCBIWWW\n"
		except Exception, e:
			print "Inside Exception\n"
			time.sleep(RETRY_DELAY)
			print "errors : ", e
		finally:
		        print 'In finally block for cleanup'
			break
#		except urllib2.URLError, e:
#			print 'In exception 104\n'
#			if hasattr(e, 'reason'): # <--
#			        print 'We failed to reach a server.'
#			        print 'Reason: ', e.reason
#				time.sleep(RETRY_DELAY)
#				print 'getting delayed\n'
#			elif hasattr(e, 'code'): # <--
#			        print 'The server couldn\'t fulfill the request.'
#				print e.code
			#	if e.reason[0] == 104:
#				if e.code == 404: # Will throw TypeError if error is local, but we probably don't care
#				      print 'getting delayed\n'
#				      time.sleep(RETRY_DELAY)						   
#				else:
#
#				      raise # re-raise any other error
##		except urllib2.HTTPError as e:
##			print('The server couldn\'t fulfill the request.')
##		        print('Error code: ', e.code)
##			time.sleep(RETRY_DELAY)
##			print 'getting delayed\n'
##		except urllib2.URLError as e:
##		        print('We failed to reach a server.')
##		        print('Reason: ', e.reason)
##			time.sleep(RETRY_DELAY)
##			print 'getting delayed\n'
##		except httplib.BadStatusLine,msg:
##			print "%r" %(msg)
##			print 'getting delayed\n'
##			time.sleep(RETRY_DELAY)
#		except IOError:
#			print "Connection problem!"			
##		else:
##			# Checking - whether this section works properly or not!
##			print "breaking down"
##	 	        break # We've got resp sucessfully, stop iteration
	# Checking - whether this section works properly or not!
	print "Already broken :)"
#	handle = result_handle
#processing prasing stuffs from the Blast result
        if (result_handle!=""):
		for x in range(5):
			try:
                		blast_record = NCBIXML.read(result_handle)
			except Exception, e:
                	        print "Inside Result Handler Exception\n"
                        	time.sleep(RETRY_DELAY)		
			finally:
                        	print 'In finally block for cleanup'
	                        break
#all_records = NCBIXML.parse(handle)
#first_record = all_records.next()
#print first_record

# Commented on May 8, 2012
# f = open('/home/mamun/Dropbox/output.txt', 'w')
# processing the values based on E-Values
	#CHANGE threshold e-value if you need
	epsilon = 1e-5
	var = ""
	value2 = ""
	ss = ""

	#accessing all the records from BLAST output
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			# if (1) and (2) and (3)
			# 1 - checking threshold e-value
			# 2 - Contains Apis mellifera or not
			# 3 - Should not contain "PREDICTED:" in title
			# last condition of the following if condition was added later!
			if( ((hsp.expect < epsilon) and ("Apis mellifera" in str(alignment.title)) and ("PREDICTED:" not in str(alignment.title))) if predict else ((hsp.expect < epsilon) and  ("Apis mellifera" in str(alignment.title)))):
				value1 = '****Alignment**** | ', seq.id
                                s = str(value1)
                                value2 = 'sequence:', alignment.title
                                ss = str(value2)
                                value3 = 'length:', alignment.length
                                s = str(value3)
                                value4 = 'e value:', hsp.expect
                                s = str(value4)
                                value5 = hsp.query[0:75] + '...'
                                s = str(value5)
                                value6 = hsp.match[0:75] + '...'
                                s = str(value6)
#                               f.write('\n')
                                value7 = hsp.sbjct[0:75] + '...'
                                s = str(value7)
                                # breaking after the first iteration!
                                if "Apis mellifera" in ss:
                                        f.write(str(value1))
                                        f.write('\n')
                                        f.write(str(value2))
                                        f.write('\n')
                                        f.write(str(value3))
                                        f.write('\n')
                                        f.write(str(value4))
                                        f.write('\n')
                                        f.write(str(value5))
                                        f.write('\n')
                                        f.write(str(value6))
					f.write('\n')
					f.write(str(value7))
                                        f.write('\n')
					break
		if "Apis mellifera" in ss:
			break
#	else:
#		continue
#searching for BLAST HIT
#best_hit =  find_best_blast_hit (blast_record)
#print value2

	best_matched_seq = ""
	handle = StringIO("")
	
	# value2 is coming from the previous for loop!
	if (value2):
		try:
			handle = Entrez.efetch ( db='nucleotide',id= get_id_from_title(str(value2)) ,rettype='fasta')
			var = str(handle.read())
#			print handle
#		print id
#		if(handle):
#			best_matched_seq = SeqIO.read(handle,format='fasta')
#		best_matched_seq = SeqIO.read(handle, format='fasta')
#			best_matched_seq = SeqIO.parse(handle,format='fasta').next()
#		print best_matched_seq
			handle.close()
#	if (best_matched_seq):
		        print 'Here we go: \n'	
		#Edited on May 11, 2012	
#		result.write(str(best_matched_seq))
			result.write('\n')
			result.write(seq.id)
			result.write("*********************************")
			result.write('\n')
			result.write(var)
#		http://www.biostars.org/post/show/17652/efetch-and-biopython/
		except urllib2.HTTPError:
			print "Bad id"
		except IOError:
			print "Problem connecting to NCBI"
#	else:
#		print "Could NOT find expected result!"
#	print "*************"
#	print "%s with %i length" % (best_matched_seq.id, len(best_matched_seq))
