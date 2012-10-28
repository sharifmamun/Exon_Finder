#!/usr/bin/env python
#	Program was made by Michael Domartzki 
#	Some time in 2011 for Barb Sharownoski
#	I received this code on April 26, 2012 from Mike
#	I will be working on this code in upcoming days
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from math import fabs
import re
import os
import sys 
try: 
  from cStringIO import StringIO 
except ImportError: 
  from StringIO import StringIO 

def _as_string(s): 
  if isinstance(s, str): 
    return s 
  return s.decode() 

#def _mike_qblast(program, blast_program, database, sequence, auto_format=None,composition_based_statistics=None, db_genetic_code=None,endpoints=None,entrez_query='(none)', expect=10.0,filter=None,gapcosts=None,genetic_code=None, hitlist_size=50,i_thresh=None,layout=None,lcase_mask=None, matrix_name=None,nucl_penalty=None,nucl_reward=None, other_advanced=None,perc_ident=None,phi_pattern=None, query_file=None,query_believe_defline=None,query_from=None, query_to=None,searchsp_eff=None,service=None,threshold=None, ungapped_alignment=None,word_size=None, alignments=500,alignment_view=None,descriptions=500, entrez_links_new_window=None,expect_low=None,expect_high=None, format_entrez_query=None,format_object=None,format_type='XML', ncbi_gi=None,results_file=None,show_overview=None, megablast=None, template_length=None,template_type=None): 
# Edited on April 29, 2012 [by Mamun]
def _mike_qblast(program, database, sequence, auto_format=None,composition_based_statistics=None, db_genetic_code=None,endpoints=None,entrez_query='(none)', expect=10.0,filter=None,gapcosts=None,genetic_code=None, hitlist_size=50,i_thresh=None,layout=None,lcase_mask=None, matrix_name=None,nucl_penalty=None,nucl_reward=None, other_advanced=None,perc_ident=None,phi_pattern=None, query_file=None,query_believe_defline=None,query_from=None, query_to=None,searchsp_eff=None,service=None,threshold=None, ungapped_alignment=None,word_size=None, alignments=500,alignment_view=None,descriptions=500, entrez_links_new_window=None,expect_low=None,expect_high=None, format_entrez_query=None,format_object=None,format_type='XML', ncbi_gi=None,results_file=None,show_overview=None, megablast=None, template_length=None,template_type=None): 
  import urllib, urllib2 
  import time 
# FOR FUTURE REFERENCE: this is the qblast from biopython, except that it has one added parameter:
# BLAST_PROGRAM.
# no other code was changed -- mike d., may 31/11
# oh wait, i added template Length and template type jun 1/11

# Format the "Put" command, which sends search requests to qblast. 
# Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node5.html on 9 July 2007 
# Additional parameters are taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node9.html on 8 Oct 2010 
  parameters = [ 
    ('AUTO_FORMAT',auto_format), 
    ('COMPOSITION_BASED_STATISTICS',composition_based_statistics), 
    ('DATABASE',database), 
    ('DB_GENETIC_CODE',db_genetic_code), 
    ('ENDPOINTS',endpoints), 
    ('ENTREZ_QUERY',entrez_query), 
    ('EXPECT',expect), 
    ('FILTER',filter), 
    ('GAPCOSTS',gapcosts), 
    ('GENETIC_CODE',genetic_code), 
    ('HITLIST_SIZE',hitlist_size), 
    ('I_THRESH',i_thresh), 
    ('LAYOUT',layout), 
    ('LCASE_MASK',lcase_mask), 
   # ('MEGABLAST',megablast), 
    ('MATRIX_NAME',matrix_name), 
    ('NUCL_PENALTY',nucl_penalty), 
    ('NUCL_REWARD',nucl_reward), 
    ('OTHER_ADVANCED',other_advanced), 
    ('PERC_IDENT',perc_ident), 
    ('PHI_PATTERN',phi_pattern), 
    ('PROGRAM',program), 
    #('PSSM',pssm), - It is possible to use PSI-BLAST via this API? 
    ('QUERY',sequence), 
    ('QUERY_FILE',query_file), 
    ('QUERY_BELIEVE_DEFLINE',query_believe_defline), 
    ('QUERY_FROM',query_from), 
    ('QUERY_TO',query_to), 
    #('RESULTS_FILE',...), - Can we use this parameter? 
    ('SEARCHSP_EFF',searchsp_eff), 
    ('SERVICE',service), 
    ('TEMPLATE_LENGTH',template_length),
    ('TEMPLATE_TYPE',template_type),
    ('THRESHOLD',threshold), 
    ('UNGAPPED_ALIGNMENT',ungapped_alignment), 
    ('WORD_SIZE',word_size), 
    ('CMD', 'Put')#, 
    #('BLAST_PROGRAM',blast_program) #Edited on April 29, 2012
    ] 
  query = [x for x in parameters if x[1] is not None] 
  message = urllib.urlencode(query) 

# Send off the initial query to qblast. 
# Note the NCBI do not currently impose a rate limit here, other 
# than the request not to make say 50 queries at once using multiple 
# threads. 
  request = urllib2.Request("http://blast.ncbi.nlm.nih.gov/Blast.cgi", 
                          message, 
                          {"User-Agent":"BiopythonClient"}) 
  handle = urllib2.urlopen(request) 

# Format the "Get" command, which gets the formatted results from qblast 
# Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node6.html on 9 July 2007     
  rid, rtoe = _parse_qblast_ref_page(handle) 
  parameters = [ 
    ('ALIGNMENTS',alignments), 
    ('ALIGNMENT_VIEW',alignment_view), 
    ('DESCRIPTIONS',descriptions), 
    ('ENTREZ_LINKS_NEW_WINDOW',entrez_links_new_window), 
    ('EXPECT_LOW',expect_low), 
    ('EXPECT_HIGH',expect_high), 
    ('FORMAT_ENTREZ_QUERY',format_entrez_query), 
    ('FORMAT_OBJECT',format_object), 
    ('FORMAT_TYPE',format_type), 
    ('NCBI_GI',ncbi_gi), 
    ('RID',rid), 
    ('RESULTS_FILE',results_file), 
    ('SERVICE',service), 
    ('SHOW_OVERVIEW',show_overview), 
    ('CMD', 'Get'), 
    ] 
  query = [x for x in parameters if x[1] is not None] 
  message = urllib.urlencode(query) 
  print "IT'S WORKING!!"
# Poll NCBI until the results are ready.  Use a 3 second wait 
  delay = 3.0 
  previous = time.time() 
  while True: 
      current = time.time() 
      wait = previous + delay - current 
      if wait > 0: 
          time.sleep(wait) 
          previous = current + wait 
      else: 
          previous = current 
  
      request = urllib2.Request("http://blast.ncbi.nlm.nih.gov/Blast.cgi", 
                                message, 
                                {"User-Agent":"BiopythonClient"}) 
      handle = urllib2.urlopen(request) 
      results = _as_string(handle.read()) 
  
      # Can see an "\n\n" page while results are in progress, 
      # if so just wait a bit longer... 
      if results=="\n\n": 
          continue 
    #   XML results don't have the Status tag when finished 
      if results.find("Status=") < 0: 
          break 
      i = results.index("Status=") 
      j = results.index("\n", i) 
      status = results[i+len("Status="):j].strip() 
      if status.upper() == "READY": 
          break 
  return StringIO(results)


def _parse_qblast_ref_page(handle): 
      """Extract a tuple of RID, RTOE from the 'please wait' page (PRIVATE). 
   
      The NCBI FAQ pages use TOE for 'Time of Execution', so RTOE is proably 
      'Request Time of Execution' and RID would be 'Request Identifier'. 
      """ 
      s = _as_string(handle.read()) 
      i = s.find("RID =") 
      if i == -1: 
          rid = None 
      else: 
          j = s.find("\n", i) 
          rid = s[i+len("RID ="):j].strip() 
   
      i = s.find("RTOE =") 
      if i == -1: 
          rtoe = None 
      else: 
          j = s.find("\n", i) 
          rtoe = s[i+len("RTOE ="):j].strip() 
   
      if not rid and not rtoe: 
          #Can we reliably extract the error message from the HTML page? 
          #e.g.  "Message ID#24 Error: Failed to read the Blast query: 
          #       Nucleotide FASTA provided for protein sequence" 
          #or    "Message ID#32 Error: Query contains no data: Query 
          #       contains no sequence data" 
          # 
          #This used to occur inside a <div class="error msInf"> entry: 
          i = s.find('<div class="error msInf">') 
          if i != -1: 
              msg = s[i+len('<div class="error msInf">'):].strip() 
              msg = msg.split("</div>",1)[0].split("\n",1)[0].strip() 
              if msg: 
                  raise ValueError("Error message from NCBI: %s" % msg) 
          #In spring 2010 the markup was like this: 
          i = s.find('<p class="error">') 
          if i != -1: 
              msg = s[i+len('<p class="error">'):].strip() 
              msg = msg.split("</p>",1)[0].split("\n",1)[0].strip() 
              if msg: 
                  raise ValueError("Error message from NCBI: %s" % msg) 
          #Generic search based on the way the error messages start: 
          i = s.find('Message ID#') 
          if i != -1: 
              #Break the message at the first HTML tag 
              msg = s[i:].split("<",1)[0].split("\n",1)[0].strip() 
              raise ValueError("Error message from NCBI: %s" % msg) 
          #We didn't recognise the error layout :( 
          #print s 
          raise ValueError("No RID and no RTOE found in the 'please wait' page, " 
                           "there was probably an error in your request but we " 
                           "could not extract a helpful error message.") 
      elif not rid: 
          #Can this happen? 
          raise ValueError("No RID found in the 'please wait' page." 
                           " (although RTOE = %s)" % repr(rtoe)) 
      elif not rtoe: 
          #Can this happen? 
          raise ValueError("No RTOE found in the 'please wait' page." 
                           " (although RID = %s)" % repr(rid)) 
   
      try: 
          return rid, int(rtoe) 
      except ValueError: 
          raise ValueError("A non-integer RTOE found in the 'please wait' page, %s" % repr(rtoe)) 

#Modified on, April 29, 2012 [by mamun]
Entrez.email = "mamun@cs.umanitoba.ca" ## replace 

## 
# get species mapping 
#

#def get_species_map (sepcies_file):
#  try: 
#    map = dict()
#    default = set() 
#    f = open ( species_file, 'r') 
#    species_list = [x.strip() for x in f]
#    for y in species_list:
#      if (y.startswith('#DEFAULT:')):
#        default.add( y.split(':')[1].strip() )
#	#comment on May 1, 2012
#	#map=y.split(':')[1].strip()
#      else:
#        map_info = y.split(':')
#        input = map_info[0]
#        output = set([ x.strip() for x in map_info[1].strip().split(',')])
#        output = output.union(default)
#        map[input] = list(output)
#    return map 
#  except IOError:
#    print species_file + ' not found'
#    sys.exit()

##
# get species in a list
#
#def get_species ( species_file ):
#  try: 
#    f = open ( species_file, 'r') 
#    species_list = [x.strip().split() for x in f]
#    return species_list 
#  except IOError:
#    print species_file + ' not found'
#    sys.exit()

## 
# ungap the seq
def ungap ( seq ) : 
  if ( seq.find('-') >= 0 ):  
    seq = Seq ( str(seq).replace('-',''))
  if ( seq.find('~') >= 0 ):  
    seq = Seq ( str(seq).replace('~',''))
  return seq 

##
# get id from title
# necessary since the blast alignments only return
# a title and not an id
# NOTE: searches specifically for WGS ids AAAAnnnnnnnn
def get_id_from_title ( title ):
  first_try = re.search ( '[A-Z]{4}\d{8}', title )
  if first_try:
    return str(first_try.group(0))
  else:
    # it BETTER be the inconsistently named D. melanogaster WGS  
	# I THINK, I NEED TO CHANGE SOMETHING HERE!!
    droso = re.search ( 'NW_\d{9}', title )  
    if droso:
      return str ( droso.group(0) )
    else: return None 


##
# get binomial name from title
# assumes that fasta record is in the form xx|xxxx|xx|ID|Genus species
def get_binom_from_title (title):
  if ( title.find('|') == -1):
    print 'expected full fasta header, did not find.'
    print 'assuming format is \">Genus species\"'  
    sp = title.strip().split()
    if (len (sp) == 1):
      print 'expected \">Genus species\", did not find.'
      print 'assuming \">Genus\" only.'
      return [ sp[0], None ]
    else:
      return sp[0:2]
  else:
    sp0 = title.strip().split('|')[4]
    sp = sp0.strip().split()
    return sp[0:2]
   
##
# find best blast hit
# takes a blast record and a genus and species
# returns None if nothing in record matches genus and species.
# o/w returns alignment 
def find_best_blast_hit ( blast_record , genus, species) :
  epsilon = 1e-30
  best_matched = None
  best_matched_seq = None
  best_score = sys.maxint 
  count = 1
  
  for (descr,align) in zip(blast_record.descriptions,blast_record.alignments):
    #print count,descr,align
    #count = count+1
    if (re.search ( genus + '( )*' + species , descr.title) and (not re.search ( 'chromosome', descr.title ))): 

      # matched the right genus / species. 
      # but maybe not the best.. check.
      if ( descr.e < best_score):
        best_score = descr.e
        best_matched = align 
      # test equality -- floating point style.
      elif ( fabs(descr.e - best_score) <= epsilon):
        #print 'must test length!'
         
        # get actual best sequence from db, if not already done.
        if not best_matched_seq:
          handle = Entrez.efetch ( db='nucleotide',id= get_id_from_title (best_matched.title) ,rettype='fasta')
          best_matched_seq = SeqIO.read(handle,format='fasta') 
	  handle.close()
       
        # get actual current sequence from db, 
        new_handle = Entrez.efetch ( db='nucleotide',id= get_id_from_title (descr.title) ,rettype='fasta')
        new_matched_seq = SeqIO.read(new_handle,format='fasta')
        new_handle.close()
        #print descr.title, best_matched.title
        #print len(new_matched_seq), len(best_matched_seq)
        if ( len(new_matched_seq) > len( best_matched_seq)) :
          best_matched_seq = new_matched_seq
          best_score = descr.e
          best_matched = align

  if (not best_matched_seq and best_matched):
    # we've made it to the end, found something, but not loaded the seq from the db yet. 
    handle = Entrez.efetch ( db='nucleotide',id= get_id_from_title (best_matched.title) ,rettype='fasta')
    best_matched_seq = SeqIO.read(handle,format='fasta') 
    handle.close()
  if (best_matched):
    return (best_matched_seq,best_matched)
  else:
    return None

def get_unique_file ( pref ) :
  unique = pref + ".fas"
  if os.path.isfile( pref + ".fas" ):
    count = 1
    unique = pref +'-'+ str(count) + ".fas" 
    while os.path.isfile ( unique ): 
      count = count + 1  
      unique = pref +'-'+ str(count) + ".fas" 
  return unique
  
if (len(sys.argv) < 3):
  print 'usage: python get_Record.py <input_file.fas> <species_list.txt>'
  sys.exit()

name = sys.argv[1]
print 'processing file ' + name
species_file = sys.argv[2]
print 'species file ' +  species_file
#edited by mamun on April 24, 2012
#species_dict = get_species_map( species_file )
#species_dict = ('Apis', 'Mellifera');
#print species_dict
master_set = set() 
master_aligns = set() 
try: 
  fasta = open (name)
except IOError:
  print name + ' not found'
  sys.exit()

for seq_record in SeqIO.parse(fasta,'fasta'):
  #seq_record=SeqIO.parse(fasta,'fasta')
  fullTitle =  seq_record.description
  #target_list = [ x  for x in species_dict.keys() if re.search(x,fullTitle) ]
  #if not (len(target_list) ==1 ): 
  #  print 'error, cannot uniquely find ' + fullTitle + ' in ' + str(species_dict.keys())
  #else:
    #target_key = target_list[0]
    # get file ready 
  save_file_name_pref = fullTitle.split('|')[0].replace(' ','') 
  save_file_name = get_unique_file( save_file_name_pref) 
  save_file = open (save_file_name,"w")
  current_species_list = []
  bad_string = ""
  print 'writing to file ' + save_file_name

    ## NOTE: E-value is hard-wired at 10^(-20) in the next line.
    #Commented on April 29, 2012
	#result_handle = _mike_qblast('blastn','discoMegaBlast','nt', ungap(seq_record.seq),expect=1e-20,hitlist_size=1000,megablast='true',template_length=18,template_type=0)
  result_handle = _mike_qblast('tblastx','nt', ungap(seq_record.seq))
  print seq_record.seq
  #print result_handle
  #text = result_handle.read()
  #print text
  read = NCBIXML.read(result_handle)
  
    #commented on May 1, 2012ttributeError:
    #for sp in species_dict[target_key]:
      ### get genus and species separately. 
    #genus_list = sp.split()
    #genus = genus_list[0]
    #species = genus_list[1]
  genus = "Apis"
  species = "Mellifera"
  best_hit =  find_best_blast_hit ( read, genus, species)
  print best_hit
  worst_hit = NCBIWWW.qblast("tblastx", "nr", "Apis Meliferra")
  print NCBIXML.read(worst_hit)
  if (best_hit):
    # we actually found something! 
    current_species_list.append ( best_hit[0])
    master_set.add ( best_hit[0])
    master_aligns.add (( fullTitle,  best_hit[1]))  
    print 'Consider it Done!'
  else:
    bad_string = bad_string + '>'+genus+' '+species+' not found\n'

  save_file.write ( bad_string) 
  SeqIO.write ( current_species_list , save_file, "fasta")  
  save_file.close()

# process master set
master_file_name = name.split('.')[0] + "_master_file.fas"
print 'writing to master file ' + master_file_name
master_file = open ( master_file_name,'w')
SeqIO.write ( list (master_set) , master_file, 'fasta') 
master_file.close()

# process master align
master_align_file_name = name.split('.')[0] + "_master_aligns.txt"
print 'writing to master alignment file ' + master_align_file_name
master_align = open (master_align_file_name,'w')
for x in master_aligns:
  for y in x[1].hsps:
    if ( y.expect < 1e-20):
      master_align.write ( 'match: ' + x[1].title + '\nvs.\n'  + x[0] + '\n' ) 
      master_align.write ( 'query position: ' + str ( y.query_start)  + ' subject position: ' + str ( y.sbjct_start ) + '\n') 
      if (y.frame and not (y.frame[0] == y.frame[1])):
        master_align.write ( '***REVERSE COMPLEMENT?***\n')
      master_align.write ( y.query +'\n') 
      master_align.write ( y.match +'\n')  
      master_align.write ( y.sbjct +'\n') 
      master_align.write ( '\n')
master_align.close()
