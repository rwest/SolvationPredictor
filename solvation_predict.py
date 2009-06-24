#!/usr/bin/env python
# encoding: utf-8
"""
solvation_predict.py to run RMG Abraham estimator for many solvent/solute combinations

Reads in three files:
  adjList.txt containing RMG-style adjacency lists of all the solutes by name
  Solvent Database.csv containing Abraham's parameters for the solvents (with a header row)
  Solvent-Solute List.csv containing a list of Solvent/Solute pairs that should be estimated (with a header row)
Creates one file:
  Solvent-Solute-Solvation_by_RMG.csv
  
Requires a working RMG Abraham.class (which may be inside a .jar)
 which in turn needs functioning GATPFit.exe in the appropriate path (see RMG instructions)

Created by Richard West on 2009-06-23.
Copyright (c) 2009 MIT. All rights reserved.
"""


# I think all these are standard libraries, so you should have them all (in Python 2.5+)
import sys
import os
import csv
import subprocess
import re

# you may need to change this to point to where your Abraham.class file is (which may be in an RMG.jar)
RMG_classpath = '%rmg%/classes'


# read in adjacency list data for all solutes
adjDict = dict()
adjList = open('adjList.txt','U') # the 'U' allows for universal line endings (ie. windows stuff)
for line in adjList:
    line=line.rstrip()+'\n' # strip the line ending and replace with a proper one
    if line.startswith('InChI='): # new solute; reset adj
        inchi = line.strip()
        adj = ''
    adj += line 
    if line.isspace(): # reached end of adj; store it
        adjDict[inchi] = adj
# should now have object like adjDict['InChI=1/C2H4O2/c1-4-2-3/h2H,1H3'] = 'multi-line adjacency list string'


# read in solvent data
solventReader = csv.DictReader(open('Solvent Database.csv'), dialect='excel') # assumes first line contains field names
solventDict = dict()
for solvent in solventReader:
    assert solvent.has_key('a'), "Solvent Database should have 'a' column"
    solventName = solvent['Solvents']
    solventDict[solventName] = solvent
# should now have object like: solventDict['water']['a']


# read in Solvent/Solute list
solventSoluteList = list() # used to make sure the output is in the same order as the input
solventSoluteReader = csv.DictReader(open('Solvent-Solute list.csv'), dialect='excel') # assumes first line contains field names
solutesDict = dict() # a list of solutes for each solvent
dataDict = dict() # a dictionary (of dictionaries) that we will store values in later
for entry in solventSoluteReader:
    solvent = entry['Solvent']
    solute = entry['Solute']
    if not solutesDict.has_key(solvent):
        solutesDict[solvent]=list()
    if not dataDict.has_key(solvent):
        dataDict[solvent]=dict()
    solutesDict[solvent].append(solute)
    solventSoluteList.append((solvent,solute))  
    
# should now have structure like: solutesDict['water'] = [list of solutes]
#                           and:  dataDict['water'] = {} # empty dictionary to store values in later
#                           and:  solventSoluteList = [list of (solvent,solute) tuples]
    


errors=[]  # we'll store them all up to summarise at the end

# now iterate over solvents 
for solvent,solutes in solutesDict.items():  # (warning: iterating over a dictionary, so order may be different each time)
    sv=solventDict[solvent] # get solvent parameters
    AbrahamFile=open('Abraham_input.txt','w')
    AbrahamFile.write('SolventParameters: %s %s %s %s %s %s   %s\n\n'%(sv['c'], sv['e'], sv['s'], sv['a'], sv['b'], sv['l'], sv['Solvents']))
    for solute in solutes:
        AbrahamFile.write(adjDict[solute]) # write the adjacency lists out to the file
    AbrahamFile.close()

    # run the java program 
    command = 'java -classpath %s Abraham'%RMG_classpath
    os.path.isdir('GATPFit') or os.makedirs('GATPFit') # apparently this dir must exist for the java to run
    print "Running java Abraham for ",solvent
    process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True   )
    output,error=process.communicate()
    if error:
        errorMessage= "Error for something(s) dissolved in %s:\n%s"%(solvent,error)
        print errorMessage
        errors.append(errorMessage)
        # raise StandardError(error) # if we want it to stop on error, uncomment this line
        
#    ## interpreting the standard output
#    # not sure this works as expected...
#    searcher=re.compile('The free energy of solvation in (?P<solvent>\S+) at 298K without reference state corrections is  = (?P<correction>\S+) J/mol')
#    matches = searcher.findall(output)
#    assert len(matches) == len(solutes), "%d answers from %d inputs to java Abraham code"%(len(matches),len(solutes))
#    for i, (solvent, correction) in enumerate(matches):
#        solute=solutes[i]
#        print solvent,solute,correction

        
    ## reading the output file
    # I think this works better
    AbrahamFile=open('Abraham_output.txt','U') # U to allow universal line endings (for Windows)
    line=AbrahamFile.readline()
    match = re.search('Solvent:\s+(?P<solvent>.*)',line)
    assert match, "Didn't find Solvent: line at top of output file"
    assert solvent.startswith( match.groupdict()['solvent'].strip()), "Solvent in output is not solvent in input!"
       # use startswith() because when solvent is 'dibutyl ether' the Abraham_output.txt will just say 'dibutyl'
    
    searcher=re.compile('(?P<solute>InChI=\S+)\s+(?P<value>\S+)')
    matches = searcher.findall(AbrahamFile.read()) # read the rest in one go
    
    for (solute, value) in matches:
        print solvent,solute,value
        dataDict[solvent][solute]=value 


# create the output file 
outfile = open('Solvent-Solute-Solvation_by_RMG.csv','wb')
outputter = csv.DictWriter(outfile,['solvent','solute','logK'],dialect='excel')
for (solvent,solute) in solventSoluteList: # the order of the input file is preserved
    answer={'solvent':solvent, 'solute':solute, 'logK':dataDict[solvent][solute]}
    outputter.writerow(answer)
outfile.close()

print "\nSummary of errors:"
for error in errors:
    print error
        
def main():
    pass
            
if __name__ == '__main__':
    main()

