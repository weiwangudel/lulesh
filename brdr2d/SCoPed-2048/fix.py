import sys                                                                       
import string
import re                                                                        
                                                                                 
f = open(sys.argv[1], 'rU')                                                      
current = 0                                                                      
fileStr = []
in_the_block = 0
the_block = []
num_parallel = 0
applied = 0
for line in f:   ## iterates over the lines of the file                          
  #prepare insert 
  if (num_parallel == 1 and applied < 2):
    if (string.find(line,"for") == -1):
      for i in the_block:
        fileStr += i,
      applied += 1
      the_block = []
      num_parallel = 0

  if (in_the_block == 0):
    if (string.find(line,"ax1;") != -1):
      in_the_block = 1
      the_block += line,
    else:
      fileStr +=  line,
  else:
    if (string.find(line,"int jj;") != -1):
      in_the_block = 0 
    the_block += line,

  if (string.find(line,"omp parallel") != -1):
    num_parallel += 1 



f.close() 

for line in fileStr:
  print line,
