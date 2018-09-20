#!/usr/bin/python
import math
import sys
from distutils.version import LooseVersion, StrictVersion

###### BEGIN CLASS BLOCK ######
# block size is the full size of the block as given in 
# cfl3d.prout but the points being kept are for the lower 
# block boundaries and so range from 2 to (size-1)
class Block:
    
    def __init__(self, idnum, loc, size):
        self.id  = idnum
        self.loc = loc
        self.size = size
        self.g1 = []
        self.g2 = []
        self.xpts = []
        self.zpts = []
        

    def addpt1(self,p):
        self.g1.append(p)
        w = p.split()
        self.xpts.append(float(w[3]))
        self.zpts.append(float(w[5]))
    
    def addpt2(self,p):
        self.g2.append(p)
    
    def writeS(self):
        print "Block " + str(self.id) + " has size: " 
        print self.size

    def writePts(self,f,jst=0):
        if not f.closed:
            if jst == 0 :
               for i in range( self.size[1]-2 ):
                   f.write(self.g2[i])
            else:
                w = self.g2[0].split() # split into words
                js = int(w[1]) # get the first j for the block
                js = jst-js+1 # shift the starting j value 
                #print js
                for i in range( self.size[1]-2 ):
                    w = self.g2[i].split()
                    s1 = " " + w[1] + " "
                    s2 = " " + str(js+int(w[1])) + " "
                    li = self.g2[i].replace(s1,s2)
                    f.write(li)
                

###### END CLASS BLOCK ######


verbose = True
if len(sys.argv) > 1:
   if len(sys.argv) == 2:
      if sys.argv[1] == 'quiet' or sys.argv[1] == 'q':
         verbose = False
   else:
      print "only one argument is accepted (q or quiet for no output)"

# get python version, I don't know if it is backwards
# compatible but whatev
pyversion = sys.version.split()[0]
if verbose:
   print(sys.version)

# Read in file cfl3d.prout
lines=[]
with open('cfl3d.prout','r') as file:
    lt = file.readlines()
    # check for unicode since this has been a problem before
    if isinstance(lt[0], unicode):
        for i in range( len(lt) ):
            lines.append(lt[i][:-2].encode("ascii"))
    else:
        lines = lt
        del lt


# find the line that starts with "BLOCK" because this has the dimension
# and create a new block using these dimensions
nbl = 0
blocks = []
varLine = ""
for i in range( len(lines) ):
    if  lines[i].find("BLOCK") != -1  :
        words = lines[i].split()
        #print words
        size = []
        size.append(int(words[-3]))
        size.append(int(words[-2]))
        size.append(int(words[-1]))
        nbl = nbl+1
        blocks.append(Block(nbl,i,size))
        #j = i + size[1] + 6 
        j = i + 4

        # print str(i) + " " + str(j)

        for k in range(j, j+size[1]):
            blocks[nbl-1].addpt1(lines[k])
            i = k
    
        j = i + 3
 
        if varLine == "":
            varLine = lines[j-1]
        
        for k in range(j, j+size[1]-2):
            blocks[nbl-1].addpt2(lines[k])
            i = k

if verbose:
   if LooseVersion(pyversion) < LooseVersion("3.1"):
      print "{0} Blocks total".format(nbl)
   else:
      print "{} Blocks total".format(nbl)

# lists to hold the lowest blocks and their id numbers
bottom = []
bottomIds = []

# Find first block, first point lies in the lowest x and lowest z
j=0
for i in range(nbl):
    if blocks[i].xpts[0] <= blocks[j].xpts[0]:
        if blocks[i].zpts[0] < blocks[j].zpts[0]:
            j=i
            
bottom.append(blocks[j])
bottomIds.append(blocks[j].id)

# Now fill out the bottom row one at a time. Find the next lego 
# piece
it = 0
while 2>1 and it<nbl:
    it = it+1 # just in case
    cdist = 1e6
    j = 0
    for i in range(nbl):
        # skip if this block is already in my list
        if blocks[i].id in bottomIds:
            continue
    
        # Distance between current last block and the trial block
        # keep the shortest so far
        dist = math.pow( (bottom[-1].xpts[-1] - blocks[i].xpts[0]) ,2)
        dist = dist + math.pow( (bottom[-1].zpts[-1] - blocks[i].zpts[0]), 2)
        dist = math.sqrt(dist)
    
        #print str(dist) + " " + str(cdist)
    
        if dist < cdist:
            cdist = dist
            j = i
            #print str(bottom[-1].xpts[-1]) + " " + str(blocks[i].xpts[0])
            #print str(bottom[-1].zpts[-1]) + " " + str(blocks[i].zpts[0])
            #print

    #print str(cdist) + " " + str(blocks[j].id)

    if cdist > .001 :
        break
    
    bottom.append(blocks[j])
    bottomIds.append(blocks[j].id)

if verbose:   
   if LooseVersion(pyversion) < LooseVersion("3.1"):
      print "lower boundary has {0} blocks".format(len(bottom))
   else:
      print "lower boundary has {} blocks".format(len(bottom))
#for i in range(len(bottom)):
#    print ( str(bottom[i].id) + 
#            " " + str(bottom[i].xpts[0]) + 
#            " " + str(bottom[i].zpts[0]) +
#            " " + str(bottom[i].xpts[-1]) + 
#            " " + str(bottom[i].zpts[-1]) )

with open("boundary.dat","w") as file:
    file.write(varLine)
    js = 1
    for i in range(len(bottom)):
        if i > 0:
            js = js + bottom[i-1].size[1]-2
        bottom[i].writePts(file,js)
