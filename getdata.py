#!/usr/bin/python
import io

# Read in file cfl3d.prout
with io.open('cfl3d.prout','r') as file:
   lines = file.readlines()

# find the line that starts with "BLOCK" because this has the dimension
# and then read the jdim
for i in range( len(lines) ):
   words = lines[i].split()
   if len(words) > 0:
      if words[0] == "BLOCK" :
         jdim = words[6]
         j = i
         break
print 'jdim: ' + jdim

# Skip the first block of information, starts at the top of the block then
# loops until the first blank line
for i in range( j+3, len(lines) ):
   if len(lines[i]) == 1:
      k = i
      break

# Now write out the second block to some file
with io.open('boundary.dat','w') as file:
   for i in range( k+1, len(lines) ):
      
      if len(lines[i]) == 1:
         break
      
      if i == k+1:
         file.write( 'variables =' + lines[i] )
      
      else:
         file.write(lines[i])
      

   
