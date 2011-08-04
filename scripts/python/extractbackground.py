
import glob
import os
import random
import bisect

nbtoextract=10000
seqlength=2000
limittry=nbtoextract*5

aligncoords=[]
weights=[]
for folder in glob.glob('[0-9]*/'):
   os.chdir(folder)
   with open("mavid.mfa") as fmavid:
      firstline=fmavid.readline().rstrip().split(" ")
      coord=[firstline[1],int(firstline[2]),int(firstline[3])]
      aligncoords.append(coord)
      weights.append(float(coord[2]-coord[1]))
   os.chdir("..")
print "Alignment coordinates imported" 

sumweights=0.0
for w in weights:
   sumweights+=w
for i in range(len(weights)):
      weights[i]/=sumweights

## Open file containing coding sequence coordinates
fcds=open("CDS-coord.dat")
cds=[i.rstrip().split(" ") for i in fcds.readlines()]
for i in cds:
   i[1]=int(i[1])
   i[2]=int(i[2])
fcds.close()
print "CDS coordinates imported" 



def weighted_choice(coords,weights):
   x=random.random()
   for i,w in enumerate(weights):
      if x<w:
         return coords[i];
      x -= w

def isnotin(coord,listcoords):
   for i in listcoords:
      if coord[0]==i[0] and i[1]<coord[2] and i[2]>coord[1]:
         return False
   return True

class ToManyTriesError(Exception):
   def __init__(self, value):
       self.value = value
   def __str__(self):
       return repr(self.value)

## Record selection functions
def posrecordleft(coordlist,coord):
   pos=bisect.bisect_left(coordlist,coord)
   if pos>0:
      pos=pos-1
   return pos

coords=[]
nbtry=0
nbextracted=0
while nbextracted<nbtoextract:
   align=weighted_choice(aligncoords,weights)
   if align[2]-align[1]-seqlength>0:
      start=random.randint(0,align[2]-align[1]-seqlength)
      potcoord=[align[0],align[1]+start,align[1]+start+2000]

      alignstart=[align[0],align[1],align[1]]
      alignstop=[align[0],align[2],align[2]]
      posstart=posrecordleft(cds,alignstart)
      posstop=posrecordleft(cds,alignstop)
      cdsinalign=[]
      for i in range(posstart,posstop):
         cdsinalign.append(cds[i])

      if isnotin(potcoord,cdsinalign) and isnotin(potcoord,coords):
         coords.append(potcoord)
         nbextracted+=1
   nbtry+=1
   if nbtry>limittry:
      print "error"
      exit()
#      raise ToManyTriesError

print coords
