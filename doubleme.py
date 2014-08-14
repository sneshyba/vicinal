import vstuff as vs; reload(vs)
import numpy as np
import copy

xcel = 4.4907312
ycel = 7.7781746
zcel = 3.6666666

# vmd: set cell [pbc set {44.907312 45.9843789846 43.9999992} -all]; pbc box ... tilt angle is 9.826 deg
filename = 'spc_10_6_12.pdb' 
xbox=xcel*10; ybox=ycel*6; zbox=zcel*12
outfilename = 'spc_20_6_12.pdb'

# open and read the input file
fin = open(filename, "r")
data = fin.readlines()
fin.close()

# open the output file and write the initial data
fout = open(outfilename, "w")
for i in range(len(data)-1):
    fout.write(data[i])
    pass

# loop over the data and count the atoms
natoms = 0
for i in range(len(data)):
    line = data[i]
    splitline = line.split()
    if splitline[0] == 'ATOM':
        natoms += 1
nresidues = natoms/3
print natoms, nresidues        

# loop over the data again and write out with shift
space = ' '
for i in range(len(data)):
    line = data[i]
    splitline = line.split()
    if splitline[0] == 'ATOM':
        natoms += 1
        #splitline[1] = str(natoms)
        if np.mod(natoms,3) == 1:
            nresidues += 1
        #splitline[4] = str(nresidues)
        #splitline[5] = str(float(splitline[5])+xbox)
        #dum = space.join(splitline)+'\n'
        dum = "%4s %6i %3s %4s %5i %11.3f %7.3f %7.3f \n" \
        %('ATOM',natoms,splitline[2],splitline[3],nresidues, \
        float(splitline[5])+xbox, float(splitline[6]), float(splitline[7]))
        fout.write(dum)
    elif splitline[0]=='REMARK':
        pass
    else:
        fout.write(line) 
   
fout.close()
print xbox*2, ybox, zbox