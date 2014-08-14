import vstuff as vs; reload(vs)
import numpy as np

xcel = 4.4907312
ycel = 7.7781746
zcel = 3.6666666

# vmd: set cell [pbc set {44.907312 45.9843789846 43.9999992} -all]; pbc box ... tilt angle is 9.826 deg
filename = 'spc_10_6_12.pdb' 
xbox=xcel*10; ybox=ycel*6; zbox=zcel*12
namestem = filename.find('.pdb')

# Specify which is the exposed surface, and load the slab
xyzO, xyzH1, xyzH2, shift, structure = vs.loaditnew(filename, xbox, ybox, zbox, 'y')

# Shift in x
#xyzO[:,0]+=xbox
#xyzH1[:,0]+=xbox
#xyzH2[:,0]+=xbox

# Save it
slab = vs.slab(filename, structure, xyzO, xyzH1, xyzH2, xbox, ybox, zbox)
slab.filename = filename[0:namestem]+'_shifted0x.pdb'
slab.saveit()

