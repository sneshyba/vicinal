import vstuff as vs; reload(vs)
import numpy as np
import copy

# Get the spc slab (makes new xyzO, xyzH1, xyzH2, and shift-related arrays)
#filename = 'spc_4_4_6_v_withdefects.pdb'; xbox=17.9629248; ybox=29.3333332974; zbox=23.3345234043

#filename = 'spc_4_4_6_v_withdefects_fixed.pdb'; xbox=17.9629248; ybox=29.3333332974; zbox=23.3345234043
# vmd: set cell [pbc set {17.9629248 29.333333 23.3345234} -all]; pbc box

filename = 'spc_4_4_6_v5_withdefects_fixed.pdb'; xbox=17.9629248; ybox=30.6376674579 ; zbox=22.3411052195
# vmd: set cell [pbc set {17.9629248 30.6376674579 22.3411052195} -all]; pbc box
namestem = filename.find('.pdb')

# Specify which is the exposed surface, and load the slab
viscinaldir = 'y'; nycel=0
xyzO, xyzH1, xyzH2, shift, structure = vs.loaditnew(filename, xbox, ybox, zbox, viscinaldir)
slab = vs.slab(filename, structure, xyzO, xyzH1, xyzH2, xbox, ybox, zbox)

# Get  nearest neighbor and defect information
nni,xyzshift = vs.getnni(xyzO,shift)
nnitol, nnltoi = vs.getnnitoletc(nni,xyzO,xyzH1,xyzH2,xyzshift)
Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
print "Number of defects = ", len(Ddefect), len(Adefect)

# Save the surface and bulk slab subsets
surfacelist, dum = np.where(nni==-1)
slab_surface = copy.deepcopy(slab)
slab_surface.filename = filename[0:namestem]+'_surface.pdb'
slab_surface.keep(surfacelist)
slab_surface.saveit()
slab_bulk = copy.deepcopy(slab)
slab_bulk.filename = filename[0:namestem]+'_bulk.pdb'
slab_bulk.cullout(surfacelist)
slab_bulk.saveit()

# Save Ddefects
Ddefectlist = np.unique(Ddefect[:,0])
slab_Ddefect = copy.deepcopy(slab)
slab_Ddefect.filename = filename[0:namestem]+'_Ddefect.pdb'
slab_Ddefect.keep(Ddefectlist)
slab_Ddefect.saveit()

# Save Adefects
Adefectlist = np.unique(Adefect[:,0])
slab_Adefect = copy.deepcopy(slab)
slab_Adefect.filename = filename[0:namestem]+'_Adefect.pdb'
slab_Adefect.keep(Adefectlist)
slab_Adefect.saveit()

# Save instances where there are three -1's in nni
a,b = vs.count_unique(surfacelist)
a3list = a[np.argwhere(b==3)]
slab_a3 = copy.deepcopy(slab)
slab_a3.filename = filename[0:namestem]+'_a3.pdb'
slab_a3.keep(a3list)
slab_a3.saveit()
