import vstuff as vs; reload(vs)
import numpy as np
import copy

# Get the spc slab (makes new xyzO, xyzH1, xyzH2, and shift-related arrays)
#filename = 'spc_4_4_6_v_withdefects.pdb'; xbox=17.9629248; ybox=29.3333332974; zbox=23.3345234043
#namestem = filename.find('.pdb')

filename = 'spc_4_4_6_v5_withdefects.pdb'; xbox=17.9629248; ybox=30.6376674579 ; zbox=22.3411052195
namestem = filename.find('.pdb')

# Specify which is the exposed surface, and load the slab
viscinaldir = 'y'; nycel=0
xyzO, xyzH1, xyzH2, shift, structure = vs.loaditnew(filename, xbox, ybox, zbox, viscinaldir)
slab = vs.slab(filename, structure, xyzO, xyzH1, xyzH2, xbox, ybox, zbox)

# Get  nearest neighbor and defect information
nni,xyzshift = vs.getnni(xyzO,shift)
nnitol, nnltoi = vs.getnnitoletc(nni,xyzO,xyzH1,xyzH2,xyzshift)

# Check for any initial defects based on the original viscinal xyzO, xyzH1, and xyzH2 arrays
nnitol, nnltoi, Ddefect, Adefect = vs.checkfordefects(nni,xyzO,xyzH1,xyzH2,xyzshift)

# Make new nnitol and,nnltoi arrays that correct for defects
nnitol_fixed, nnltoi_fixed = vs.fixit(nni, nnitol, nnltoi, 500)

# Reconstruct & rotate a viscinal slab with defects fixed
xyzO_fixed, xyzH1_fixed, xyzH2_fixed = vs.reconstructit(xyzO, xyzH1, xyzH2, nni, nnitol_fixed, xyzshift)

# Save the good vicinal slab 
outfilename = filename[0:namestem]+'_fixed.pdb'
slab_v = vs.slab(outfilename,structure,xyzO_fixed, xyzH1_fixed, xyzH2_fixed)
slab_v.saveit()
