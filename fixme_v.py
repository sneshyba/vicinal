import vstuff as vs; reload(vs)
import numpy as np
import copy

# Get the spc slab (makes new xyzO, xyzH1, xyzH2, and shift-related arrays)
#filename = 'spc_4_4_6_v_orig.pdb'
#xbox=17.9629; ybox=29.3333; zbox=23.334523
#outfilename = 'spc_4_4_6_v.pdb'

#filename = 'spc_10_6_12_v_orig.pdb'
#filename = 'spc_10_6_12_v_origplus1.pdb'
#filename = 'spc_10_6_12_v_origplus2.pdb'
#filename = 'spc_10_6_12_v_origplus3.pdb'
#outfilename = 'spc_10_6_12_v_06131.pdb'
#xbox=44.90725; ybox=45.9565; zbox=44.6822104389

#filename = 'spc_10_6_12_vx_orig.pdb'
#outfilename = 'spc_10_6_12_vx06271.pdb'
#xbox=45.57594; ybox=45.98438; zbox=43.99999

filename = 'spc_20_6_12_vx_orig.pdb'
outfilename = 'spc_20_6_12_vx07081.pdb'
xbox=90.1507996879; ybox=46.4950169842; zbox=43.99999

# Specify which is the exposed surface, and load the slab
vicinaldir = 'y'; nycel=0
xyzO, xyzH1, xyzH2, shift, structure = \
vs.loaditnew(filename, xbox, ybox, zbox, vicinaldir)
slab = vs.slab(filename, structure, xyzO, xyzH1, xyzH2, xbox, ybox, zbox)

# Report dipole
print slab.getdipole()


# Get  nearest neighbor and defect information
nni,xyzshift = vs.getnni(xyzO,shift)
nnitol, nnltoi = vs.getnnitoletc(nni,xyzO,xyzH1,xyzH2,xyzshift)

# Check for any initial defects based on the original vicinal xyzO, xyzH1, 
# and xyzH2 arrays
nnitol, nnltoi, Ddefect, Adefect = \
vs.checkfordefects(nni,xyzO,xyzH1,xyzH2,xyzshift)

# Make new nnitol and,nnltoi arrays that correct for defects
nnitol_fixed, nnltoi_fixed = vs.fixitfast(\
slab, xyzO, xyzH1, xyzH2, xyzshift, nni, nnitol, nnltoi, 10000)

# Reconstruct the slab with defects fixed
xyzO_fixed, xyzH1_fixed, xyzH2_fixed = \
vs.reconstructall(xyzO, xyzH1, xyzH2, nni, nnitol_fixed, xyzshift)

# Save the good vicinal slab 
slab_v = vs.slab(outfilename,structure,xyzO_fixed, xyzH1_fixed, xyzH2_fixed)
slab_v.saveit()

# Report dipole
vdipole = slab_v.getdipole()
print vdipole, np.sqrt(np.sum(vdipole*vdipole))

