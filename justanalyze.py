import vstuff as vs; reload(vs)
import numpy as np
import copy

# Get the spc slab (makes new xyzO, xyzH1, xyzH2, and shift-related arrays)
filename = 'spc_4_4_6.pdb'; nx = 4; ny = 4; nz = 6
#filename = 'spc_4_4_2.pdb'; nx = 4; ny = 4; nz = 2
#filename = 'spc_10_6_12.pdb'; nx = 10; ny = 6; nz = 12
#filename = 'spc_10_6_14.pdb'; nx = 10; ny = 6; nz = 14

# Naming the output file
dum = filename.find('.pdb')
surfacefilename = filename[0:dum]+'_surface.pdb'
bulkfilename = filename[0:dum]+'_bulk.pdb'

# Not going to do anything to this but add space above the y-facet
viscinaldir = 'y'; nycel=0
xyzO, xyzH1, xyzH2, viscinaldir, shift, vshift, xbox, ybox, zbox, structure = vs.loadit(filename, nx, ny, nz, viscinaldir, nycel)
NR = len(xyzO)

# Get the nearest neighbor index and nearest-neighbor shift array
nni,xyzshift = vs.getnni(xyzO,shift)

# Check for any initial defects based on the original viscinal xyzO, xyzH1, and xyzH2 arrays
nnitol, nnltoi, Ddefect, Adefect = vs.checkfordefects(nni,xyzO,xyzH1,xyzH2,xyzshift)

# Get residues at the surface & bulk
surfacelist, dum = np.where(nni==-1)
alllist = range(NR)
bulklist = np.squeeze(np.argwhere(~np.in1d(alllist,surfacelist)))

# Get slabs of bulk and surface residues
slab_surface = vs.slab(surfacefilename,structure,xyzO, xyzH1, xyzH2)
slab_surface.cullout(bulklist)
slab_surface.saveit()

slab_bulk = vs.slab(bulkfilename,structure,xyzO, xyzH1, xyzH2)
slab_bulk.cullout(surfacelist)
slab_bulk.saveit()


'''
# Reconstruct & rotate a viscinal slab with defects fixed
xyzO_fixed, xyzH1_fixed, xyzH2_fixed = vs.reconstructit(xyzO, xyzH1, xyzH2, nni, nnitol_fixed, xyzshift)
xyzO_rot, xyzH1_rot, xyzH2_rot, xboxp, yboxp, zboxp = vs.rotateit(xyzO_fixed, xyzH1_fixed, xyzH2_fixed, viscinaldir, shift, vshift, xbox, ybox, zbox)

# Save the good vicinal slab 
slab_v = vs.slab(outfilename,structure,xyzO_rot, xyzH1_rot, xyzH2_rot)
slab_v.saveit()

'''