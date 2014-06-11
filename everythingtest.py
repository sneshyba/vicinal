import vstuff as vs; reload(vs)
import time

# Get the spc slab (makes new xyzO, xyzH1, xyzH2, and shift-related arrays)
filename = 'spc_4_4_6.pdb'; nx = 4; ny = 4; nz = 6
#filename = 'spc_4_4_2.pdb'; nx = 4; ny = 4; nz = 2
#filename = 'spc_10_6_12.pdb'; nx = 10; ny = 6; nz = 12
#filename = 'spc_10_6_14.pdb'; nx = 10; ny = 6; nz = 14

# Naming the output file
dum = filename.find('.pdb')
outfilename = filename[0:dum]+'_v.pdb'
badfilename = filename[0:dum]+'_v_orig.pdb'

# Specify which viscinal surface to generate, and load the slab
viscinaldir = 'y'; nycel=1
xyzO, xyzH1, xyzH2, viscinaldir, shift, vshift, xbox, ybox, zbox, structure = vs.loadit(filename, nx, ny, nz, viscinaldir, nycel)

# Get the nearest neighbor index and nearest-neighbor shift array
nni,xyzshift = vs.getnni(xyzO,shift) 

# Check for any initial defects based on the original viscinal xyzO, xyzH1, and xyzH2 arrays
nnitol, nnltoi, Ddefect, Adefect = vs.checkfordefects(nni,xyzO,xyzH1,xyzH2,xyzshift)

#Output bad rotated slab
