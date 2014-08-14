import vstuff as vs; reload(vs)

# Get the spc slab (makes new xyzO, xyzH1, xyzH2, and shift-related arrays)
#filename = 'spc_4_4_6.pdb'; nx = 4; ny = 4; nz = 6 
#filename = 'spc_4_4_2.pdb'; nx = 4; ny = 4; nz = 2
#filename = 'spc_10_6_12.pdb'; nx = 10; ny = 6; nz = 12
#filename = 'spc_10_6_14.pdb'; nx = 10; ny = 6; nz = 14
filename = 'spc_20_6_12.pdb'; nx = 20; ny = 6; nz = 12

# Naming the output file
dum = filename.find('.pdb')
outfilename = filename[0:dum]+'_vx_orig.pdb'

# Specify which vicinal surface to generate, and load the slab
vicinaldir = 'yx'; nycel=1
xyzO, xyzH1, xyzH2, vicinaldir, shift, vshift, xbox, ybox, zbox, structure = vs.loadit(filename, nx, ny, nz, vicinaldir, nycel)
slab = vs.slab(filename,structure,xyzO, xyzH1, xyzH2, xbox, ybox, zbox)

# Reconstruct & rotate it, then save it
slab.rotateit(vicinaldir,vshift)
slab.saveit(outfilename)
print slab.xbox, slab.ybox, slab.zbox

