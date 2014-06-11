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
testfilename = filename[0:dum]+'_shift.pdb'

# Specify which viscinal surface to generate, and load the slab
viscinaldir = 'y'; nycel=1
xyzO, xyzH1, xyzH2, viscinaldir, shift, vshift, xbox, ybox, zbox, structure = vs.loadit(filename, nx, ny, nz, viscinaldir, nycel)

# Save the original slab 
slab = vs.slab(testfilename,structure,xyzO, xyzH1, xyzH2)
slab.saveit()


# Get the nearest neighbor index and nearest-neighbor shift array
nni,xyzshift = vs.getnni(xyzO,shift) 

# Check for any initial defects based on the original viscinal xyzO, xyzH1, and xyzH2 arrays
nnitol, nnltoi, Ddefect, Adefect = vs.checkfordefects(nni,xyzO,xyzH1,xyzH2,xyzshift)

# Make new nnitol and,nnltoi arrays that correct for defects
nnitol_fixed, nnltoi_fixed = vs.fixit(nni, nnitol, nnltoi, 3000)

# Reconstruct & rotate a viscinal slab with defects fixed
xyzO_fixed, xyzH1_fixed, xyzH2_fixed = vs.reconstructit(xyzO, xyzH1, xyzH2, nni, nnitol_fixed, xyzshift)
xyzO_rot, xyzH1_rot, xyzH2_rot, xboxp, yboxp, zboxp = vs.rotateit(xyzO_fixed, xyzH1_fixed, xyzH2_fixed, viscinaldir, shift, vshift, xbox, ybox, zbox)

# Save the good vicinal slab 
slab_v = vs.slab(outfilename,structure,xyzO_rot, xyzH1_rot, xyzH2_rot)
slab_v.saveit()



# Rotate a viscinal slab that has the original defects
xyzO_rot_orig, xyzH1_rot_orig, xyzH2_rot_orig, xboxp, yboxp, zboxp = vs.rotateit(xyzO, xyzH1, xyzH2, viscinaldir, shift, vshift, xbox, ybox, zbox)

# Save it
slab_vorig = vs.slab(badfilename,structure,xyzO_rot_orig, xyzH1_rot_orig, xyzH2_rot_orig)
slab_vorig.saveit()
#vs.saveit(badfilename,structure,xyzO_rot_orig, xyzH1_rot_orig, xyzH2_rot_orig)
#slab2 = vs.slab(badfilename,structure,xyzO_rot_orig, xyzH1_rot_orig, xyzH2_rot_orig)
#slab2.saveit()


#tic = time.time()
#toc = time.time()
#print (toc-tic)/60, "minutes"
#execfile("plotit.py")
