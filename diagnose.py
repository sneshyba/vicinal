import vstuff as vs; reload(vs)
import time

# Get the spc slab (makes new xyzO, xyzH1, xyzH2, and shift-related arrays)
#filename = 'spc_4_4_6.pdb'; nx = 4; ny = 4; nz = 6
filename = 'spc_4_4_2.pdb'; nx = 4; ny = 4; nz = 2
#filename = 'spc_10_6_12.pdb'; nx = 10; ny = 6; nz = 12
#filename = 'spc_10_6_14.pdb'; nx = 10; ny = 6; nz = 14

# Naming the output file
dum = filename.find('.pdb')
surfacefilename = filename[0:dum]+'_surface.pdb'
Ddefectfilename = filename[0:dum]+'_Ddefects.pdb'
Adefectfilename = filename[0:dum]+'_Adefects.pdb'

# Specify which viscinal surface to generate, and load the slab
viscinaldir = 'y'; nycel=1
xyzO, xyzH1, xyzH2, viscinaldir, shift, vshift, xbox, ybox, zbox, structure = vs.loadit(filename, nx, ny, nz, viscinaldir, nycel)

# Get the nearest neighbor index and nearest-neighbor shift array
nni,xyzshift = vs.getnni(xyzO,shift) 

# Check for any initial defects based on the original viscinal xyzO, xyzH1, and xyzH2 arrays
nnitol, nnltoi, Ddefect, Adefect = vs.checkfordefects(nni,xyzO,xyzH1,xyzH2,xyzshift)

# Identify surface defects
Ddefectsurf, Adefectsurf = vs.idsurfacedefects(nni, nnitol, nnltoi)


#Ddefectsurface, Adefectsurface = vs.idsurfacedefects(nni,nnltoi,nnitol)


# Isolate Ddefects with no -1 in nni

# Isolate Adefects with no -1 in nni

# Reconstruct, rotate, & save a viscinal slab with defects fixed
#xyzO_fixed, xyzH1_fixed, xyzH2_fixed = vs.reconstructit(xyzO, xyzH1, xyzH2, nni, nnitol_fixed, xyzshift)
#xyzO_rot, xyzH1_rot, xyzH2_rot, xboxp, yboxp, zboxp = vs.rotateit(xyzO_fixed, xyzH1_fixed, xyzH2_fixed, viscinaldir, shift, vshift, xbox, ybox, zbox)
#vs.saveit(outfilename,structure,xyzO_rot, xyzH1_rot, xyzH2_rot)

# Rotate & save a viscinal slab that has the original defects
#xyzO_rot_orig, xyzH1_rot_orig, xyzH2_rot_orig, xboxp, yboxp, zboxp = vs.rotateit(xyzO, xyzH1, xyzH2, viscinaldir, shift, vshift, xbox, ybox, zbox)
#vs.saveit(badfilename,structure,xyzO_rot_orig, xyzH1_rot_orig, xyzH2_rot_orig)

#tic = time.time()
#toc = time.time()
#print (toc-tic)/60, "minutes"
#execfile("plotit.py")
