import numpy as np
import copy
import Bio
import Bio.PDB 
import pdb

class slab:
    def __init__(self,filename,structure,xyzO, xyzH1, xyzH2, xbox=0, ybox=0, zbox=0):
        self.xyzO = copy.deepcopy(xyzO)
        self.xyzH1 = copy.deepcopy(xyzH1)
        self.xyzH2 = copy.deepcopy(xyzH2)
        self.filename = copy.deepcopy(filename)
        self.structure = copy.deepcopy(structure)
        self.xbox = xbox
        self.ybox = ybox
        self.zbox = zbox
        self.NR = len(xyzO)


    def saveit(self,filename=None):
        if filename == None:
            filename = self.filename
        print 'saving ', filename
        structure = self.structure 
        xyzO = self.xyzO
        xyzH1 = self.xyzH1
        xyzH2 = self.xyzH2

        # Set the pdb structure
        IO = Bio.PDB.PDBIO()
        IO.set_structure(structure)

        # Get all the coordinates
        i = -1 # residue counter
        for model in structure:
            for chain in model:
                j=0 # atom counter
                for residue in chain:
                        i = i+1 # residue counter
                        for atom in residue:
                            test=np.mod(j,3)
                    	    if (test == 0):
                                temp = xyzO[i]
                            elif (test ==1):
                                temp = xyzH1[i]
                            elif (test == 2):
                       	        temp = xyzH2[i]
                            atom.set_coord(temp)
                            j = j+1
        print len(chain)
        IO.save(filename)


    def cullout(self,badlist):
        structure = self.structure 
        xyzO = self.xyzO
        xyzH1 = self.xyzH1
        xyzH2 = self.xyzH2

        # Set the pdb structure
        IO = Bio.PDB.PDBIO()
        IO.set_structure(structure)

        # Get all the coordinates
        i = -1 # residue counter
        NR = 0
        for model in structure:
            for chain in model:
                oldchain = copy.deepcopy(chain)
                j=0 # atom counter
                for residue in oldchain:
                    i = i+1 # residue counter
                    if i in badlist: 
                        chain.detach_child(residue.id)
                    else:
                        NR += 1
        self.xyzO  = np.delete(xyzO,badlist,0)
        self.xyzH1 = np.delete(xyzH1,badlist,0)
        self.xyzH2 = np.delete(xyzH2,badlist,0)
        self.structure = structure
        self.NR = NR

    def keep(self,keeplist):
        alllist = range(self.NR)
        otherlist = np.squeeze(np.argwhere(~np.in1d(alllist,keeplist)))
        self.cullout(otherlist)

    def rotateit(self,viscinaldir,vshift):
        xyzO_step1 = copy.deepcopy(self.xyzO)
        xyzH1_step1 = copy.deepcopy(self.xyzH1)
        xyzH2_step1 = copy.deepcopy(self.xyzH2)
        xyzO_step2 = np.zeros(np.shape(self.xyzO))
        xyzH1_step2 = np.zeros(np.shape(self.xyzH1))
        xyzH2_step2 = np.zeros(np.shape(self.xyzH2))
        xbox = self.xbox
        ybox = self.ybox
        zbox = self.zbox
        nR,dum = xyzO_step1.shape

        if viscinaldir == 'y':  
    
            phi = np.arctan(vshift/zbox); print "phi =", phi*180/np.pi
            zboxp = zbox/np.cos(phi)
            yboxp = ybox*np.cos(phi)
            xboxp = xbox
            Rmat = np.array([[1,0,0],[0,np.cos(phi),np.sin(phi)],[0,-np.sin(phi),np.cos(phi)]])
            line1_m = -np.tan(phi)
            line1_b = ybox
            line2_m = line1_m
            line2_b = 0.0
        
            #print Rmat

            for i in range(nR):
   	        if (xyzO_step1[i,1] > line1_b + line1_m*xyzO_step1[i,2]):
	           xyzO_step1[i,1] = xyzO_step1[i,1] - ybox
	           xyzH1_step1[i,1] = xyzH1_step1[i,1] - ybox
	           xyzH2_step1[i,1] = xyzH2_step1[i,1] - ybox
	        if (xyzO_step1[i,1] < line2_b + line2_m*xyzO_step1[i,2]):
	           xyzO_step1[i,1] = xyzO_step1[i,1] + ybox
	           xyzH1_step1[i,1] = xyzH1_step1[i,1] + ybox
	           xyzH2_step1[i,1] = xyzH2_step1[i,1] + ybox
	        xyzO_step2[i,:] = np.dot(Rmat,xyzO_step1[i,:])
	        xyzH1_step2[i,:] = np.dot(Rmat,xyzH1_step1[i,:])
	        xyzH2_step2[i,:] = np.dot(Rmat,xyzH2_step1[i,:])
	        if (xyzO_step2[i,2]<0):
	           xyzO_step2[i,2] = xyzO_step2[i,2] + zboxp
	           xyzH1_step2[i,2] = xyzH1_step2[i,2] + zboxp
	           xyzH2_step2[i,2] = xyzH2_step2[i,2] + zboxp
	        if (xyzO_step2[i,2]<0):
	           xyzO_step2[i,2] = xyzO_step2[i,2] + zboxp
	           xyzH1_step2[i,2] = xyzH1_step2[i,2] + zboxp
	           xyzH2_step2[i,2] = xyzH2_step2[i,2] + zboxp
	        if (xyzO_step2[i,2]>zboxp):
	           xyzO_step2[i,2] = xyzO_step2[i,2] - zboxp
	           xyzH1_step2[i,2] = xyzH1_step2[i,2] - zboxp
	           xyzH2_step2[i,2] = xyzH2_step2[i,2] - zboxp
	        if (xyzO_step2[i,2]>zboxp):
	           xyzO_step2[i,2] = xyzO_step2[i,2] - zboxp
	           xyzH1_step2[i,2] = xyzH1_step2[i,2] - zboxp
	           xyzH2_step2[i,2] = xyzH2_step2[i,2] - zboxp

        else:
            print "Not implemented yet"
 
        #return xyzO_step2, xyzH1_step2, xyzH2_step2, xboxp, yboxp, zboxp
        self.xyzO = xyzO_step2
        self.xyzH1 = xyzH1_step2
        self.xyzH2 = xyzH2_step2
        self.xbox = xboxp
        self.ybox = yboxp
        self.zbox = zboxp
        return




def loadit(filename, nx, ny, nz, viscinaldir='y', nycel=1):

    # These are cell dimensions
    xcel = 4.4907312
    ycel = 7.7781746
    zcel = 3.6666666

    # Read in the pdb structure & specify the box size
    parser = Bio.PDB.PDBParser()
    #pdb.set_trace()
    #nx = 4; ny = 4, nz = 2
    structure = parser.get_structure('pdb', filename); xbox = xcel*nx; ybox = ycel*ny; zbox = zcel*nz

    # Set the shift information according to the viscinal surface we want
    if viscinaldir == 'y':
        vshift = ycel*nycel # This is the viscinal shift
        shift = np.array([\
            [ xbox,       0,        0      ], \
            [  0,        ybox+10,   0,     ], \
            [  0,        vshift,     zbox    ]])
    else:
        print "Not implemented yet"
        #break
                
    # A clumsy way to count the number of residues
    nR = 0
    for model in structure:
        for chain in model:
            j=0
            for residue in chain:
                for atom in residue:
                    j = j+1
    nR = j/3
    print "nR = ", nR
    
    # Allocate space for the various arrays
    xyz=np.zeros((nR*3,3))
    xyzO=np.zeros((nR,3))
    xyzH1=np.zeros((nR,3))
    xyzH2=np.zeros((nR,3))

    # Get all the coordinates
    for model in structure:
        for chain in model:
            j=0
            for residue in chain:
                for atom in residue:
                    temp1 = atom.get_coord()
                    #print temp1, j
                    xyz[j] = temp1
                    j = j+1

    # Sort the coordinates into O, H1, H2
    i=-1
    for j in range(nR*3):
        test=np.mod(j,3)
        if (test == 0):
            i = i+1
            xyzO[i]=xyz[j]
        elif (test ==1):
            xyzH1[i]=xyz[j]
        elif (test == 2):
            xyzH2[i]=xyz[j]
    
    #Return variables
    return xyzO,xyzH1,xyzH2, viscinaldir, shift, vshift, xbox, ybox, zbox, structure

def loaditnew(filename, xbox, ybox, zbox, viscinaldir='y'):

    # Read in the pdb structure & specify the box size
    parser = Bio.PDB.PDBParser()
    #pdb.set_trace()
    #nx = 4; ny = 4, nz = 2
    structure = parser.get_structure('pdb', filename)

    # Set the shift information according to the viscinal surface we want
    if viscinaldir == 'y':
        shift = np.array([\
            [ xbox,       0,        0      ], \
            [  0,        ybox+10,   0,     ], \
            [  0,         0,       zbox    ]])
    else:
        print "Not implemented yet"
        #break
                
    # A clumsy way to count the number of residues
    nR = 0
    for model in structure:
        for chain in model:
            j=0
            for residue in chain:
                for atom in residue:
                    j = j+1
    nR = j/3
    print "nR = ", nR
    
    # Allocate space for the various arrays
    xyz=np.zeros((nR*3,3))
    xyzO=np.zeros((nR,3))
    xyzH1=np.zeros((nR,3))
    xyzH2=np.zeros((nR,3))

    # Get all the coordinates
    for model in structure:
        for chain in model:
            j=0
            for residue in chain:
                for atom in residue:
                    temp1 = atom.get_coord()
                    #print temp1, j
                    xyz[j] = temp1
                    j = j+1

    # Sort the coordinates into O, H1, H2
    i=-1
    for j in range(nR*3):
        test=np.mod(j,3)
        if (test == 0):
            i = i+1
            xyzO[i]=xyz[j]
        elif (test ==1):
            xyzH1[i]=xyz[j]
        elif (test == 2):
            xyzH2[i]=xyz[j]
    
    #Return variables
    return xyzO,xyzH1,xyzH2,shift, structure


def saveit(filename,structure,xyzO,xyzH1,xyzH2):


    # Read in the pdb structure & specify the box size
    #parser = Bio.PDB.PDBParser()
    IO = Bio.PDB.PDBIO()
    IO.set_structure(structure)

    # Get all the coordinates
    i = -1 # residue counter
    for model in structure:
        for chain in model:
            j=0 # atom counter
            for residue in chain:
                for atom in residue:
                    test=np.mod(j,3)
                    if (test == 0):
                        i = i+1
                        temp = xyzO[i]
                    elif (test ==1):
                        temp = xyzH1[i]
                    elif (test == 2):
                	temp = xyzH2[i]
                    atom.set_coord(temp)
                    #print i,np.mod(j,3), temp
                    j = j+1
    IO.save(filename)
    
def findnbad(nni):
    nR, nk = nni.shape
    nbad = 0
    for i in range (nR):
        if np.size(np.where(nni[i]<0))>0:
            nbad = nbad+1
    #print "nbad = ", nbad
    return nbad

def finddefects(nni,nnltoi,nnitol):
    nR, nk = nni.shape; #print nR, nk
    Ddefect=np.zeros((nR*10,2)).astype(np.int32); nDdefect = 0
    Adefect=np.zeros((nR*10,2)).astype(np.int32); nAdefect = 0
    for i in range(nR):
        for k in range(4):
            if (nni[i,k]>=0):
                if (nnltoi[i,k]!=0) & (nnitol[i,k]!=0):
                    Ddefect[nDdefect,0]=i
                    Ddefect[nDdefect,1]=nni[i,k]
                    nDdefect += 1
                elif (nnltoi[i,k]==0) & (nnitol[i,k]==0):
                    Adefect[nAdefect,0]=i
                    Adefect[nAdefect,1]=nni[i,k]
                    nAdefect += 1
        #print i, nni[i], nnltoi[i], nnitol[i], Ddefect, Adefect
    Ddefect_ret = Ddefect[0:nDdefect]
    Adefect_ret = Adefect[0:nAdefect]
    return Ddefect_ret, Adefect_ret
    
def findthreezeros(nni,nnitol):
    nR, nk = nni.shape
    problem = False
    for i in range (nR):
        test = np.argwhere(nnitol[i]==0)
        if len(test)>2:
            #print "There's a problem at i, nni, nnitol =", i, nni[i], nnitol[i]
            problem = True
    return problem

def getnnitoletc(nni,xyzO,xyzH1,xyzH2,xyzshift):

    # Pre-allocate the nnitol and nnltoi arrays
    nR, nk = nni.shape; #print nR, nk
    nnitol=np.zeros((nR,4), dtype='int32')
    nnltoi=np.zeros((nR,4), dtype='int32')

    # Minimium projection required in order call this a donor
    HBproject = 2.7
    
    # Loop over all the residues
    for i in range(nR):
    
        #Is l(k) donating to i?
        for k in range(4):
            if nni[i,k]>=0:
                l= nni[i,k]
                viO=xyzO[i]
                vkO=xyzO[l]+xyzshift[i,k]
                vkH1=xyzH1[l]+xyzshift[i,k]
                vkH2=xyzH2[l]+xyzshift[i,k]
                vOO=viO-vkO
                vH1O=vkH1-vkO
                vH2O=vkH2-vkO
                #print i,l,np.dot(vOO, vH1O), np.dot(vOO, vH2O)
                if np.dot(vOO, vH1O)>HBproject:
                    nnltoi[i,k]=1
                if np.dot(vOO, vH2O)>HBproject:
                    nnltoi[i,k]=2
    
        #Is i donating to l(k)?
        for k in range(4):
	    if nni[i,k]>=0:
                l= nni[i,k]
                viO=xyzO[i]
                vkO=xyzO[l]+xyzshift[i,k]
                viH1=xyzH1[i]
                viH2=xyzH2[i]
                vOO=vkO-viO
            	vH1O=viH1-viO
                vH2O=viH2-viO
                #print i,l,np.dot(vOO, vH1O), np.dot(vOO, vH2O)
                if np.dot(vOO, vH1O)>HBproject:
                    nnitol[i,k]=1
                if np.dot(vOO, vH2O)>HBproject:
                    nnitol[i,k]=2
    
    return nnitol, nnltoi

            
def donorsandacceptors(nni,xyzO,xyzH1,xyzH2,xyzshift):

    # Pre-allocate the nnitol and nnltoi arrays
    nR, nk = nni.shape; #print nR, nk
    
    # Get the nnitol etc arrays
    nnitol, nnltoi = getnnitoletc(nni,xyzO,xyzH1,xyzH2,xyzshift)
    
    
    # Fixing the missing H2 cases
    for i in range (nR):
        
        # Find index of H2 for ith residue
        H2ofi=np.argwhere(nnitol[i]==2)
    
        # See if we didn't find H2
        if len(H2ofi) == 0:
        
            # OK, let's find a good place for it
            test = np.argwhere(nni[i]==-1)
            if len(test) > 0:
            
                # Find a spot that's not being used
                for j in range(len(test)): 
            
                    if nnitol[i,test[j]] == 0:
                        nnitol[i,test[j]]=2
                        #print nni[i], nnitol[i], nnitol[i]
                        break
            else:
                print "Internal inconsistency ..."
                f = np.sqrt(-1.0)

    # Fixing the missing H1 cases
    for i in range (nR):

        # Find index of H1 for ith residue
        H1ofi=np.argwhere(nnitol[i]==1)
    
        # See if we didn't find H2
        if len(H1ofi) == 0:
        
            # OK, let's find a good place for it
            test = np.argwhere(nni[i]==-1)
            if len(test) > 0:
            
                # Find a spot that's not being used
                for j in range(len(test)): 
            
                    if nnitol[i,test[j]] == 0:
                        nnitol[i,test[j]]=1
                        #print nni[i], nnitol[i], nnitol[i]
                        break
            else:
                print "Internal inconsistency ..."
                f = np.sqrt(-1.0)

    
    return nnitol, nnltoi
    
    
def fixsurface(nni,nnitol,nnltoi):

    # Check for defects
    Ddefect, Adefect = finddefects(nni,nnltoi,nnitol)

    #fix Ddefect with a "-1" nearest neighbor:
    for m in range(len(Ddefect)):
    
        # Pull out the index to the next residue that has a donor defect
        i = Ddefect[m,0]
        l = Ddefect[m,1]

        # See if this residue has any "-1" which means doesn't have a nearest neighbor    
        test = np.size(np.where(nni[i]<0))
    
        # If this defective residue is missing a nearest neighbor, point its hydrogen toward the space
        if (test>0):
            
            kzeroofi = np.squeeze(np.argwhere(nni[i]<0)[0])
            test2 = nnitol[i,kzeroofi]
            
            if (test2==0):
        
            
                kzeroofi = np.argwhere(nni[i]<0)[0]

            
                # Get the mutual pointer positions
                klofi = np.squeeze(np.argwhere(nni[i]==l))
                kiofl = np.squeeze(np.argwhere(nni[l]==i))

                # Fix it by changing nnitol
                #print 'fixing ', i, l, ' donor defect'
                #print 'before: ', i,l,nni[i],nni[l],nnitol[i],nnitol[l],nnltoi[l]
                temp = nnitol[i,klofi] # Save which of i's Hydrogens is the donor defective guy
                nnitol[i,klofi] = 0 # Point i's lone pair to l
                nnitol[i,kzeroofi] = temp # Point i's Hydrogen to what was -1 
                nnltoi[l,kiofl] = 0 # Confirms that i no longer points its H to l
                #print 'after: ', i,l,nni[i],nni[l],nnitol[i],nnitol[l],nnltoi[l]
        
    # Check for defects
    Ddefect, Adefect = finddefects(nni,nnltoi,nnitol)
                
    #Fix Adefect with a "-1" as a nearest neighbor(as long as nnltoi[i] does not have four 0's)
    #Same logic as fixing Ddefect with "-1" as nearest neighbor

    #print "Now we are fixing Adefects ..."
    for n in range(len(Adefect)):
        i = Adefect[n,0]
        l = Adefect[n,1]
        position = np.where(nni[i]<0)
        test = np.size(position)
        if test>0:
            klofi = np.squeeze(np.argwhere(nni[i]==l))
            kiofl = np.squeeze(np.argwhere(nni[l]==i))
        
            # Figure out which hydrogen is available (if any)
            timesH1isused = np.size(np.argwhere(nnitol[i]==1))
            timesH2isused = np.size(np.argwhere(nnitol[i]==2))
            whichone = 0
            if timesH1isused == 0:
                whichone = 1
            elif timesH2isused == 0:
                whichone = 2
            if whichone != 0:
                #print 'fixing ', i, l, ' acceptor defect'
                #print 'before: ', i,l,nni[i],nni[l],nnitol[i],nnitol[l],nnltoi[l]
                nnitol[i,klofi] = whichone
                nnltoi[l,kiofl] = whichone
                #print 'after:  ', i,l,nni[i],nni[l],nnitol[i],nnitol[l],nnltoi[l]
  

    return nnitol,nnltoi
    

def idsurfacedefects(nni,nnitol,nnltoi):
    nR, nk = nni.shape
    Ddefect, Adefect = finddefects(nni,nnitol,nnltoi)
    
    Ddefectsurf=np.zeros((nR*10,1)).astype(np.int32); nDdefectsurf=0
    Adefectsurf=np.zeros((nR*10,1)).astype(np.int32); nAdefectsurf=0
    
    for m in range(len(Ddefect)):
        
        i=Ddefect[m,0]
        test=np.size(np.where(nni[i]<0))
        
        if (test>0):
            Ddefectsurf[nDdefectsurf,0]=i
            nDdefectsurf +=1
            
    for n in range(len(Adefect)):
        
        i=Adefect[n,0]
        test=np.size(np.where(nni[i]<0))
        
        if (test>0):
            Adefectsurf[nAdefectsurf,0]=i
            nAdefectsurf +=1
            
    Ddefectsurf_ret = Ddefectsurf[0:nDdefectsurf]
    Adefectsurf_ret = Adefectsurf[0:nAdefectsurf]
    return Ddefectsurf_ret, Adefectsurf_ret


def getrefcoords(xyzO,xyzH1,xyzH2):
    lOH = np.sqrt(np.sum((xyzO-xyzH1)**2)); #print lOH
    lOH2 = np.sqrt(np.sum((xyzO-xyzH2)**2)); #print lOH2
    lHH = np.sqrt(np.sum((xyzH2-xyzH1)**2)); #print lHH
    ycoord = lHH/2
    theta = np.arccos(ycoord/lOH); #print ycoord/lOH
    zcoord = lOH*np.sin(theta)
    phi = 2*(np.pi/2-theta)
    vH1 = np.array([0.,-ycoord,-zcoord])
    vH2 = np.array([0.,ycoord,-zcoord])
    vO = np.array([0.,0.,0.])
    return vO, vH1, vH2
    
def getnni(xyzO,shift):

    # Calculate O-O distances
    nR,dum = xyzO.shape
    dist=np.zeros((nR,nR))
    for i in range(nR):
        xi=xyzO[i][0]; yi=xyzO[i][1]; zi=xyzO[i][2]
        for j in range(nR):
            xj=xyzO[j][0]; yj=xyzO[j][1]; zj=xyzO[j][2] 
            dist[i,j]=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
            #print i,j
        
    # Construct an initial nearest neighbor list, and arrays for related information 
    OOdist = 3.0
    nni=np.zeros((nR,4),dtype='int32')
    nnd=np.zeros((nR,4))
    xyzshift=np.zeros((nR,4,3))
    for i in range (nR):
        onecol = np.argsort(dist[i])
        x= dist[i,onecol[1:5]]
        nnd[i,]=x
        nni[i,]=onecol[1:5]
        for k in range(4):
            if nnd[i,k]>=OOdist:
                nni[i,k]=-1  # This is flag saying there is a missing nearest neighbor

    # Check out this list
    print "nbad before searching across periodic boundaries = ", findnbad(nni)

    # Find nearest neighbors across periodic boundaries & fix the nearest neighbor list accordingly
    for i in range (nR):
        if np.size(np.where(nni[i]<0))>0:
            xi=xyzO[i][0]; yi=xyzO[i][1]; zi=xyzO[i][2]
            for k in range (4):
                if nni[i,k]<0:
                    for j in range (nR):
                    
                        # Look across x
                        xj=xyzO[j][0]-shift[0][0]; yj=xyzO[j][1]-shift[0][1]; zj=xyzO[j][2]-shift[0][2]; disttest1=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                        xj=xyzO[j][0]+shift[0][0]; yj=xyzO[j][1]+shift[0][1]; zj=xyzO[j][2]+shift[0][2]; disttest2=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5

                        # Look across y
                        xj=xyzO[j][0]-shift[1][0]; yj=xyzO[j][1]-shift[1][1]; zj=xyzO[j][2]-shift[1][2]; disttest3=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                        xj=xyzO[j][0]+shift[1][0]; yj=xyzO[j][1]+shift[1][1]; zj=xyzO[j][2]+shift[1][2]; disttest4=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                    
                        # Look across z
                        xj=xyzO[j][0]-shift[2][0]; yj=xyzO[j][1]-shift[2][1]; zj=xyzO[j][2]-shift[2][2]; disttest5=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                        xj=xyzO[j][0]+shift[2][0]; yj=xyzO[j][1]+shift[2][1]; zj=xyzO[j][2]+shift[2][2]; disttest6=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                    
                        # Look across x and z
                        xj=xyzO[j][0]-shift[0][0]-shift[2][0]; 
                        yj=xyzO[j][1]-shift[0][1]-shift[2][1]; 
                        zj=xyzO[j][2]-shift[0][2]-shift[2][2]; 
                        disttest7=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                        
                        xj=xyzO[j][0]+shift[0][0]-shift[2][0]; 
                        yj=xyzO[j][1]+shift[0][1]-shift[2][1]; 
                        zj=xyzO[j][2]+shift[0][2]-shift[2][2];
                        disttest8=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                        
                        xj=xyzO[j][0]-shift[0][0]+shift[2][0]; 
                        yj=xyzO[j][1]-shift[0][1]+shift[2][1]; 
                        zj=xyzO[j][2]-shift[0][2]+shift[2][2];
                        disttest9=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                        
                        xj=xyzO[j][0]+shift[0][0]+shift[2][0]; 
                        yj=xyzO[j][1]+shift[0][1]+shift[2][1]; 
                        zj=xyzO[j][2]+shift[0][2]+shift[2][2];
                        disttest10=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                        
                        # Identify which cross-boundary search resulted in a nearest neighbor catch
                        T1 = disttest1<OOdist
                        T2 = disttest2<OOdist
                        T3 = disttest3<OOdist
                        T4 = disttest4<OOdist
                        T5 = disttest5<OOdist
                        T6 = disttest6<OOdist
                        T7 = disttest7<OOdist
                        T8 = disttest8<OOdist
                        T9 = disttest9<OOdist
                        T10 = disttest10<OOdist                   
                        whichone = np.where([T1,T2,T3,T4,T5,T6,T7,T8,T9,T10])
                        test = np.size(whichone)
                    
                        # Record the shift information; this logic must be consistent with other shift logic in this loop
                        if (test==0):
                            xyzshift[i,k] = np.array([0.,0.,0.])
                        elif T1:
                            xyzshift[i,k] = -shift[0] 
                        elif T2:
                            xyzshift[i,k] = shift[0]
                        elif T3:
                            xyzshift[i,k] = -shift[1]
                        elif T4:
                            xyzshift[i,k] = shift[1]
                        elif T5:
                            xyzshift[i,k] = -shift[2]
                        elif T6:
                            xyzshift[i,k] = shift[2]
                        elif T7:
                            xyzshift[i,k] = -shift[0]-shift[2]
                        elif T8:
                            xyzshift[i,k] = shift[0]-shift[2]
                        elif T9:
                            xyzshift[i,k] = -shift[0]+shift[2]
                        elif T10:
                            xyzshift[i,k] = shift[0]+shift[2]
                        
                            
                        
                        
                        # Record nearest neighbor information in the nni matrix, taking care to eliminate redundancies
                        if test>0:
                            if test>1:
                                print "Caught a double boundary case; the algorithm may not be valid"
                            taken = 0
                            for kp in range(4):
                                jp = nni[i,kp]
                                if jp==j:
                                    #print "already taken ..."
                                    taken = 1
                            if taken==0:
                                temp = copy.deepcopy(nni[i])
                                nni[i,k]=j
                                #print i, temp, "was replaced by", nni[i], "because of ", xyzshift[i,k] 
                                break

    # Check out this list
    print "nbad after searching across periodic boundaries = ", findnbad(nni)

    return nni, xyzshift
    
    
def checkfordefects(nni,xyzO,xyzH1,xyzH2,xyzshift):
    # Get an initial list of donors and acceptors
    nnitol, nnltoi = donorsandacceptors(nni,xyzO,xyzH1,xyzH2,xyzshift)
    Ddefect, Adefect = finddefects(nni,nnltoi,nnitol)
    print "Found defects:", len(Ddefect), len(Adefect)
    return nnitol, nnltoi, Ddefect, Adefect 
    #Eliminate second half of Ddefect and Adefect 

def fixit(nni, nnitol_in, nnltoi_in, nprop):
    import random
    #import pdb
    #import copy

    # Fix any surface defects
    nnitol,nnltoi = fixsurface(nni,nnitol_in,nnltoi_in)
    Ddefect, Adefect = finddefects(nni,nnltoi,nnitol)
    print "After initial surface fix:", len(Ddefect), len(Adefect)

    # Propagate a defect and check
    icount = 0
    for iprop in range(nprop):
        if len(Ddefect) > 0:
            icount += 1
            #print "****"
            #print "Starting iprop = ", iprop, " w/NDdefect, NAdefect = ", len(Ddefect), len(Adefect)
            #m = random.randint(0,len(Ddefect)-1)
            m = 0
            i = Ddefect[m,0]
            l = Ddefect[m,1]
            #print "working on ", i, l
            klofi = np.squeeze(np.argwhere(nni[i]==l))
            kiofl = np.squeeze(np.argwhere(nni[l]==i))
    
            # Figure out which Hydrogen of i (1 or 2) is pointing to l
            Hofi = nnitol[i,klofi]; #print Hofi    

            # Decide on a new nearest neighbor to point this Hydrogen to
            first=np.squeeze(np.argwhere(nnitol[i]==0)[0]) 
            second=np.squeeze(np.argwhere(nnitol[i]==0)[1])
            zeroorone = random.randint(0,1) #sloppy way to get random index
            if zeroorone == 0:
                klofip = first
            else:
                klofip = second
            #print "first, second, klofip = ", first, second, klofip
            lp = nni[i,klofip]; #print lp
            kioflp = np.squeeze(np.argwhere(nni[lp]==i)); #print kioflp

            # Point a lone pair to l instead
            nnitol[i,klofi] = 0
            nnltoi[l,kiofl] = 0
    
            # Point the offending hydrogen of i to next victim (nearest neighbor)
            nnitol[i,klofip] = Hofi
            nnltoi[lp,kioflp] = Hofi
    
            # Check for surface defects
            nnitol,nnltoi = fixsurface(nni,nnitol,nnltoi)
            Ddefect, Adefect = finddefects(nni,nnltoi,nnitol)
    print "After fixing Ddefects:", len(Ddefect), len(Adefect), "which took", icount, "iterations"

    #Propagate Adefects
    icount = 0
    for iprop in range(nprop):
        if len(Adefect) > 0:
            icount += 1
            #print "****"
        
            #print "Starting iprop = ", iprop, " w/NDdefect, NAdefect = ", len(Ddefect), len(Adefect)
            m = random.randint(0,len(Adefect)-1)
            i = Adefect[m,0]
            l = Adefect[m,1]
            #print "working on ", i, l
            klofi = np.squeeze(np.argwhere(nni[i]==l))
            kiofl = np.squeeze(np.argwhere(nni[l]==i))  

            # Decide on a new nearest neighbor to point this lone pair to
            knonzeros=np.squeeze(np.argwhere(nnitol[i]!=0))
            ikrandom = random.randint(0,len(knonzeros)-1) #sloppy way to get random index
            klofip = knonzeros[ikrandom]
            #print "nnitol[i], klofip = ", nnitol[i], klofip
            lp = nni[i,klofip]; #print lp
            kioflp = np.squeeze(np.argwhere(nni[lp]==i)); #print kioflp

            # Point a hydrogen to l instead
            Hofi=nnitol[i,klofip]
            nnitol[i,klofi] = Hofi
            nnltoi[l,kiofl] = Hofi
    
            # Point the offending lone pair of i to next victim (nearest neighbor)
            nnitol[i,klofip] = 0
            nnltoi[lp,kioflp] = 0
        
            #pdb.set_trace()
            nnitol,nnltoi = fixsurface(nni,nnitol,nnltoi)
            Ddefect, Adefect = finddefects(nni,nnltoi,nnitol)
    print "After fixing Adefects:", len(Ddefect), len(Adefect), "which took", icount, "iterations"
    return nnitol, nnltoi


'''#FIXIT REDONE
def fixit(nni, nnitol_in, nnltoi_in, nprop):
    import random
    #import pdb
    #import copy

    # Fix any surface defects
    nnitol,nnltoi = fixsurface(nni,nnitol_in,nnltoi_in)
    Ddefect, Adefect = finddefects(nni,nnltoi,nnitol)
    print "After initial surface fix:", len(Ddefect), len(Adefect)

    # Propagate a defect and check
    icount = 0
    for iprop in range(nprop):
        if len(Ddefect) > 0:
            icount += 1
            #print "****"
            #print "Starting iprop = ", iprop, " w/NDdefect, NAdefect = ", len(Ddefect), len(Adefect)
            #m = random.randint(0,(len(Ddefect)/2)-1)
            i = Ddefect[0,0]
            l = Ddefect[0,1]
            #print "working on ", i, l
            klofi = np.squeeze(np.argwhere(nni[i]==l))
            kiofl = np.squeeze(np.argwhere(nni[l]==i))
    
            # Figure out which Hydrogen of i (1 or 2) is pointing to l
            Hofi = nnitol[i,klofi]; #print Hofi    

            # Decide on a new nearest neighbor to point this Hydrogen to
            first=np.squeeze(np.argwhere(nnitol[i]==0)[0]) 
            second=np.squeeze(np.argwhere(nnitol[i]==0)[1])
            zeroorone = random.randint(0,1) #sloppy way to get random index
            if zeroorone == 0:
                klofip = first
            else:
                klofip = second
            #print "first, second, klofip = ", first, second, klofip
            lp = nni[i,klofip]; #print lp
            kioflp = np.squeeze(np.argwhere(nni[lp]==i)); #print kioflp

            # Point a lone pair to l instead
            nnitol[i,klofi] = 0
            nnltoi[l,kiofl] = 0
    
            # Point the offending hydrogen of i to next victim (nearest neighbor)
            nnitol[i,klofip] = Hofi
            nnltoi[lp,kioflp] = Hofi
            
            
            nR, nk = nni.shape
            for i in range(nR):
                if nnitol[i,klofip]>=1 and nnitol[lp,kioflp]>=1:
                
                    #change the the value of i and l to trace the defect
                    temp = i
                    i = lp
                    l = temp
                    
                    #See if new i is at the surface
                    test = np.size(np.where(nni[i]<0))
                
    
                    # If this defective residue is missing a nearest neighbor, point its hydrogen toward the space
                    if (test>0):
            
                        kzeroofi = np.squeeze(np.argwhere(nni[i]<0)[0])
                        test2 = nnitol[i,kzeroofi]
                
                        if (test2==0):
        
            
                            kzeroofi = np.argwhere(nni[i]<0)[0]

            
                            # Get the mutual pointer positions
                            klofi = np.squeeze(np.argwhere(nni[i]==l))
                            kiofl = np.squeeze(np.argwhere(nni[l]==i))
    
                            # Fix it by changing nnitol
                            #print 'fixing ', i, l, ' donor defect'
                            #print 'before: ', i,l,nni[i],nni[l],nnitol[i],nnitol[l],nnltoi[l]
                            temp = nnitol[i,klofi] # Save which of i's Hydrogens is the donor defective guy
                            nnitol[i,klofi] = 0 # Point i's lone pair to l
                            nnitol[i,kzeroofi] = temp # Point i's Hydrogen to what was -1 
                            nnltoi[l,kiofl] = 0 # Confirms that i no longer points its H to l
            
                    #if i is not at the surface then create a new defect to fix existing one
                    else: 
                        klofi = np.squeeze(np.argwhere(nni[i]==l))
                        kiofl = np.squeeze(np.argwhere(nni[l]==i))
    
                        # Figure out which Hydrogen of i (1 or 2) is pointing to l
                        Hofi = nnitol[i,klofi]; #print Hofi    

                        # Decide on a new nearest neighbor to point this Hydrogen to
                        first=np.squeeze(np.argwhere(nnitol[i]==0)[0]) 
                        second=np.squeeze(np.argwhere(nnitol[i]==0)[1])
                        zeroorone = random.randint(0,1) #sloppy way to get random index
                        if zeroorone == 0:
                            klofip = first
                        else:
                            klofip = second
                        #print "first, second, klofip = ", first, second, klofip
                        lp = nni[i,klofip]; #print lp
                        kioflp = np.squeeze(np.argwhere(nni[lp]==i)); #print kioflp

                        # Point a lone pair to l instead
                        nnitol[i,klofi] = 0
                        nnltoi[l,kiofl] = 0
    
                        # Point the offending hydrogen of i to next victim (nearest neighbor)
                        nnitol[i,klofip] = Hofi
                        nnltoi[lp,kioflp] = Hofi

            # Check for surface defects
            nnitol,nnltoi = fixsurface(nni,nnitol,nnltoi)
            Ddefect, Adefect = finddefects(nni,nnltoi,nnitol)
    print "After fixing Ddefects:", len(Ddefect), len(Adefect), "which took", icount, "iterations"

    #Propagate Adefects
    icount = 0
    for iprop in range(nprop):
        if len(Adefect) > 0:
            icount += 1
            #print "****"
        
            #print "Starting iprop = ", iprop, " w/NDdefect, NAdefect = ", len(Ddefect), len(Adefect)
            #m = random.randint(0,(len(Adefect)/2)-1)
            i = Adefect[0,0]
            l = Adefect[0,1]
            #print "working on ", i, l
            klofi = np.squeeze(np.argwhere(nni[i]==l))
            kiofl = np.squeeze(np.argwhere(nni[l]==i))  

            # Decide on a new nearest neighbor to point this lone pair to
            knonzeros=np.squeeze(np.argwhere(nnitol[i]!=0))
            ikrandom = random.randint(0,len(knonzeros)-1) #sloppy way to get random index
            klofip = knonzeros[ikrandom]
            #print "nnitol[i], klofip = ", nnitol[i], klofip
            lp = nni[i,klofip]; #print lp
            kioflp = np.squeeze(np.argwhere(nni[lp]==i)); #print kioflp

            # Point a hydrogen to l instead
            Hofi=nnitol[i,klofip]
            nnitol[i,klofi] = Hofi
            nnltoi[l,kiofl] = Hofi
    
            # Point the offending lone pair of i to next victim (nearest neighbor)
            nnitol[i,klofip] = 0
            nnltoi[lp,kioflp] = 0
            
            
            #trace defect through crystal
            nR, nk = nni.shape
            for i in range(nR):
                if nnitol[i,klofip]==0 and nnitol[lp,kioflp]==0:
                
                    #change the the value of i and l to trace the defect
                    temp = i
                    i = lp
                    l = temp
                
                    test = np.size(np.where(nni[i]<0))
                
                    #check if surface residue, if so then fix it
                    if test>0:
                        klofi = np.squeeze(np.argwhere(nni[i])==l)
                        kiofl=np.squeeze(np.argwhere(nni[l]==i))
                    
                        # Figure out which hydrogen is available (if any)
                        timesH1isused = np.size(np.argwhere(nnitol[i]==1))
                        timesH2isused = np.size(np.argwhere(nnitol[i]==2))
                        whichone = 0
                        if timesH1isused == 0:
                            whichone = 1
                        elif timesH2isused == 0:
                            whichone = 2
                        if whichone != 0:
                            nnitol[i,klofi] = whichone
                            nnltoi[l,kiofl] = whichone
                        
                    #if i is not at surface - create new defect and check/trace
                    else:
                        klofi = np.squeeze(np.argwhere(nni[i]==l))
                        kiofl = np.squeeze(np.argwhere(nni[l]==i))  
    
                        # Decide on a new nearest neighbor to point this lone pair to
                        knonzeros=np.squeeze(np.argwhere(nnitol[i]!=0))
                        ikrandom = random.randint(0,len(knonzeros)-1) #sloppy way to get random index
                        klofip = knonzeros[ikrandom]
                        #print "nnitol[i], klofip = ", nnitol[i], klofip
                        lp = nni[i,klofip]; #print lp
                        kioflp = np.squeeze(np.argwhere(nni[lp]==i)); #print kioflp

                        # Point a hydrogen to l instead
                        Hofi=nnitol[i,klofip]
                        nnitol[i,klofi] = Hofi
                        nnltoi[l,kiofl] = Hofi
            
                        # Point the offending lone pair of i to next victim (nearest neighbor)
                        nnitol[i,klofip] = 0
                        nnltoi[lp,kioflp] = 0  
                    
            nnitol,nnltoi = fixsurface(nni,nnitol,nnltoi)
            Ddefect, Adefect = finddefects(nni,nnltoi,nnitol)
    print "After fixing Adefects:", len(Ddefect), len(Adefect), "which took", icount, "iterations"
    return nnitol, nnltoi'''


def reconstructit(xyzO, xyzH1, xyzH2, nni, nnitol, xyzshift):
    # Orient residues in a desired config

    # Get water residue in the reference configuration;
    # This is oxygen at the origin, H1 and H2 in the y-z plane, both with z<0, but H1 has y<0 and H2 has y>0 
    i=0; vO, vH1, vH2 = getrefcoords(xyzO[i],xyzH1[i],xyzH2[i])

    # Create new coordinate arrays
    xyzH1_new = copy.deepcopy(xyzH1); #np.zeros(np.shape(xyzO))
    xyzH2_new = copy.deepcopy(xyzH2); #np.zeros(np.shape(xyzO))
    xyzO_new = copy.deepcopy(xyzO) # For convenience; this will be two handles on the same data

    #Reference residue vectors ... not needed because in this representation, Vnot=identity matrix
    #v21not=vH2-vH1; v21not = v21not/np.linalg.norm(v21not)
    #vBnot=(2*vO)-(vH1+vH2); vBnot = vBnot/np.linalg.norm(vBnot)
    #vNnot=np.cross(v21not,vBnot); vNnot = vNnot/np.linalg.norm(vNnot)
    #Vnot = np.transpose([vNnot,v21not,vBnot])

    nR, dum = nni.shape
    for i in range (nR):
    
        #Find index of H1 and H2 for ith residue
        H1ofi=np.squeeze(np.argwhere(nnitol[i]==1))
        H2ofi=np.squeeze(np.argwhere(nnitol[i]==2))
        L1ofi=np.squeeze(np.argwhere(nnitol[i]==0))[0]
        L2ofi=np.squeeze(np.argwhere(nnitol[i]==0))[1]
    
        #Find corresponding residues that H1 and H2 point at
        k1ofi=nni[i,H1ofi]
        k2ofi=nni[i,H2ofi]
        j1ofi=nni[i,L1ofi]
        j2ofi=nni[i,L2ofi]
    
        # if two hydrogens have nearest neighbors, use them directly
        if (k2ofi>=0) & (k1ofi>=0):

            # Get the tensor using hydrogen nearest neighbors
            v21=(xyzO[k2ofi]+xyzshift[i,H2ofi]) -(xyzO[k1ofi]+xyzshift[i,H1ofi]);            v21 = v21/np.linalg.norm(v21)
            vB=(2*xyzO[i])-(xyzO[k1ofi]+xyzshift[i,H1ofi]+xyzO[k2ofi]+xyzshift[i,H2ofi]);    vB  = vB/np.linalg.norm(vB)
            vN=np.cross(v21,vB);                                                             vN = vN/np.linalg.norm(vN)

            # Construct the rotation matrix        
            R = np.transpose([vN,v21,vB])
        
            # Apply the rotation matrix to the hydrogens
            xyzH1_new[i] = np.dot(R,vH1)+xyzO[i]
            xyzH2_new[i] = np.dot(R,vH2)+xyzO[i]
            #vOdum, vH1dum, vH2dum = vs.getrefcoords(xyzO_new[i],xyzH1_new[i],xyzH2_new[i])

        # if two lone pairs have nearest neighbors, use them (with a tensor correction)
        elif (j2ofi>=0) & (j1ofi>=0):

            # Get the tensor using lone pair nearest neighbors
            v21=(xyzO[j2ofi]+xyzshift[i,L2ofi]) - (xyzO[j1ofi]+xyzshift[i,L1ofi]);           v21 = v21/np.linalg.norm(v21)
            vB=(2*xyzO[i])-(xyzO[j1ofi]+xyzshift[i,L1ofi]+xyzO[j2ofi]+xyzshift[i,L2ofi]);    vB  = vB/np.linalg.norm(vB)
            vN=np.cross(v21,vB);                                                             vN = vN/np.linalg.norm(vN)

            # Use lone pair tensor to construct the hydrogen tensor         
            v21_p = vN
            vB_p = -vB
            vN_p = np.cross(v21_p,vB_p)

            # Construct the rotation matrix                
            R = np.transpose([vN_p,v21_p,vB_p])
        
            # Apply the rotation matrix to the hydrogens
            xyzH1_new[i] = np.dot(R,vH1)+xyzO[i]
            xyzH2_new[i] = np.dot(R,vH2)+xyzO[i]
            #vOdum, vH1dum, vH2dum = vs.getrefcoords(xyzO_new[i],xyzH1_new[i],xyzH2_new[i])

        # if we got here, there must be one of each, and two free slots ...
        # One possibility is H1 and a lone pair
        elif (j1ofi>=0) & (k1ofi>=0):
        
            # Just checking ... there had better be two free slots
            test = np.squeeze(np.argwhere(nni[i]==-1))
            if len(test)!=2:
                print "Unexpected combination in Reconstructit"

            # Get the tensor using hydrogen & lone pair nearest neighbors ... H is "H1" and the 1st lone pair is "H2"
            v21=(xyzO[j1ofi]+xyzshift[i,L1ofi])-(xyzO[k1ofi]+xyzshift[i,H1ofi]);                v21 = v21/np.linalg.norm(v21)
            vB=(2*xyzO[i])-(xyzO[k1ofi]+xyzshift[i,H1ofi]+xyzO[j1ofi]+xyzshift[i,L1ofi]);       vB  = vB/np.linalg.norm(vB)
            vN=np.cross(v21,vB);                                                                vN = vN/np.linalg.norm(vN)

            # Construct the rotation matrix        
            R = np.transpose([vN,v21,vB])

            # Apply the rotation matrix to H1
            xyzH1_new[i] = np.dot(R,vH1)+xyzO[i]

            # Get the other hydrogen & lone pair tensor to construct the other tensor ... H is "H2" and the 1st lone pair is "H1"  
            v21_p = vN
            vB_p = -vB
            vN_p = np.cross(v21_p,vB_p)
           
            # Construct the rotation matrix        
            R = np.transpose([vN_p,v21_p,vB_p])
                                             
            # Apply the rotation matrix to H2
            xyzH2_new[i] = np.dot(R,vH2)+xyzO[i]
            #vOdum, vH1dum, vH2dum = vs.getrefcoords(xyzO_new[i],xyzH1_new[i],xyzH2_new[i])
                                         
        # if we got here, there must be one of each, and two free slots ...
        # The other possibility is H2 and a lone pair
        elif (j1ofi>=0) & (k2ofi>=0):
        
            # Just checking ... there had better be two free slots
            test = np.squeeze(np.argwhere(nni[i]==-1))
            if len(test)!=2:
                print "Unexpected combination in Reconstructit"

            # Get the tensor using hydrogen & lone pair nearest neighbors ... H is "H2" and the 1st lone pair is "H1"
            v21=(xyzO[k2ofi]+xyzshift[i,H2ofi])-(xyzO[j1ofi]+xyzshift[i,L1ofi]);                v21 = v21/np.linalg.norm(v21)
            vB=(2*xyzO[i])-(xyzO[k2ofi]+xyzshift[i,H2ofi]+xyzO[j1ofi]+xyzshift[i,L1ofi]);       vB  = vB/np.linalg.norm(vB)
            vN=np.cross(v21,vB);                                                                vN = vN/np.linalg.norm(vN)

            # Construct the rotation matrix        
            R = np.transpose([vN,v21,vB])
        
            # Apply the rotation matrix to H2
            xyzH2_new[i] = np.dot(R,vH2)+xyzO[i]

            # Get the other hydrogen & lone pair tensor to construct the other tensor ... H is "H2" and the 1st lone pair is "H1"  
            v21_p = vN
            vB_p = -vB
            vN_p = np.cross(v21_p,vB_p)
           
            # Construct the rotation matrix        
            R = np.transpose([vN_p,v21_p,vB_p])
                                             
            # Apply the rotation matrix to H1
            xyzH1_new[i] = np.dot(R,vH1)+xyzO[i]
            vOdum, vH1dum, vH2dum = getrefcoords(xyzO_new[i],xyzH1_new[i],xyzH2_new[i])
                                         
        else:
            print "Unexpected combination in Reconstructit", i, nni[i], nnitol[i]
    
    return xyzO_new, xyzH1_new, xyzH2_new
    
def rotateit(xyzO_new, xyzH1_new, xyzH2_new, viscinaldir, shift, vshift, xbox, ybox, zbox):
    xyzO_step1 = copy.deepcopy(xyzO_new)
    xyzH1_step1 = copy.deepcopy(xyzH1_new)
    xyzH2_step1 = copy.deepcopy(xyzH2_new)
    xyzO_step2 = np.zeros(np.shape(xyzO_new))
    xyzH1_step2 = np.zeros(np.shape(xyzO_new))
    xyzH2_step2 = np.zeros(np.shape(xyzO_new))

    nR,dum = xyzO_new.shape

    if viscinaldir == 'y':  
    
        phi = np.arctan(vshift/zbox); print "phi =", phi*180/np.pi
        zboxp = zbox/np.cos(phi)
        yboxp = ybox*np.cos(phi)
        xboxp = xbox
        Rmat = np.array([[1,0,0],[0,np.cos(phi),np.sin(phi)],[0,-np.sin(phi),np.cos(phi)]])
        line1_m = -np.tan(phi)
        line1_b = ybox
        line2_m = line1_m
        line2_b = 0.0
        
        #print Rmat

        for i in range(nR):
   	   if (xyzO_step1[i,1] > line1_b + line1_m*xyzO_step1[i,2]):
	       xyzO_step1[i,1] = xyzO_step1[i,1] - ybox
	       xyzH1_step1[i,1] = xyzH1_step1[i,1] - ybox
	       xyzH2_step1[i,1] = xyzH2_step1[i,1] - ybox
	   if (xyzO_step1[i,1] < line2_b + line2_m*xyzO_step1[i,2]):
	       xyzO_step1[i,1] = xyzO_step1[i,1] + ybox
	       xyzH1_step1[i,1] = xyzH1_step1[i,1] + ybox
	       xyzH2_step1[i,1] = xyzH2_step1[i,1] + ybox
	   xyzO_step2[i,:] = np.dot(Rmat,xyzO_step1[i,:])
	   xyzH1_step2[i,:] = np.dot(Rmat,xyzH1_step1[i,:])
	   xyzH2_step2[i,:] = np.dot(Rmat,xyzH2_step1[i,:])
	   if (xyzO_step2[i,2]<0):
	       xyzO_step2[i,2] = xyzO_step2[i,2] + zboxp
	       xyzH1_step2[i,2] = xyzH1_step2[i,2] + zboxp
	       xyzH2_step2[i,2] = xyzH2_step2[i,2] + zboxp
	   if (xyzO_step2[i,2]<0):
	       xyzO_step2[i,2] = xyzO_step2[i,2] + zboxp
	       xyzH1_step2[i,2] = xyzH1_step2[i,2] + zboxp
	       xyzH2_step2[i,2] = xyzH2_step2[i,2] + zboxp
	   if (xyzO_step2[i,2]>zboxp):
	       xyzO_step2[i,2] = xyzO_step2[i,2] - zboxp
	       xyzH1_step2[i,2] = xyzH1_step2[i,2] - zboxp
	       xyzH2_step2[i,2] = xyzH2_step2[i,2] - zboxp
	   if (xyzO_step2[i,2]>zboxp):
	       xyzO_step2[i,2] = xyzO_step2[i,2] - zboxp
	       xyzH1_step2[i,2] = xyzH1_step2[i,2] - zboxp
	       xyzH2_step2[i,2] = xyzH2_step2[i,2] - zboxp
	    


    else:
        print "Not implemented yet"
 
    return xyzO_step2, xyzH1_step2, xyzH2_step2, xboxp, yboxp, zboxp

def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)