import copy
xyzO_step1 = copy.deepcopy(xyzO_new)
xyzH1_step1 = copy.deepcopy(xyzH1_new)
xyzH2_step1 = copy.deepcopy(xyzH2_new)
xyzO_step2 = np.zeros(np.shape(xyzO))
xyzH1_step2 = np.zeros(np.shape(xyzO))
xyzH2_step2 = np.zeros(np.shape(xyzO))

if viscinaldir == 'y':  
    
    phi = np.arctan(yshift/zbox)    
    zboxp = zbox/np.cos(phi)
    yboxp = ybox*np.cos(phi)
    Rmat = np.array([[1,0,0],[0,np.cos(phi),np.sin(phi)],[0,-np.sin(phi),np.cos(phi)]])
    line1_m = -np.tan(phi)
    line1_b = ybox
    line2_m = line1_m
    line2_b = 0.0

      
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
	if (xyzO_step2[i,2]>zboxp):
	    xyzO_step2[i,2] = xyzO_step2[i,2] - zboxp
	    xyzH1_step2[i,2] = xyzH1_step2[i,2] - zboxp
	    xyzH2_step2[i,2] = xyzH2_step2[i,2] - zboxp
	    


else:
    print "Not implemented yet"
    #break

