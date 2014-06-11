from mpl_toolkits.mplot3d import *
import matplotlib.pyplot as plt
import numpy as np
fig = plt.figure(1)
fig.clf()
ax = fig.gca(projection='3d')   
ax.scatter(xyzO[:,0],xyzO[:,1],xyzO[:,2]);                        # plot a 3d scatter plot
ax.set_xlabel('x label')
ax.set_ylabel('y label')
ax.set_zlabel('z label')
ax.set_xlim3d(np.min(xyzO[:,0])-5,np.max(xyzO[:,0]+5))
ax.set_ylim3d(np.min(xyzO[:,1])-5,np.max(xyzO[:,1]+5))
ax.set_zlim3d(np.min(xyzO[:,2])-5,np.max(xyzO[:,2]+5))
plt.title('Original')
plt.show()

fig = plt.figure(2)
fig.clf()
ax = fig.gca(projection='3d')   
ax.scatter(xyzO_rot[:,0],xyzO_rot[:,1],xyzO_rot[:,2]);                        # plot a 3d scatter plot
ax.set_xlabel('x label')
ax.set_ylabel('y label')
ax.set_zlabel('z label')
ax.set_xlim3d(np.min(xyzO_rot[:,0])-5,np.max(xyzO_rot[:,0]+5))
ax.set_ylim3d(np.min(xyzO_rot[:,1])-5,np.max(xyzO_rot[:,1]+5))
ax.set_zlim3d(np.min(xyzO_rot[:,2])-5,np.max(xyzO_rot[:,2]+5))
plt.title('Viscinal')
plt.show()