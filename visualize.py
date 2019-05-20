#%%
import numpy as np
import os
import json
import glob
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import sys


#%%
def bar_range(col, r):
    return np.arange(min(col), max(col) + r, r)

#%%
sys.path.append('spatialtree')
points = json.load(open('spatialtree/points_annotated.json'))['points']
points = np.array(points)

xdata = points[:,0]
ydata = points[:,2]
zdata = points[:,1]
cdata = points[:,3]
fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection='3d')
ax.scatter3D(xdata, ydata, zdata, c=cdata, s=1)

#%%
