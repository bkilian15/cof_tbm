from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.visualize import view
from ase.io import read, write
import os
import numpy as np
from math import *


def get_total_distance(pts, numk):
    
    ds = []
    for i, pt in enumerate(pts):
        p1 = pt
        p2 = pts[i+1]
        d = sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)
        ds.append(d)
        
        if i == len(pts) - 2:
            break
    td = sum(ds)
    
    nks = []
    for prop in ds:
        nks.append(round((prop/td)*numk))
        

    return nks

def create_path(fullpath, pts, numk=100):
    """
    takes the points (including name and coords) and the number of desired total k-points
    and returns all k-points along the path
    """
    pt_coords = []
    for pt in pts:
        pt_coords.append(fullpath[pt])
        
    num_steps = get_total_distance(pt_coords, numk=numk)
    #num_steps = total_distance/(len(pts) -1)
    print(num_steps, sum(num_steps))
    
    
    kpath = []
    for i, pt in enumerate(pts):
        p1 = pt_coords[i]
        p2 = pt_coords[i+1]
        a = np.linspace(p1,p2,num_steps[i])
        for p in a:
            kpath.append(p)
            
        if i == len(pts) - 2:
            break
        
    return kpath
