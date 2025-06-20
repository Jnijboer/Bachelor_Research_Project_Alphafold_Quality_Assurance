#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 19:02:44 2025

@author: mariabucheru
"""

import numpy as np
from pymol import cmd
import time
import matplotlib.pyplot as plt


def objlist():
    '''
     This function retrieves a list of objects from a set of structures loaded 
     into Pymol and performs several steps to manipulate the list using iterations. 
     It prints the entire list,individual elements, and pairs of elements.
     Additionally, it performs and times an alignment operation between the 
     first two objects in the list.
    '''
    objects = cmd.get_object_list()
    print(objects)
    
    for i in objects: print(i)
    for idx, obj1 in enumerate(objects): print(idx, obj1)
    for idx, obj1 in enumerate(objects):
        for obj2 in objects[idx+1:]: print(idx, obj1, obj2)
    print(obj1, obj2)
    print(cmd.align(objects[0], objects[1], cycles=0))
    print(-time.time() + (cmd.align(objects[0], objects[1], cycles=0) and time.time()))
    return objects


def align_scores(objects):
    '''
    This function computes pairwise alignments between each pair of objects in the given list.
    It iterates over each object and aligns it with every subsequent object in the list.
    The alignment results are stored in a list, and the time taken for each alignment is printed.
    '''
    alignments = []
    for idx, obj1 in enumerate(objects):
        for obj2 in objects[idx+1:]:
            print(obj1, obj2, -time.time() + (alignments.append(cmd.align(obj1, obj2)) or time.time()))
    return alignments


def trimatrix(objects, alignments):
    '''
    This function creates a square matrix initialized with zeros and populates the upper triangular part
    (excluding the diagonal) with alignment scores as RMSD extracted from the alignments list. The matrix is then
    made symmetric by adding its transpose. The resulting matrix is printed, displayed as a heatmap, and returned.
    '''
    A= np.zeros((len (objects), len(objects)))
    A[np.triu_indices(len(objects),1) ]= [aln[3] for aln in alignments]
    A += A.T
    print (A)
    plt.imshow(A)
    plt.show(A)
    return A


def clustering(A):
    '''
    This function computes the principal coordinates of a given matrix A by first centering it
    and then performing an eigen decomposition. The resulting principal coordinates are used to
    create scatter plots for visualization. The first plot is a 2D scatter plot of the first two
    principal coordinates, and the second plot is colored by the third principal coordinate.
    '''
    m = A.mean(axis=0)
    vals, vecs = np.linalg.eigh(-0.5 * (A - m[:,None] - m[None,:] + m.mean()))
    print(vals)
    pcoords = vecs[:, ::-1] * vals[::-1] **0.5
    
    plt.scatter(pcoords[:, 0], pcoords[:, 1])
    plt.gca().set_aspect('equal')
    plt.show()
    
    plt.scatter(pcoords[:, 0], pcoords[:, 1], c=pcoords[:, 2])
    plt.show()
    
    
def full_clustering_analysis():
    '''
    This function combines all the functions made before, so the analysis can be
    executed with one step to get the results.
    '''
    objects = objlist()
    alignments = align_scores(objects)
    A = trimatrix(objects, alignments)
    clustering(A)
