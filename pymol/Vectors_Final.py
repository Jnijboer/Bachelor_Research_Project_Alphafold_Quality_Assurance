#!/usr/bin/env python
# coding: utf-8

# In[1]:


import trr
import MDAnalysis as mda
import numpy as np

u = mda.Universe(tpr, xtc)
ref = mda.Universe(pdb)
protein = u.select_atoms('protein')
T = trr.TRR(trr_file, selection=protein.ix)
F = T.forces
X = np.zeros((len(u.trajectory), len(protein), 3))
for ix, frame in enumerate(u.trajectory): X[ix] = protein.positions


# In[ ]:


def dim2pbc(arr: np.ndarray) -> np.ndarray:
    '''
    Convert unit cell definition from PDB CRYST1 format to lattice definition.
    '''
    lengths = arr[:, :3]
    angles = arr[:, 3:] * (np.pi / 180)

    cosa = np.cos(angles)
    sing = np.sin(angles[:, 2])    
    
    pbc = np.zeros((len(arr), 9))
    pbc[:, 0] = lengths[:, 0]
    pbc[:, 3] = lengths[:, 1] * cosa[:, 2]
    pbc[:, 4] = lengths[:, 1] * sing    
    pbc[:, 6] = lengths[:, 2] * cosa[:, 1]
    pbc[:, 7] = lengths[:, 2] * (cosa[:, 0] - cosa[:, 1] * cosa[:, 2]) / sing
    pbc[:, 8] = (lengths[:, 2] ** 2 - (pbc[:, 6:8] ** 2).sum(axis=1)) ** 0.5
    
    return pbc.reshape((-1, 3, 3))


# In[ ]:


def boxdim(u, selection):
    '''
    Saving and processing the box dimensions over the entire directory using the dim2pbc function.
    '''
    B = np.empty((len(u.trajectory), 6))
    for idx, f in enumerate(u.trajectory):
        B[idx] = selection.dimensions.copy()
    B = dim2pbc(B)
    return B


# In[ ]:


def nojump(X, B, ref=None):
    """
    Remove periodic boundary condition (PBC) jumps from a molecular dynamics trajectory.

    This function unwraps the trajectory `X` using the provided simulation box definitions `B`,
    correcting for discontinuities caused by particles crossing periodic boundaries. The result
    is a continuous trajectory in real space. The first frame or a user-defined reference frame
    is used as the starting point for unwrapping.

    Parameters
    ----------
    X : ndarray, shape (n_frames, n_atoms, 3)
        Molecular dynamics trajectory positions in Cartesian coordinates.
    B : ndarray, shape (n_frames, 3, 3)
        Simulation box matrices corresponding to each frame in `X`.
    ref : ndarray, shape (n_atoms, 3), optional
        Reference frame for unwrapping. If None, the first frame of `X` is used.

    Returns
    -------
    X_unwrapped : ndarray, shape (n_frames, n_atoms, 3)
        The unwrapped trajectory with periodic jumps removed and positions expressed
        in real Cartesian space.

    Notes
    -----
    - The trajectory is converted to fractional coordinates relative to the box,
      corrected for jumps, and then converted back to Cartesian coordinates.
    - This function assumes orthorhombic or triclinic boxes defined per-frame in `B`.
    """
    if ref is None:
        ref = X[0]
    ref = ref @ np.linalg.inv(B[0])
    
    X = X @ np.linalg.inv(B)
    X[1:] = np.diff(X, axis=0) + 0.5
    X[0] -= ref - 0.5
    X -= np.floor(X)
    X[0] += ref
    X = np.cumsum(X - 0.5, axis=0) @ B

    return X


# In[ ]:


def force_vectors(X,F):
    '''
    This function connects the coordinate data from the trajectory to the coordinates of the original structure in the pdb file,
    after rotating the frames onto a reference set to the first frame the covariance, eigenvectors and eigenvalues are determined.
    For the starting points of the vectors the atoms from the pdb file are taken minus the PC of interest * the corresponding eigenvalue
    (divided by a number to keep it in the right perspective) This is then loaded ont 
    '''
    Q = np.array(cmd.get_model('sele').get_coord_list())
    X -= X.mean(axis=1, keepdims=True)
    ref = X[0]
    U, L, V = np.linalg.svd(ref.T @ X / len(ref))
    X = X @ (U @ V).transpose((0, 2, 1))
    F = F @ (U @ V).transpose((0, 2, 1))
    
    C = (F[:,:,:,None] * F[:,:, None, :]).mean(axis=0)
    vals, vecs = np.linalg.eigh(C)
    
    starts1 = Q - vecs[:,:,2] * (vals[2] / 200000) ** 0.5
    ends1 = Q + vecs[:,:,2] * (vals[2] / 200000) ** 0.5
    forcevectors1 = [u for start,end in zip(starts1, ends1) for u in [9.0, *start, *end, 0.1 , 1,1,1,1,0,0]]
    
    starts2 = Q - vecs[:,:,1] * (vals[1] / 200000) ** 0.5
    ends2 = Q + vecs[:,:,1] * (vals[1] / 200000) ** 0.5
    forcevectors2 = [u for start,end in zip(starts2, ends2) for u in [9.0, *start, *end, 0.1 , 1,1,1,0,1,0]]
    
    starts3 = Q - vecs[:,:,0] * (vals[0] / 200000) ** 0.5
    ends3 = Q + vecs[:,:,0] * (vals[0] / 200000) ** 0.5
    forcevectors3 = [u for start,end in zip(starts3, ends3) for u in [9.0, *start, *end, 0.1 , 1,1,1,0,0,1]]
    
    cmd.load_cgo(forcevectors1, 'forcevectors1')
    cmd.load_cgo(forcevectors2, 'forcevectors2')
    cmd.load_cgo(forcevectors3, 'forcevectors3')


# In[ ]:


def movement_vectors(X):
    Q = np.array(cmd.get_model('sele').get_coord_list())
    X -= X.mean(axis=1, keepdims=True)
    ref = X[0]
    U, L, V = np.linalg.svd(ref.T @ X / len(ref))
    X = X @ (U @ V).transpose((0, 2, 1))
    
    Y = X - X.mean(axis=0)
    A = ((Y[:,:,:,None] * Y[:,:, None, :]).mean(axis=0) / len(X)) *10
    vals2, vecs2 = np.linalg.eigh(A)
    
    starts4 = Q - vecs2[:,:,2] * (vals2[2]) ** 0.5
    ends4 = Q + vecs2[:,:,2] * (vals2[2]) ** 0.5
    movement1 = [u for start,end in zip(starts4, ends4) for u in [9.0, *start, *end, 0.1 , 1,1,1,0.5,0,0]]
    
    starts5 = Q - vecs2[:,:,1] * (vals2[1]) ** 0.5
    ends5 = Q + vecs2[:,:,1] * (vals2[1]) ** 0.5
    movement2 = [u for start,end in zip(starts5, ends5) for u in [9.0, *start, *end, 0.1 , 1,1,1,0,0.5,0]]
    
    starts6 = Q - vecs2[:,:,0] * (vals2[0]) ** 0.5
    ends6 = Q + vecs2[:,:,0] * (vals2[0]) ** 0.5
    movement3 = [u for start,end in zip(starts6, ends6) for u in [9.0, *start, *end, 0.1 , 1,1,1,0,0,0.5]]
    
    cmd.load_cgo(movement1, 'Movements1')
    cmd.load_cgo(movement2, 'Movements2')
    cmd.load_cgo(movement3, 'Movements3')


# In[ ]:


B = boxdim(u, protein)
X = nojump(X,B, ref.atoms.positions)
movement_vectors(X)
force_vectors(X,F)

