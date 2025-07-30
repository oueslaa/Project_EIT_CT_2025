import numpy as np
import os
import scipy
import nibabel as nib
from nibabel.affines import apply_affine
import pymeshlab as mlab
import plotly.graph_objects as go
from visualise_MC import view_mesh_3D
import trimesh
import warnings
try:
    from pymeshlab import Percentage as Percentage
except:
    from pymeshlab import PercentageValue as Percentage
import matplotlib.pyplot as plt

### was given ##
def make_mesh_alpha(case_id,organ_name):
    """"loads a .nii.gz file and produces a mesh to visualise it.
    Outputs:
    ms: Surface mesh to visualise the volume
    data: Class from nibabel containing all data from the file
    p: array of coordinate points at which there is a non-zero value
    f: array of non-zero values from the image
    """

    base_dir=r'Data_set'   
    subject_path = os.path.join(base_dir, case_id)
    file=os.path.join(subject_path, "segmentations", organ_name)




    # Load file and get data
    data = nib.load(file)
    f = data.get_fdata()
    hdr =data.header
    M = data.affine

    # Obtain voxel coordinates
    nx,ny,nz = f.shape
    i = np.arange(0,nx,1)
    j = np.arange(0,ny,1)
    k = np.arange(0,nz,1)
    [i,j,k] = np.meshgrid(i,j,k,indexing='ij') #small error corrected
    sz = list(i.shape)
    sz.append(3)
    index_coords = np.ones(tuple(sz))
    index_coords[:,:,:,0] = i
    index_coords[:,:,:,1] = j
    index_coords[:,:,:,2] = k

    # Obtain RAS coordinates
    p = apply_affine(M,index_coords)
    x = p[:,:,:,0]
    y = p[:,:,:,1]
    z= p[:,:,:,2]
    R = M[0:3,0:3]

    # Use minimum eigenvalue of Affine transformation as length scale for alpha shape
    eig = np.linalg.eigvals(R)
    alpha = np.min(np.max(np.abs(eig)))

    # Filter out non-zero entries
    i,j,k = np.where(f > 0)
    f = f[i,j,k]
    x = x[i,j,k]
    y = y[i,j,k]
    z = z[i,j,k]
    p = np.vstack((x,y,z))
    p = np.array(p).T

    if np.var(f) != 0:
        warnings.warn(f'{file} is not a binary image')


    # Obtain alpha shape
    ms = mlab.MeshSet()
    m = mlab.Mesh(vertex_matrix = p)
    ms.add_mesh(m)
    bbox = ms.get_geometric_measures()['bbox']
    diag = bbox.diagonal()
    ms.generate_alpha_shape(alpha = Percentage(100*alpha/diag),filtering =1) 


    #the mesh obtained untill now could contain small holes, and edges shared by more than two triangles.
    #this prevents from calculation the volume directly,and calculating the volume manually is not preferable
    #The next filter creates a "shell"over the surface to deal with these irrigularities.with MC, it is automatically done with the scikit-image function



    alpha_frac=0.02 #Lets the wrap settle into concavities down to about 2% of the object’s size so we don’t lose small anatomical details
    offset_frac=0.0001 #Gives the shell a thickness to prevent any faces from folding  on themselve

    
    #with a lot  the parameters tested, MC is visually still better at capturing anaomical details
    
    ms.apply_filter('generate_alpha_wrap',alpha_fraction=alpha_frac,offset_fraction=offset_frac)
    
    m0=ms.current_mesh()
    faces=m0.face_matrix()    
    pts=m0.vertex_matrix()  
    mesh=trimesh.Trimesh(vertices=pts,faces=faces)
    vol_cm3=abs(mesh.volume)/1000
    surf_cm2=mesh.area/100



    return mesh,pts,faces,vol_cm3,surf_cm2

