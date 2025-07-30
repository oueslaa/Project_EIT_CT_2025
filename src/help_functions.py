from shapely.geometry import Polygon, Point
import numpy as np
import pyeit.eit.protocol as protocol


from pyeit.mesh.shape import *
from scipy.spatial import Delaunay, cKDTree
import nibabel as nib


## Help functions for EIT##


def make_fd_body(body_poly: Polygon):
    """
    Returns a function fd_body(pts) that takes an (M,2) array of (x,y) points
    and returns a length‐M array of signed distances to the polygon boundary:
       fd_body[i] < 0 if point i is iniside body_poly,
       fd_body[i] = 0 if point i is exactly on the boundary,
       fd_body[i] > 0 if point i is outside body_poly.
       To be used directly in mesh.creaate, as the "fd" variable
    """
    
    def fd_body(pts):
        """return signed distance of polygon"""
        pts_ = [Point(p) for p in pts]
        # calculate signed distance
        dist = [body_poly.exterior.distance(p) for p in pts_]
        sign = np.sign([-int(body_poly.contains(p)) + 0.5 for p in pts_])

        return sign * dist
    

    return fd_body

def compute_element_labels_affine(mask, mesh, case_id, slice_index=None):
    """
    Pour chaque centre de triangle (en mm), retourne le label du mask 2D.

    - mask        : array 2D (H x W)
    - mesh.elem_centers : array (Ne,2) en coordonnées réelles (mm)
    - case_id     : pour charger Data_set/{case_id}/ct.nii.gz
    - slice_index : si elem_centers est (Ne,2), donne le z (int) de coupe

    Sortie :
    - labels : array (Ne,) de labels issus de mask[row, col]
    """
    # 1) Récupère les centres en mm
    centers = mesh.elem_centers  # shape (Ne,2)

    # 2) Construis les homogènes (Ne,4)
    N = centers.shape[0]
    ones = np.ones((N,1))
    if centers.shape[1] == 2:
        if slice_index is None:
            raise ValueError("slice_index requis si elem_centers est (N,2)")
        zs = np.full((N,1), slice_index)
        homog = np.hstack([centers, zs, ones])
    elif centers.shape[1] == 3:
        homog = np.hstack([centers, ones])
    else:
        raise ValueError("elem_centers doit être (N,2) ou (N,3)")

    # 3) Applique l’inverse de l’affine du CT pour retomber en voxels
    img     = nib.load(f"Data_set/{case_id}/ct.nii.gz")
    inv_aff = np.linalg.inv(img.affine)
    vox     = homog.dot(inv_aff.T)  # (Ne,4) -> (i,j,k,1)
    xi      = np.round(vox[:, 0]).astype(int)
    yi      = np.round(vox[:, 1]).astype(int)

    # 4) Comme imshow(mask) met l’origine en haut à gauche,
    #    il faut inverser l’indice de ligne :
    rows = mask.shape[0] - 1 - yi
    cols = xi

    # 5) Clip et retourne
    rows = np.clip(rows, 0, mask.shape[0]-1)
    cols = np.clip(cols, 0, mask.shape[1]-1)
    return mask[rows, cols]

def set_protocol(n_el,dist_exc,step_meas):
    """
    Set the protocol for EIT measurements.
    Returns a protocol object with the specified parameters.
    
    n_el =    # Number of electrodes
    dist_exc  # Distance between excitation electrodes
    step_meas # Step size for measurements

    """
    return protocol.create(n_el, dist_exc=dist_exc, step_meas=step_meas,parser_meas="std")

def set_condu_eit(mask2d,mesh,case_id,slice_index,len_organs):
    """
    Set the conductivity values for the EIT mesh.
    Returns an array of conductivity values for each element in the mesh.
    mesh: EIT mesh object, of type pyeit.mesh.Mesh
    condu_body: Conductivity value for the body
    condu_lung: Conductivity value for the lung
    labels_elems: Array of labels for each element in the mesh
    Ne = labels_elems.shape[0](number of triangles in the mesh)
    Output:
    perm: Array of conductivity values for each element in the mesh
    perm0: Array of initial conductivity values (all ones)(Reference conductivity)
    """
    tri_labels=compute_element_labels_affine(mask2d,mesh,case_id, slice_index=slice_index)
    Ne=mesh.element.shape[0]
    
    perm= np.ones(Ne)

    # applique les valeurs
    for j in range(1,150):

        perm[tri_labels==j] = j
    return perm

def compute_node_labels_affine(mask, nodes, case_id, slice_index):
    """
    Pour chaque nœud (en mm), retourne le label du mask 2D.

    - mask       : array 2D (H x W)
    - nodes      : array (Nnodes,2) ou (Nnodes,3) coords réelles (mm)
    - case_id    : pour charger Data_set/{case_id}/ct.nii.gz
    - slice_index: indice de coupe Z, requis si nodes.shape[1]==2

    Sortie :
    - labels_nœuds : array (Nnodes,) de labels issus de mask[row, col]
    """
    # 1) Récupère l’affine inverse du CT
    img     = nib.load(f"Data_set/{case_id}/ct.nii.gz")
    inv_aff = np.linalg.inv(img.affine)

    # 2) Prépare les coordonnées homogènes (Nnodes,4)
    nodes_mm = nodes[:, :2]
    N = nodes_mm.shape[0]
    ones = np.ones((N,1))
    if nodes.shape[1] == 2:
        zs = np.full((N,1), slice_index)
        homog = np.hstack([nodes_mm, zs, ones])
    elif nodes.shape[1] == 3:
        homog = np.hstack([nodes_mm, nodes[:,2:3], ones])
    else:
        raise ValueError("nodes must be (N,2) or (N,3)")

    # 3) Passage en voxels floating
    vox = homog.dot(inv_aff.T)  # shape (N,4)
    xi  = np.round(vox[:, 0]).astype(int)
    yi  = np.round(vox[:, 1]).astype(int)

    # 4) Inversion de l’axe Y pour matcher imshow(origin='lower')
    rows = mask.shape[0] - 1 - yi
    cols = xi

    # 5) Clip pour rester dans l’image
    rows = np.clip(rows, 0, mask.shape[0]-1)
    cols = np.clip(cols, 0, mask.shape[1]-1)

    # 6) Extraction
    return mask[rows, cols]


def set_condu_nodes(mask, mesh,case_id, slice_index,len_organs):
    """
    Retourne un vecteur de conductivité par nœud (taille mesh.n_nodes),
    à partir du masque et des valeurs pour poumon et corps.
    """
    # récupère les labels 0/1/2 pour chaque nœud
    nodes_mm     = mesh.node[:, :2]  # (Nnodes,2)
    node_labels  = compute_node_labels_affine(mask, nodes_mm, case_id, slice_index)
    Nnodes       = mesh.node.shape[0]

    # initialise à 1 (milieu de référence)
    perm_nodes = np.ones(Nnodes)

    # applique les valeurs
    for j in range(1,150):

        perm_nodes[node_labels==j] = j

    return perm_nodes











def sample_electrodes(poly, n_el):
    perim = poly.exterior.length
    dists = np.linspace(0, perim, num=n_el, endpoint=False)
    return np.array([poly.exterior.interpolate(d).coords[0] for d in dists])

def find_closest_node_indices(nodes, pts):
    tree = cKDTree(nodes)
    _, idx = tree.query(pts)
    return idx
