
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
import nibabel as nib

from pyeit.mesh import PyEITMesh
import pyeit.mesh as mesh
from pyeit.eit.interp2d import sim2pts
from pyeit.eit.protocol import PyEITProtocol
from pyeit.visual.plot    import create_mesh_plot, create_plot
from pyeit.mesh.shape import *
from pyeit.mesh.utils import edge_list,check_ccw, check_order,to_polar
import sys, os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from Extract_skin_mask import *
from help_functions import *
from EIT_sim import *

import pickle




class ThinPlateSpline2D:
    def __init__(self):
        self.w = None
        self.a = None
        self.src = None

    def _kernel(self, r):
        with np.errstate(divide='ignore', invalid='ignore'):
            K = r**2 * np.log(r)
        K[np.isnan(K)] = 0
        return K

    def fit(self, src, dst):
        M = src.shape[0]
        # assemble
        D = cdist(src, src)
        K = self._kernel(D)
        P = np.hstack([np.ones((M,1)), src])
        A = np.zeros((M+3, M+3))
        A[:M,:M]   = K
        A[:M,M:]   = P
        A[M:,:M]   = P.T
        Y = np.vstack([dst, np.zeros((3,2))])
        sol = np.linalg.solve(A, Y)
        self.w = sol[:M]
        self.a = sol[M:]
        self.src = src

    def transform(self, pts):
        D = cdist(pts, self.src)
        U = self._kernel(D)
        non_affine = U.dot(self.w)
        P = np.hstack([np.ones((pts.shape[0],1)), pts])
        affine     = P.dot(self.a)
        return non_affine + affine


def extract_ordered_boundary(mesh,bars):
    """
    1) Repère les arêtes de contour (edge_list)
    2) Construit la liste ordonnée des nœuds sur ce contour.
    Retourne :
      - bnd_idx : array de shape (K,) des indices des nœuds,  
      - bnd_pts : array de shape (K,2) leurs coordonnées (x,y).
    """
 
    print("Hello")
    # 2) construis le dictionnaire d'adjacence sur le contour
    # chaque nœud y apparaît exactement dans 2 arêtes (boucle fermée)
    uniq = np.unique(bars)
    neigh = {i: [] for i in uniq}
    for i,j in bars:
        neigh[i].append(j)
        neigh[j].append(i)

    # 3) parcours la boucle
    start = bars[0,0]
    boundary = [start]
    prev = None
    curr = start
    while True:
        nbrs = neigh[curr]
        # choisis le suivant : celui qui n'est pas le précédent
        nxt = nbrs[0] if nbrs[0] != prev else nbrs[1]
        if nxt == start:
            break
        boundary.append(nxt)
        prev, curr = curr, nxt

    bnd_idx = np.array(boundary, dtype=int)
    bnd_pts = mesh.node[bnd_idx, :2]
    return bnd_idx, bnd_pts


def resample_contour_by_arclength(pts, N):
    """
    Échantillonne N points régulièrement selon la longueur d'arc
    le long d'une chaîne fermée de points pts (shape K×2).

    Retour
    ------
    new_pts : array (N×2)
        Les N points rééchantillonnés.
    """
    # 1) Assurer que le contour est bien fermé en re-ajoutant le 1er point à la fin
    pts_closed = np.vstack([pts, pts[0]])
    
    # 2) Calculer les longueurs de chaque segment
    seg_vecs    = np.diff(pts_closed, axis=0)              # (K,2)
    seg_lengths = np.linalg.norm(seg_vecs, axis=1)         # (K,)
    
    # 3) Distance cumulée le long du contour
    cumdist = np.concatenate([[0], np.cumsum(seg_lengths)])  # (K+1,)
    total_length = cumdist[-1]
    
    # 4) Positions cibles régulièrement réparties sur [0, total_length)
    sample_d = np.linspace(0, total_length, N, endpoint=False)
    
    # 5) Pour chaque position d, trouver le segment et interpoler
    new_pts = np.zeros((N, 2))
    for i, d in enumerate(sample_d):
        # indice du segment précédent le point d
        idx = np.searchsorted(cumdist, d) - 1
        idx = np.clip(idx, 0, len(seg_lengths)-1)
        # paramètre t dans [0,1] sur ce segment
        t = (d - cumdist[idx]) / seg_lengths[idx]
        # interpolation linéaire
        new_pts[i] = (1 - t) * pts_closed[idx] + t * pts_closed[idx+1]
    
    return new_pts

def align_contours(src, dst):
    """
    Aligne src sur dst par rotation circulaire (et éventuellement inversion)
    pour minimiser la somme des distances au carré point-à-point.

    Parameters
    ----------
    src : (N×2) array
        Contour source échantillonné.
    dst : (N×2) array
        Contour destination échantillonné.

    Returns
    -------
    src_aligned : (N×2) array
        src après éventuel renversement et rotation pour s'aligner sur dst.
    best_shift : int
        Nombre de pas de rotation appliqué (np.roll(src, -best_shift)).
    reversed_flag : bool
        True si on a utilisé src[::-1] pour l’alignement.
    """
    N = src.shape[0]
    best_err = np.inf
    best_shift = 0
    reversed_flag = False

    # tester les deux orientations : normale et renversée
    for rev in (False, True):
        candidate = src if not rev else src[::-1]
        # tester toutes les rotations circulaires
        for k in range(N):
            rolled = np.roll(candidate, -k, axis=0)
            err = np.sum((rolled - dst)**2)
            if err < best_err:
                best_err = err
                best_shift = k
                reversed_flag = rev

    # construire la version alignée finale
    aligned = src[::-1] if reversed_flag else src
    src_aligned = np.roll(aligned, -best_shift, axis=0)
    return src_aligned, best_shift, reversed_flag


############# Post-treatment de mesh PyEIT ####################

def boundary_loop_indices_from_edges(bars):
    # renvoie les indices ordonnés du contour (CCW)
    # construire l’adjacence
    adj = {}
    for u, v in bars:
        u, v = int(u), int(v)
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)
    # parcourir le cycle
    start = bars[0, 0]
    prev, cur = None, int(start)
    loop = [cur]
    while True:
        nbrs = adj[cur]
        nxt = nbrs[0] if nbrs[0] != prev else nbrs[1]
        if nxt == start:
            break
        loop.append(nxt)
        prev, cur = cur, nxt
    return np.array(loop, dtype=int)

def ensure_ccw(nodes_xy, loop_idx):
    # calcul de l’aire par la formule du polygone (shoelace)
    pts = nodes_xy[loop_idx]
    x, y = pts[:,0], pts[:,1]
    area2 = np.dot(x, np.roll(y,-1)) - np.dot(y, np.roll(x,-1))
    # si orientation horaire on inverse
    return loop_idx[::-1] if area2 < 0 else loop_idx

def smooth_mesh_boundary(mesh: PyEITMesh, iterations:int=20, alpha:float=0.5,bars:list=None):
    """
    Lisse les nœuds frontières de `mesh` en appliquant un smoothing
    itératif le long de la boucle frontière.

    - iterations : nombre d’itérations de lissage
    - alpha      : poids du laplacian smoothing (0<α≤1)
    """
    nodes = mesh.node.copy()
    # 1) extraire et ordonner la boucle frontière
    loop = boundary_loop_indices_from_edges(bars)
    loop = ensure_ccw(nodes[:,:2], loop)

    for _ in range(iterations):
        new_nodes = nodes.copy()
        coords = nodes[loop, :2]
        N = len(loop)
        for k, idx in enumerate(loop):
            prev = coords[k-1] if k>0    else coords[-1]
            curr = coords[k]
            nxt  = coords[k+1] if k<N-1 else coords[0]
            avg  = (prev + curr + nxt)/3
            # on fait une interpolation curr → avg
            new_nodes[idx, 0:2] = (1-alpha)*curr + alpha*avg
        nodes = new_nodes

    # reconstruire le mesh
    mesh_smooth = PyEITMesh(
        node    = nodes,
        element = mesh.element,
        el_pos  = mesh.el_pos
    )
    mesh_smooth.perm = mesh.perm
    return mesh_smooth




def compute_z_bounds(ct_path, seg_dir, top, bottom, margin):
    """
    Crop en Z entre `bottom`  et `top` avec une marge de `margin` voxels.
    """
    # 
    ct_img   = nib.load(ct_path)
    ct       = ct_img.get_fdata()       
    nz       = ct.shape[2]

    
    bottom_mask_z = np.zeros(nz, dtype=bool)
    for fname in os.listdir(seg_dir):
        if fname.startswith(bottom) and fname.endswith(".nii.gz"):
            seg = nib.load(os.path.join(seg_dir, fname)).get_fdata()
            # True pour chaque coupe Z contenant bottom
            bottom_mask_z |= np.any(seg > 0, axis=(0,1))
    if not bottom_mask_z.any():
        raise RuntimeError(f"Aucun repère bottom `{bottom}` trouvé.")
    z_min = bottom_mask_z.argmax()       # premier True

    
    top_mask_z = np.zeros(nz, dtype=bool)
    for fname in os.listdir(seg_dir):
        if fname.startswith(top) and fname.endswith(".nii.gz"):
            seg = nib.load(os.path.join(seg_dir, fname)).get_fdata()
            top_mask_z |= np.any(seg > 0, axis=(0,1))
    if not top_mask_z.any():
        raise RuntimeError(f"Aucun repère top `{top}` trouvé.")
    z_max = nz - 1 - top_mask_z[::-1].argmax()  # dernier True

    
    z_min = max(z_min - margin, 0)
    z_max = min(z_max + margin, nz - 1)
    if z_min >= z_max:
        raise RuntimeError(f"z_min ({z_min}) >= z_max ({z_max})")

    
    crop_ct = ct[:, :, z_min : z_max + 1]
    return crop_ct, z_min, z_max,ct



def organs_present_in_crop(seg_dir, z_min, z_max, structures):
    """
    Détermine quels organes sont présents dans la tranche Z [z_min, z_max]
    INPUT:  le dictionnaire `structures` de toutes les structures.
    
    """
    present = []
    
    for organ, parts in structures.items():
        found = False
        for part in parts:
            seg_path = os.path.join(seg_dir, part)
            if not os.path.exists(seg_path):
                continue
            seg = nib.load(seg_path).get_fdata()
            
            crop_seg = seg[:, :, z_min:z_max+1]
            if np.any(crop_seg > 0):
                found = True
                break
        if found:
            present.append(organ)
    return present
#####################################




def Get_soft_tissue_mask(case_id, organs_presnt, ORGANS):
    """
    Récupère le masque des tissus mous à partir des organes donnés.
    """
    selected_files = [
    fname
    for struct in organs_presnt
    for fname in ORGANS.get(struct, [])
    ]
    soft_tissue, outside, _ = Create_skin_mask_bis(case_id,selected_files)
    mask3d=[soft_tissue, outside]
    for organ in organs_presnt :
        organ_mask = Create_organ_mask(case_id, ORGANS[organ])
        mask3d.append(organ_mask)
    return mask3d



def build_and_store_meshes(
    case_id,            # 
    base_dir,            # dossier racine où sont Data_set/{case_id}/ct.nii.gz
    organs,         # liste de noms de segmentation à charger
    slice_index,         # coupe à extraire
    mask3d,
    present_organs,
    n_el,                # nb électrodes pour pyeit.mesh.create
    h0,                  # résolution pour pyeit.mesh.create
    output_dir           # où écrire les .pkl
):
    # print (len(organs))
    os.makedirs(output_dir, exist_ok=True)

    print(f"==> Traitement de {case_id} slice {slice_index}")
    mask2d = Create_mask_2D(mask3d, slice_index)
    
    print(f"mask de la slice {slice_index}crée !")
    # 2) contour et fonction distance
    body_poly = Extract_contour(mask2d)
    fd_body   = make_fd_body(body_poly)

    # 3) création du mesh PyEIT en coords normalisées
    mesh_obj = mesh.create(
         n_el=n_el, h0=h0, fd=fd_body, fh=area_uniform
    )

    # 4) passe en coordonnées réelles (mm)
    mesh_obj = scale_pyeit_mesh_to_mm(mesh_obj, case_id, slice_index, mask2d)

    # 5) étiquettes éléments & conductivités
        
    perm0 = set_condu_eit(mask2d,mesh_obj,case_id,slice_index,len(organs))
    mesh_obj.perm = perm0


    perm=set_condu_nodes(mask2d, mesh_obj, case_id, slice_index,len(organs))

     

    # 6) conversion mesh.node mm→m pour la simulation SI
    mesh_obj.node[:,0] /= 1000.0
    mesh_obj.node[:,1] /= 1000.0


    print(f"mesh de la slice {slice_index} géneré !")
    edges=edge_list(mesh_obj.element   )
    # 7) sérialisation
    out_path = os.path.join(output_dir, f"{case_id}_mesh_slice_{slice_index}.pkl")
    with open(out_path, "wb") as f:
        pickle.dump({
            "mask":      mask2d,
            "node":    mesh_obj.node,
            "element": mesh_obj.element,
            "perm":    mesh_obj.perm,
            "perm_OOEIT":    perm,
            "el_pos":  mesh_obj.el_pos,
            "boundary_edges": edges,
            "present_organs": present_organs,
        }, f, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"Mesh de la slice {slice_index} enregistré dans {out_path}\n")