




import nibabel as nib
import numpy as np
from scipy.ndimage import zoom
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from scipy.ndimage import binary_closing, binary_fill_holes,label
from scipy.interpolate import interp1d
from shapely.geometry import Polygon, Point
from skimage import measure
import imageio.v2 as imageio
from scipy.sparse import csr_matrix, bmat
from scipy.spatial import Delaunay
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix
from scipy.sparse import vstack, eye
from scipy.sparse import diags


from pyeit.mesh import PyEITMesh
import pyeit.mesh as mesh
from pyeit.eit.interp2d import sim2pts
from pyeit.eit.protocol import PyEITProtocol
from pyeit.visual.plot    import create_mesh_plot, create_plot
from pyeit.mesh.shape import *
from pyeit.eit.fem    import EITForward,calculate_ke,assemble

import sys, os
import pickle
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from Extract_skin_mask import *
from help_functions import *
from EIT_sim import *

from pyeit.mesh.utils import edge_list
from pyeit.mesh import PyEITMesh




def mse(a, b):
    return np.mean((a - b) ** 2)

def solve_and_plot(mesh_obj, mask2d, perm, protocol, zeta_val, mode, title):
    """
    Calcule le forward pour mesh_obj, affiche le mask, le mesh et le signal,
    et renvoie le vecteur de solution sol_all.
    """
    # Forward
    solver  = CEMFEM(mesh_obj, protocol, ref_node=0)
    sol_all = solver.solve_CEM(perm=perm,
                               zeta=np.ones(mesh_obj.n_el) * zeta_val,
                               mode=mode)

    # Affichage
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(14, 5))
    ax0.imshow(mask2d, cmap='viridis')
    ax0.set_title(f"{title}\nMask 2D")
    ax0.axis("off")

    create_mesh_plot(
        ax1, mesh_obj,
        electrodes=mesh_obj.el_pos,
        coordinate_labels="radiological",
        marker_text_kwargs={"color": "red", "fontsize": 6}
    )
    ax1.set_title(f"{title}\nMesh")
    ax1.axis("off")
    ax1.set_aspect("equal")

    ax2.plot(sol_all, '-o', markersize=4)
    ax2.set_title(f"{title}\nSignal ({mode})")
    ax2.set_xlabel("Index d'électrode")
    ax2.set_ylabel("Potentiel (V)" if mode == "current" else "Courant (A)")
    ax2.grid(True)

    plt.tight_layout()
    plt.show()

    return sol_all

def process_caches(
    ref_cache_dir,      # dossier meshes_cache_sujets
    target_cache_dir,   # dossier meshes_cache
    zeta_val,
    mode,
    dist_exc,
    step_meas=1
):
    # --- 1) Choix du mesh de référence (le premier .pkl de ref_cache_dir) ---
    ref_files = sorted(f for f in os.listdir(ref_cache_dir) if f.endswith("_mesh.pkl"))
    if not ref_files:
        raise RuntimeError("Aucun mesh de référence trouvé dans " + ref_cache_dir)
    ref_name = ref_files[0]
    ref_path = os.path.join(ref_cache_dir, ref_name)
    case_ref = ref_name.replace("_mesh.pkl", "")

    # Chargement
    with open(ref_path, "rb") as f:
        data = pickle.load(f)
    mask_ref = data["mask"]
    node_ref    = data["node"]
    elem_ref    = data["element"]
    perm0_ref   = data["perm"]
    elpos_ref   = data["el_pos"]

    # Reconstruire PyEITMesh
    mesh_ref = PyEITMesh(
        node=node_ref,
        element=elem_ref,
        el_pos=elpos_ref
    )
    mesh_ref.perm = perm0_ref

    # Protocole
    protocol_ref = set_protocol(n_el=mesh_ref.n_el,
                                dist_exc=dist_exc,
                                step_meas=step_meas)

    # Solve & Plot de la référence
    print(f"=== Référence : {case_ref} ===")
    sol_ref = solve_and_plot(
        mesh_ref, mask_ref,
        perm0_ref, protocol_ref,
        zeta_val, mode,
        title=f"{case_ref} (réf)"
    )

    # --- 2) Parcours de tous les meshes du dossier target_cache_dir ---
    target_files = sorted(f for f in os.listdir(target_cache_dir) if f.endswith("_mesh.pkl"))
    for tf in target_files:
        case_id = tf.replace("_mesh.pkl", "")
        path    = os.path.join(target_cache_dir, tf)
        with open(path, "rb") as f:
            data = pickle.load(f)

        mask2d   = data["mask"]
        node     = data["node"]
        element  = data["element"]
        perm0    = data["perm"]
        el_pos   = data["el_pos"]

        mesh_obj = PyEITMesh(node=node, element=element, el_pos=el_pos)
        mesh_obj.perm = perm0

        protocol = set_protocol(n_el=mesh_obj.n_el,
                                dist_exc=dist_exc,
                                step_meas=step_meas)

        print(f"--- Cas {case_id} vs {case_ref} ---")
        sol_cur = solve_and_plot(
            mesh_obj, mask2d,
            perm0, protocol,
            zeta_val, mode,
            title=case_id
        )

        # Comparaison
        err = mse(sol_cur, sol_ref)
        print(f"MSE entre {case_id} et {case_ref} : {err:.3e}\n")

if __name__ == "__main__":
    process_caches(
        ref_cache_dir    = "meshes_cache_sujets",
        target_cache_dir = "meshes_cache",
        zeta_val         = 1e-6,
        mode             = "current",
        dist_exc         = 1,
        step_meas        = 1
    )