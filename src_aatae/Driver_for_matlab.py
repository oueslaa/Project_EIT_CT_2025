import os, pickle
import numpy as np
from scipy.io import savemat
from pyeit.mesh.utils import edge_list

def export_meshes_to_mat(cache_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for fname in os.listdir(cache_dir):
        if not fname.endswith("_mesh.pkl"):
            continue
        case_id = fname.replace("_mesh.pkl", "")
        with open(os.path.join(cache_dir, fname), "rb") as f:
            data = pickle.load(f)
        node    = data["node"]      # (Nnodes,2)
        element = data["element"]   # (Nelements,3)
        el_pos  = data["el_pos"]    # (n_el,) indices   
        sigma_OOEIT = data["perm_OOEIT"]     # (Nelements,) conductivités
        # Reconstruire E comme précédemment
        boundary = edge_list(element)
        E = []
        for center in el_pos:
            segs = [e for e in boundary if center in e]
            E.append(np.vstack(segs))

        # --- Ici, on construit matdict et on y injecte sigma ---
        matdict = {
            "g":       node,
            "H":       element,
            "sigma":   data["perm_OOEIT"]  # <— on ajoute σ par élément
        }
        # Ajout dynamique des E1, E2, ..., En
        for ℓ, segs in enumerate(E, start=1):
            matdict[f"E{ℓ}"] = segs

        # Sauvegarde
        out_path = os.path.join(out_dir, f"{case_id}_mesh.mat")
        savemat(out_path, matdict, oned_as="column")
        print(f"→ Exporté {out_path}")

if __name__=="__main__":
    export_meshes_to_mat(cache_dir="meshes_cache_sujets",out_dir="meshes_mat_sujets")
    export_meshes_to_mat(cache_dir="meshes_cache",out_dir="meshes_mat")