
import numpy as np
from pyeit.mesh import PyEITMesh
import pyeit.mesh as mesh
from pyeit.eit.interp2d import sim2pts
from pyeit.eit.protocol import PyEITProtocol
from pyeit.visual.plot    import create_mesh_plot, create_plot
from pyeit.mesh.shape import *
from pyeit.eit.fem    import EITForward,calculate_ke,assemble
import argparse

import sys, os
import pickle
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from Exctract_skin_mask import *
from help_functions import *
from EIT_CEM import *
from help_functions_post_process import *
from pyeit.mesh.utils import edge_list
from pyeit.mesh import PyEITMesh

from skimage.morphology import disk
from scipy.ndimage import binary_closing, binary_fill_holes




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
    case_id,            # liste de strings e.g. ["s0011","s0012",...]
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
    body_poly = Extract_contour(
    mask2d,
    closing_radius=5,
    opening_radius=3,
    simplify_tol=0.01)
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


    print(f"mesh de la slice {slice_index} géneré !! ")
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

    print(f"  → Mesh de la slice {slice_index} enregistré dans {out_path}\n")

# exemple
# zslice_values=[318,92,324,359,180,346,325,317,347,369]
# CASES_LIBRARY =["s0028","s0691","s0708","s0720","s0070","s0613","s0561","s0456","s0507","s0369"]
if __name__ == "__main__":
    ORGANS = {
    "adrenal_gland": [
        "adrenal_gland_left.nii.gz",
        "adrenal_gland_right.nii.gz",
    ],
    "atrial_appendage": [
        "atrial_appendage_left.nii.gz",
        "atrial_appendage_right.nii.gz",
    ],
    "autochthon": [
        "autochthon_left.nii.gz",
        "autochthon_right.nii.gz",
    ],
    "brain": [
        "brain.nii.gz",
    ],
    "aorta": [
        "aorta.nii.gz",
    ],
    "brachiocephalic_trunk": [
        "brachiocephalic_trunk.nii.gz",
    ],
    "brachiocephalic_vein": [
        "brachiocephalic_vein_left.nii.gz",
        "brachiocephalic_vein_right.nii.gz",
    ],
    "clavicula": [
        "clavicula_left.nii.gz",
        "clavicula_right.nii.gz",
    ],
    "colon": [
        "colon.nii.gz",
    ],
    "common_carotid_artery": [
        "common_carotid_artery_left.nii.gz",
        "common_carotid_artery_right.nii.gz",
    ],
    "costal_cartilages": [
        "costal_cartilages.nii.gz",
    ],
    "duodenum": [
        "duodenum.nii.gz",
    ],
    "esophagus": [
        "esophagus.nii.gz",
    ],
    "femur": [
        "femur_left.nii.gz",
        "femur_right.nii.gz",
    ],
    "gallbladder": [
        "gallbladder.nii.gz",
    ],
    "gluteus_maximus": [
        "gluteus_maximus_left.nii.gz",
        "gluteus_maximus_right.nii.gz",
    ],
    "gluteus_medius": [
        "gluteus_medius_left.nii.gz",
        "gluteus_medius_right.nii.gz",
    ],
    "gluteus_minimus": [
        "gluteus_minimus_left.nii.gz",
        "gluteus_minimus_right.nii.gz",
    ],
    "heart": [
        "heart.nii.gz",
    ],
    "hip": [
        "hip_left.nii.gz",
        "hip_right.nii.gz",
    ],
    "humerus": [
        "humerus_left.nii.gz",
        "humerus_right.nii.gz",
    ],
    "iliac_artery": [
        "iliac_artery_left.nii.gz",
        "iliac_artery_right.nii.gz",
    ],
    "iliac_vena": [
        "iliac_vena_left.nii.gz",
        "iliac_vena_right.nii.gz",
    ],
    "iliopsoas": [
        "iliopsoas_left.nii.gz",
        "iliopsoas_right.nii.gz",
    ],
    "inferior_vena_cava": [
        "inferior_vena_cava.nii.gz",
    ],
    "kidney": [
        "kidney_left.nii.gz",
        "kidney_right.nii.gz",
    ],
    "kidney_cyst": [
        "kidney_cyst_left.nii.gz",
        "kidney_cyst_right.nii.gz",
    ],
    "lung": [
        "lung_lower_lobe_left.nii.gz",
        "lung_lower_lobe_right.nii.gz",
        "lung_middle_lobe_right.nii.gz",
        "lung_upper_lobe_left.nii.gz",
        "lung_upper_lobe_right.nii.gz",
    ],
    "pancreas": [
        "pancreas.nii.gz",
    ],
    "pulmonary_vein": [
        "pulmonary_vein.nii.gz",
    ],
    "portal_and_splenic_vein": [
        "portal_vein_and_splenic_vein.nii.gz",
    ],
    "prostate": [
        "prostate.nii.gz",
    ],
    "ribs": [
        "rib_left_1.nii.gz",  "rib_left_2.nii.gz",  "rib_left_3.nii.gz",
        "rib_left_4.nii.gz",  "rib_left_5.nii.gz",  "rib_left_6.nii.gz",
        "rib_left_7.nii.gz",  "rib_left_8.nii.gz",  "rib_left_9.nii.gz",
        "rib_left_10.nii.gz", "rib_left_11.nii.gz", "rib_left_12.nii.gz",
        "rib_right_1.nii.gz", "rib_right_2.nii.gz", "rib_right_3.nii.gz",
        "rib_right_4.nii.gz", "rib_right_5.nii.gz", "rib_right_6.nii.gz",
        "rib_right_7.nii.gz", "rib_right_8.nii.gz", "rib_right_9.nii.gz",
        "rib_right_10.nii.gz","rib_right_11.nii.gz","rib_right_12.nii.gz",
    ],
    "sacrum": [
        "sacrum.nii.gz",
    ],
    "scapula": [
        "scapula_left.nii.gz",
        "scapula_right.nii.gz",
    ],
    "small_bowel": [
        "small_bowel.nii.gz",
    ],
    "spinal_cord": [
        "spinal_cord.nii.gz",
    ],
    "spleen": [
        "spleen.nii.gz",
    ],
    "sternum": [
        "sternum.nii.gz",
    ],
    "stomach": [
        "stomach.nii.gz",
    ],
    "subclavian_artery": [
        "subclavian_artery_left.nii.gz",
        "subclavian_artery_right.nii.gz",
    ],
    "superior_vena_cava": [
        "superior_vena_cava.nii.gz",
    ],
    "thyroid_gland": [
        "thyroid_gland.nii.gz",
    ],
    "trachea": [
        "trachea.nii.gz",
    ],
    "urinary_bladder": [
        "urinary_bladder.nii.gz",
    ],
    "vertebrae": [
        "vertebrae_C1.nii.gz",  "vertebrae_C2.nii.gz",  "vertebrae_C3.nii.gz",
        "vertebrae_C4.nii.gz",  "vertebrae_C5.nii.gz",  "vertebrae_C6.nii.gz",
        "vertebrae_C7.nii.gz",  "vertebrae_L1.nii.gz",  "vertebrae_L2.nii.gz",
        "vertebrae_L3.nii.gz",  "vertebrae_L4.nii.gz",  "vertebrae_L5.nii.gz",
        "vertebrae_S1.nii.gz",  "vertebrae_T1.nii.gz",  "vertebrae_T2.nii.gz",
        "vertebrae_T3.nii.gz",  "vertebrae_T4.nii.gz",  "vertebrae_T5.nii.gz",
        "vertebrae_T6.nii.gz",  "vertebrae_T7.nii.gz",  "vertebrae_T8.nii.gz",
        "vertebrae_T9.nii.gz",  "vertebrae_T10.nii.gz", "vertebrae_T11.nii.gz",
        "vertebrae_T12.nii.gz",
    ],
    }
    p = argparse.ArgumentParser()
    p.add_argument("--case",           required=True)
    p.add_argument("--slice-global",   type=int, required=True)
    p.add_argument("--slice-relative", type=int, required=True)
    args = p.parse_args()

    case_id     = args.case
    slice_rel   = args.slice_relative   # <— ici on passe REL et non GLOBAL

    # 
    ct_path = os.path.join("Data_set", case_id, "ct.nii.gz")
    seg_dir = os.path.join("Data_set", case_id, "segmentations")
    crop_ct, z_min, z_max, ct = compute_z_bounds(
        ct_path, seg_dir,
        top="vertebrae_T3", bottom="vertebrae_L3", margin=0
    )

    
    present_organs = organs_present_in_crop(seg_dir, z_min, z_max, ORGANS)
    mask3d = Get_soft_tissue_mask(case_id, present_organs, ORGANS)

    
    build_and_store_meshes(
        case_id         = case_id,
        base_dir        = "Data_set",
        organs          = ORGANS,
        slice_index     = slice_rel,    # <— ici on passe REL et non GLOBAL
        mask3d          = mask3d,
        present_organs  = present_organs,
        n_el            = 16,
        h0              = 0.04,
        output_dir      = "meshes_cache_sujets",
    )   

    # for case_id, slice_index in zip(CASES_LIBRARY,zslice_values):

    #     build_and_store_meshes(                         ##
    #         case_id=case_id,
    #         base_dir="Data_set",
    #         organs=ORGANS,
    #         slice_index=slice_index,
    #         n_el=16,
    #         h0=0.04,
    #         output_dir="meshes_cache"
    #     )



    print(" ")