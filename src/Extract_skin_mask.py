import nibabel as nib
import numpy as np
from skimage import measure

import matplotlib.pyplot as plt
from scipy.ndimage import binary_closing, binary_fill_holes,label
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
from shapely.geometry import Polygon, Point
import napari      # Uncomment if you have napari installed, but it is not necessary for the main functionality

from scipy.ndimage import binary_closing, binary_fill_holes,binary_opening



##Help function 
def keep_largest_component(mask):
    """
    Conserve uniquement la plus grande composante connexe d'un masque binaire 2D.
    """
    labeled, num = label(mask)
    if num == 0:
        return mask  # rien à garder
    sizes = np.bincount(labeled.ravel())
    sizes[0] = 0  # Ignore le fond
    largest = sizes.argmax()
    mask_clean = (labeled == largest)
    return mask_clean.astype(np.uint8)

def Create_skin_mask(case_id):

    """
    INPUT:
    - case_id: ID du cas,
    - organ_parts: liste des parties de l'organe à considérer,
    - z: indice de la tranche à extraire (par exemple, 340 pour la tranche 340).
    OUTPUT:
    - lungs_mask: lungs mask,     (en 3D)
    - rest_of_body_clean: masque du reste du corps,(en 3D)
    - outside_mask_clean: masque de l'extérieur,(en 3D)
    - body: masque du corps,(en 3D)
    - ct: données de la CT, (en 3D)
    - mask_eit: la coupe en z de la CT avec  masques, pour EIT, (en 2D)
    """         
    base_dir=r'Data_set'   
    ct_path=os.path.join(base_dir,case_id,"ct.nii.gz")
    subject_path = os.path.join(base_dir, case_id)
    organ_parts = [
    "adrenal_gland_left.nii.gz",
    "adrenal_gland_right.nii.gz",
    "aorta.nii.gz",
    "atrial_appendage_left.nii.gz",
    "autochthon_left.nii.gz",
    "autochthon_right.nii.gz",
    "brachiocephalic_trunk.nii.gz",
    "brachiocephalic_vein_left.nii.gz",
    "brachiocephalic_vein_right.nii.gz",
    "brain.nii.gz",
    "clavicula_left.nii.gz",
    "clavicula_right.nii.gz",
    "colon.nii.gz",
    "common_carotid_artery_left.nii.gz",
    "common_carotid_artery_right.nii.gz",
    "costal_cartilages.nii.gz",
    "duodenum.nii.gz",
    "esophagus.nii.gz",
    "femur_left.nii.gz",
    "femur_right.nii.gz",
    "gallbladder.nii.gz",
    "gluteus_maximus_left.nii.gz",
    "gluteus_maximus_right.nii.gz",
    "gluteus_medius_left.nii.gz",
    "gluteus_medius_right.nii.gz",
    "gluteus_minimus_left.nii.gz",
    "gluteus_minimus_right.nii.gz",
    "heart.nii.gz",
    "hip_left.nii.gz",
    "hip_right.nii.gz",
    "humerus_left.nii.gz",
    "humerus_right.nii.gz",
    "iliac_artery_left.nii.gz",
    "iliac_artery_right.nii.gz",
    "iliac_vena_left.nii.gz",
    "iliac_vena_right.nii.gz",
    "iliopsoas_left.nii.gz",
    "iliopsoas_right.nii.gz",
    "inferior_vena_cava.nii.gz",
    "kidney_cyst_left.nii.gz",
    "kidney_cyst_right.nii.gz",
    "kidney_left.nii.gz",
    "kidney_right.nii.gz",
    "liver.nii.gz",
    "lung_lower_lobe_left.nii.gz",
    "lung_lower_lobe_right.nii.gz",
    "lung_middle_lobe_right.nii.gz",
    "lung_upper_lobe_left.nii.gz",
    "lung_upper_lobe_right.nii.gz",
    "pancreas.nii.gz",
    "portal_vein_and_splenic_vein.nii.gz",
    "prostate.nii.gz",
    "pulmonary_vein.nii.gz",
    "rib_left_1.nii.gz",
    "rib_left_2.nii.gz",
    "rib_left_3.nii.gz",
    "rib_left_4.nii.gz",
    "rib_left_5.nii.gz",
    "rib_left_6.nii.gz",
    "rib_left_7.nii.gz",
    "rib_left_8.nii.gz",
    "rib_left_9.nii.gz",
    "rib_left_10.nii.gz",
    "rib_left_11.nii.gz",
    "rib_left_12.nii.gz",
    "rib_right_1.nii.gz",
    "rib_right_2.nii.gz",
    "rib_right_3.nii.gz",
    "rib_right_4.nii.gz",
    "rib_right_5.nii.gz",
    "rib_right_6.nii.gz",
    "rib_right_7.nii.gz",
    "rib_right_8.nii.gz",
    "rib_right_9.nii.gz",
    "rib_right_10.nii.gz",
    "rib_right_11.nii.gz",
    "rib_right_12.nii.gz",
    "sacrum.nii.gz",
    "scapula_left.nii.gz",
    "scapula_right.nii.gz",
    "skull.nii.gz",
    "small_bowel.nii.gz",
    "spinal_cord.nii.gz",
    "spleen.nii.gz",
    "sternum.nii.gz",
    "stomach.nii.gz",
    "subclavian_artery_left.nii.gz",
    "subclavian_artery_right.nii.gz",
    "superior_vena_cava.nii.gz",
    "thyroid_gland.nii.gz",
    "trachea.nii.gz",
    "urinary_bladder.nii.gz",
    "vertebrae_C1.nii.gz",
    "vertebrae_C2.nii.gz",
    "vertebrae_C3.nii.gz",
    "vertebrae_C4.nii.gz",
    "vertebrae_C5.nii.gz",
    "vertebrae_C6.nii.gz",
    "vertebrae_C7.nii.gz",
    "vertebrae_L1.nii.gz",
    "vertebrae_L2.nii.gz",
    "vertebrae_L3.nii.gz",
    "vertebrae_L4.nii.gz",
    "vertebrae_L5.nii.gz",
    "vertebrae_S1.nii.gz",
    "vertebrae_T1.nii.gz",
    "vertebrae_T2.nii.gz",
    "vertebrae_T3.nii.gz",
    "vertebrae_T4.nii.gz",
    "vertebrae_T5.nii.gz",
    "vertebrae_T6.nii.gz",
    "vertebrae_T7.nii.gz",
    "vertebrae_T8.nii.gz",
    "vertebrae_T9.nii.gz",
    "vertebrae_T10.nii.gz",
    "vertebrae_T11.nii.gz",
    "vertebrae_T12.nii.gz"
    ]

    shape = None
    mask_total = None

    for part in organ_parts:
        seg_path=os.path.join(subject_path, "segmentations",part)
        if not os.path.exists(seg_path):
            print(f"{seg_path} n'existe pas, on saute.")
            continue
        mask = nib.load(seg_path).get_fdata().astype(np.uint8)
        if mask_total is None:
            mask_total = np.zeros_like(mask)
        mask_total = np.logical_or(mask_total, mask)


    mask_total=mask_total.astype(np.uint8)

    ct=nib.load(ct_path).get_fdata()
    body_mask=(ct>-200)# we condider that the outside has small density, so we can use a threshold
                            # could be adjusted based on the subject                        
    body_mask=keep_largest_component(body_mask).astype(np.uint8)

    body_mask_filled = binary_fill_holes(body_mask).astype(np.uint8)  #some zones inside the body may not be filled,in case we only work one or few organs
    rest_of_body_mask_filled=(body_mask_filled & (~mask_total)).astype(np.uint8)
    rest_of_body_mask_filled=binary_fill_holes(rest_of_body_mask_filled).astype(np.uint8)       #Same
    outside_mask_filled=(~np.logical_or(rest_of_body_mask_filled,mask_total)).astype(np.uint8)


    outside_mask_clean = keep_largest_component(outside_mask_filled)
    rest_of_body_clean = (~np.logical_or(outside_mask_clean,mask_total)).astype(np.uint8)
    body=np.logical_or(rest_of_body_clean,mask_total).astype(np.uint8)

    return rest_of_body_clean,outside_mask_clean,ct






def Extrat_skin(case_id):
    """
    INPUT:
    - case_id: ID du cas,
    OUTPUT:
    -ct_skin: données de la CT qui ne contient que la peau, (en 3D), de type .nii.gz
    - Enregistre le masque de la peau dans le dossier segmentations du cas.
    """
    base_dir=r'Data_set'   
    ct_path=os.path.join(base_dir,case_id,"ct.nii.gz")
    subject_path = os.path.join(base_dir, case_id,"segmentations")
    output_path = os.path.join(subject_path, "skin.nii.gz")
    ct=nib.load(ct_path)
    ct_data = ct.get_fdata()
    header = ct.header.copy()
    affine = nib.load(ct_path).affine


    skin_mask,outside_mask_clean,_= Create_skin_mask(case_id)
    ct_skin = nib.Nifti1Image(skin_mask, affine, header=header)  # Create a new Nifti1Image with the skin mask
    nib.save(ct_skin, output_path)
    return skin_mask,outside_mask_clean


# def Create_organ_mask(case_id,organ_parts):
#     """
#     INPUT:
#     - case_id: ID du cas,
#     - organ_parts: liste des parties de l'organe à considérer,
#     OUTPUT:
#     - mask_total: masque de l'organe, (en 3D)
#     """

#     base_dir=r'Data_set'   
#     ct_path=os.path.join(base_dir,case_id,"ct.nii.gz")
#     subject_path = os.path.join(base_dir, case_id)

#     shape = None
#     mask_total = None

#     for part in organ_parts:
#         seg_path=os.path.join(subject_path, "segmentations",part)
#         data=nib.load(seg_path).get_fdata()
#         if not os.path.exists(seg_path):
#             print(f"{seg_path} n'existe pas, on saute")
#             continue
#         if not np.any(data > 0):
#             # return None or raise a controlled exception
#             raise ValueError("Empty segmentation (no organ present)")
#         mask =data.astype(np.uint8)
#         if mask_total is None:
#             mask_total = np.zeros_like(mask)
#         mask_total = np.logical_or(mask_total, mask)
#     return mask_total.astype(np.uint8)


def Create_organ_mask(case_id, organ_parts):
    """
    Crée un masque 3D pour un organe à partir de ses parties,  
    en ne déclenchant une erreur QUE si TOUTES les parties sont vides.

    INPUT:
     - case_id     : ID du cas
     - organ_parts : liste des fichiers de segmentation (.nii.gz)
     - z_min, z_max: indices pour recadrer en Z

    OUTPUT:
     - mask_total  : masque 3D (uint8) fusionnant toutes les parties non-vides

    Lève ValueError si aucun voxel trouvé pour TOUTES les parties.
    """
    base_dir    = "Data_set"
    subject_dir = os.path.join(base_dir, case_id, "segmentations")

    mask_total     = None
    any_part_found = False

    for part in organ_parts:
        seg_path = os.path.join(subject_dir, part)
        if not os.path.exists(seg_path):
            # fichier absent, on ignore
            continue

        seg_vol = nib.load(seg_path).get_fdata()
        # recadrage Z
        seg_crop = seg_vol

        if np.any(seg_crop > 0):
            any_part_found = True
            part_mask = (seg_crop > 0).astype(np.uint8)
            if mask_total is None:
                # initialisation à la forme de la première partie valide
                mask_total = np.zeros_like(part_mask, dtype=np.uint8)
            # fusion OR
            mask_total = np.logical_or(mask_total, part_mask)

    if not any_part_found:
        # aucune partie n'a de voxel : on déclenche l'erreur
        raise ValueError(f"Empty segmentation for all parts of organ ({organ_parts}) in case {case_id}")

    return mask_total.astype(np.uint8)



def Create_skin_mask_bis(case_id,organ_parts):
    """
    INPUT:
    - case_id: ID du cas,
    - organ_parts: liste des parties de l'organe à considérer,
    OUTPUT:
    - rest_of_body_clean: masque du reste du corps,(en 3D)
    - outside_mask_clean: masque de l'extérieur,(en 3D)
    - mask_total: masque de l'organe, (en 3D)
    -----La différence avec Create_skin_mask est que le masqque de la peau créé ici contient aussi les organes non inclus dans organ_parts,
    """


    base_dir=r'Data_set'   
    ct_path=os.path.join(base_dir,case_id,"ct.nii.gz")
    subject_path = os.path.join(base_dir, case_id)
    

    mask_total = Create_organ_mask(case_id,organ_parts)



    mask_total=mask_total.astype(np.uint8)

    ct=nib.load(ct_path)
    ct_data=ct.get_fdata()

    affine=nib.load(ct_path).affine
    
    body_mask=(ct_data>-300)# we condider that the outside has small density, so we can use a threshold
                            # could be adjusted based on the subject                        
    body_mask=keep_largest_component(body_mask).astype(np.uint8)

    body_mask_filled = binary_fill_holes(body_mask).astype(np.uint8)  #some zones inside the body may not be filled,in case we only work one or few organs
    rest_of_body_mask_filled=(body_mask_filled & (~mask_total)).astype(np.uint8)
    radius =1
    size = 2 * radius + 1  # 9

    # créer les grilles de coordonnées centrées en 0
    # np.ogrid est plus mémoire-efficace que meshgrid pour ce cas
    Z, Y, X = np.ogrid[-radius:radius+1, -radius:radius+1, -radius:radius+1]

    # masque sphérique : inclusion si distance <= radius
    structure = (X**2 + Y**2 + Z**2) <= radius**2

    rest_of_body_mask_filled = binary_closing(rest_of_body_mask_filled, structure=structure).astype(np.uint8)  #Ferme les trous
    rest_of_body_mask_filled=binary_fill_holes(rest_of_body_mask_filled).astype(np.uint8)       #Same
    outside_mask_filled=(~np.logical_or(rest_of_body_mask_filled,mask_total)).astype(np.uint8)



    outside_mask_clean = keep_largest_component(outside_mask_filled)
    outside_mask_clean =binary_opening(outside_mask_clean, structure=structure).astype(np.uint8)  #Ferme les trous
    rest_of_body_clean = (~np.logical_or(outside_mask_clean,mask_total)).astype(np.uint8)

    
    
    # body=np.logical_or(rest_of_body_clean,mask_total).astype(np.uint8)
    return rest_of_body_clean,outside_mask_clean,mask_total

    





def Create_mask_2D(masks,z):
    """
    INPUT:
    - masks: liste de masques 3D à combiner, chaque masque doit être de la même taille,
    - masks[i]: masque i, de forme (H,W,D), où H est la hauteur, W la largeur et D la profondeur,
    - z: indice de la tranche à extraire (par exemple, 340 pour la tranche 340).
    OUTPUT:
    - mask_2D: masque 2D de la coupe en z,pour valeur value là où est le masque, 0 ailleurs.
    """
    mask_eit = np.zeros_like(masks[0][:,:,z])
    for i,mask in zip(range(0,len(masks)),masks):
        mask_eit[mask[:,:,z]>0]=i
    # mask_closed = binary_closing(mask_eit, structure=disk(3))
    # mask_filled = binary_fill_holes(mask_closed)
    return mask_eit.astype(np.uint8)


def show_masks_eit(mask_eit):
    """""
    INPUT:
    - mask_eit: la coupe en z de la CT avec  masques, pour EIT,en 2D
    OUTPUT:
    - Affiche la coupe en z de la CT avec les masques pour EIT.
    """""
    plt.figure(figsize=(6,6))
    plt.imshow(mask_eit, interpolation='nearest')
    plt.title('Tranche EIT (mask_eit)')
    plt.axis('off')
    
    plt.show()


"this function only works with napari installed, which you probably struggle to do (see README.md file), so I commented it out"
def show_masks_3D(organ_mask,body_mask,outside_mask,ct_data, case_id):
    viewer = napari.Viewer()
    viewer.add_image(ct_data, name="CT of "+case_id, colormap="gray", contrast_limits=[-200, 500])

    viewer.add_labels(organ_mask, name='lungs mask', opacity=0.5)
    viewer.add_labels(body_mask, name='body mask Mask', opacity=0.5)

    viewer.add_labels(outside_mask, name='outside Mask', opacity=0.5)


    napari.run()

from scipy.ndimage import binary_closing, binary_fill_holes
from skimage import measure
from shapely.geometry import Polygon
from skimage.morphology import disk

def Extract_contour(
    mask2d: np.ndarray,
    closing_radius: int = 5,
    opening_radius: int = 3,
    simplify_tol: float = 0.01,
):
    """
    Extrait le contour extérieur du 'body total' (mask2d==1) uniquement.

    - mask2d : 2D array int
        0 = fond, 1 = corps total, puis >=2 = organes internes
    - closing_radius : rayon en px pour enlever les troux
    - opening_radius : rayon en px pour lisser les petits pics

    Retourne un shapely.Polygon CCW en en coord normalisées.
    """
    # 1) Isolation du corps total
    body = (mask2d == 1)

    # 2) Combler les trous éventuels
    body = binary_fill_holes(body)

    # 3) Closing puis Opening pour lisser
    struct_c = disk(closing_radius)
    body = binary_closing(body, structure=struct_c)
    struct_o = disk(opening_radius)
    body = binary_opening(body, structure=struct_o)

    # 4) Extraction des contours
    contours = measure.find_contours(body.astype(float), level=0.5)
    if not contours:
        raise RuntimeError("Aucun contour trouvé pour mask2d==1 !")
    cnt = max(contours, key=lambda c: c.shape[0])

    # 5) Pixel→coordonnées normalisées [-1,1]
    h, w = mask2d.shape
    xs = (cnt[:,1] / (w-1)) * 2 - 1
    ys = 1 - (cnt[:,0] / (h-1)) * 2

    poly = Polygon(np.vstack([xs, ys]).T)

    # 6) Simplification & orientation
    poly = poly.simplify(simplify_tol, preserve_topology=True)
    if not poly.exterior.is_ccw:
        poly = Polygon(poly.exterior.coords[::-1])

    return poly