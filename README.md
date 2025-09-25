# Project_EIT_CT_2025 — CT → EIT Forward & Inverse

> Toolbox MATLAB complète pour transformer des **coupes CT** et des **masques d’organes** en **maillages FEM** (PDE Toolbox), simuler le **problème direct EIT** avec OOEIT, et résoudre le **problème inverse EIT** par différentes stratégies d’initialisation (morphologique ou homogène).
> Comprend extraction de contours, génération de maillages, placement d’électrodes, simulation forward, reconstruction inverse et plotting automatisé (batch).

---

## 📁 Structure du dépôt

```
OOEIT-main/                          # (optionnel) backend EIT FEM (OOEIT)
PROJECT_EIT_CT_2025-MAIN/
├─ BatchLogs/                        # logs & snapshots d’itérations
├─ Data_set/                         # CT/segmentations par patient
├─ meta/                             # configs (conductivités, mappings, etc.)
├─ src_anis/                         # (ce dossier) tout le pipeline MATLAB
├─ Outputs/                          # sorties par patient/slice
├─ Batch.m                           # exécution batch forward + plots
├─ Plot_Problem_Forward_fn.m         # helpers plot forward
├─ Plot_Problem_Forward.m            # démo/driver forward
├─ Problem_Inverse.m                 # driver inverse (init morphologique)
└─ Plot_Problem_Inverse.m            # plot inverse (init morphologique)
```

*(il existe aussi des variantes « ConstInit » pour l’inverse avec init homogène)*

---

## 🔧 Prérequis

* **MATLAB** R2023b+ recommandé, avec **PDE Toolbox**.
* **OOEIT** (ou équivalent) si tu utilises le backend `EITFEM`, `ForwardMesh1st`, etc.

  > Place-le dans `OOEIT-main/` ou ajoute son chemin au path MATLAB.
* Volumes **NIfTI** (`.nii`/`.nii.gz`) pour CT & segmentations 3D.
* OS testés : macOS (Apple Silicon), Linux, Windows.

**Installation rapide (MATLAB)**

```matlab
addpath(genpath('PROJECT_EIT_CT_2025-MAIN/src_anis'));
addpath('PROJECT_EIT_CT_2025-MAIN');
% (si OOEIT)
addpath(genpath('OOEIT-main'));
savepath;
```

---

## 📦 Données attendues

```
PROJECT_EIT_CT_2025-MAIN/Data_set/
└─ s0011/
   ├─ ct.nii.gz
   └─ segmentations/
      ├─ heart.nii.gz
      ├─ lung_left.nii.gz
      ├─ lung_right.nii.gz
      ├─ aorta.nii.gz
      └─ ... (masques binaires par organe)
```

Le mapping **fichier segmentation → groupe** est géré dans `extract_slice_polys.m` (tables d’équivalences, ex. “aorta” ∈ Blood).
Le CT et chaque masque 3D doivent couvrir la **même grille**.

---

## 🚶‍♀️ Workflow (de bout en bout)

### 1. Extraction contours & organes (slice z)

```matlab
[domain, shapes, ext_smooth, info] = extract_slice_polys(patient, z);
```

* `ext_smooth` : contour externe lissé
* `shapes` : polyshape par groupe d’organes
* `domain` : soft tissue

### 2. Visualiser contours

```matlab
plot_slice_segmentee(...);
plot_contours_avant_mesh(...);
```

### 3. Générer le maillage

```matlab
mb = MeshBuilder.from_polys(shapes, domain, ext_smooth, z, params);
mb.saveMesh('Outputs/.../mesh', z);
```

* `g, H, triGroup` : données FEM

### 4. Placer les électrodes

```matlab
[E, el_centers] = buildElectrodesFromContour(g, ext_smooth, Ne, el_width, H);
```

### 5. Forward EIT

```matlab
sim = EITSim(g, H, triGroup, domain, params, E, el_centers);
sim.simulate_forward_full();
sim.save_results('Outputs/...'); % écrit eit_pack.mat
```

* `Vmat` : Ne×K
* `Imeas` : Ne*K×1

### 6. Inverse EIT (slice t cible, init morpho depuis slice t−1)

```matlab
rec = simulate_inverse_slice(patient_t, z_t, patient_prev, z_prev, opts, iterDir);
```

* `sigma_rec` nodale
* `sigma_rec_tri` par triangle
* `misfit`
* snapshots itérations dans `iterDir/`

---

## ⚡️ Quick start complet

```matlab
% Batch complet forward + plots
Batch('RequirePlots', true);

% Exemple patient/slices
Batch('OnlyPatients', 's0011', 'OnlyZSlices', 280:300);

% Inverse morphologique
Problem_Inverse; 
Plot_Problem_Inverse;

% Inverse init homogène
Problem_Inverse_ConstInit;
Plot_Problem_Inverse_ConstInit;
```

---

## 🧠 API principale (src_anis)

* `EITSim` — simulation EIT (forward + Jacobien)
* `buildElectrodesFromContour` — placement d’électrodes
* `MeshBuilder` — maillage PDE
* `extract_slice_polys` — extraction polyshapes d’organes
* `simulate_inverse_slice` — pipeline inverse complet
* `init_from_prev_slice_discrete` — init nodale discrète depuis slice précédente
* `plot_slice_segmentee`, `plot_contours_avant_mesh` — visualisations

---

## ⚙️ Formats & conventions

* **Unités** : mm côté pré-processing/mesh ; conversion mm→m automatique dans `EITSim`.
* **Groupes d’organes** typiques : `{'SoftTissue','Heart','Lung','Trachea','Bone','Cartilage','Blood','LiverSpleenPancreas','Kidney','Other'}`.
* **Conductivités** : S/m (SI), via `params.cond.<Group>`.

---

## 🧪 Reproductibilité

* Chaque run `Batch` écrit dans `BatchLogs/` :

  * config complète (JSON/MAT)
  * hash git de l’état du repo
  * timings & métriques
  * seed aléatoire
* Snapshots d’itération inverse dans `Outputs/.../iters/`.

---

## 🐛 Dépannage rapide

* **Contours incohérents / mesh fail** → ajuster `opts.smooth_K/step`, `min_area_mm2`, orientation (radiological vs neurological).
* **Électrodes mal réparties** → vérifier `el_width` (mm, arc length) vs périmètre, et `Ne`.
* **Forward diverge / voltages bizarres** → vérifier unité (mm→m), `z_contact`, patterns (somme nulle).
* **Inverse instable** → augmenter `noser_only.lambda`, réduire `max_rel_step`, régler `iters_per_stage`.

---


Code généré manuellement avec aide de la LLM d'Open AI ChatGpt5

