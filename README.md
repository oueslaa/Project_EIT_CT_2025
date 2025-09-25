# Project_EIT_CT_2025 — CT → EIT Forward & Inverse

> Toolbox MATLAB complète pour transformer des **coupes CT** et des **masques d’organes** en **maillages FEM** (PDE Toolbox), simuler le **problème direct EIT** avec OOEIT, et résoudre le **problème inverse EIT** par différentes stratégies d’initialisation (morphologique ou homogène).
> Comprend extraction de contours, génération de maillages, placement d’électrodes, simulation forward, reconstruction inverse et plotting automatisé (batch).

---

## 📥 Téléchargements

Pour tester rapidement la toolbox sans avoir à préparer tes propres données :

* **Data_set.zip** – petit dataset CT + segmentations d’exemple
* **Outputs.zip** – sorties complètes générées par `Batch.m` sur ce dataset (maillages, packs forward, PNG)

[⬇️ Télécharger Data_set.zip et Outputs.zip](https://filesender.renater.fr/?s=download&token=3ee6e5b5-e57e-4e11-82ef-5acff3fd2a3d)

Décompresse les archives dans la racine du projet, tu obtiendras :

```
PROJECT_EIT_CT_2025-MAIN/
├─ Data_set/       # issu de Data_set.zip
└─ Outputs/        # issu de Outputs.zip
```

Ensuite, lance directement les scripts de plotting pour voir les résultats sans recalculer.

---

## 📁 Structure du dépôt

```
OOEIT-main/                          # (optionnel) backend EIT FEM (OOEIT)
PROJECT_EIT_CT_2025-MAIN/
├─ BatchLogs/                        # logs & snapshots d’itérations
├─ Data_set/                         # CT/segmentations par patient (exemple fourni)
├─ meta/                             # configs (conductivités, mappings, etc.)
├─ src_anis/                         # (ce dossier) tout le pipeline MATLAB
├─ Outputs/                          # sorties par patient/slice (exemple fourni)
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


# 📂 Fichiers principaux : rôle, comment les lancer, et ce qu’ils génèrent

## Batch.m — exécution batch du **forward** + **plots** (idempotent & checkpointé)

**Ce que fait le script**

* Lit `meta/meta.csv`, garde les lignes `study_type == "ct neck-thorax-abdomen-pelvis"`, filtre éventuellement par patients/slices.
* Pour chaque `(patient_id, z_slice)` où il y a **cœur ET poumons** (détectés depuis `segmentations/*.nii*`) :

  * saute si le forward/plots existent déjà (**idempotent** via `Outputs/...` + `BatchLogs/completed.tsv`),
  * sinon lance `Problem_Forward_fn` (forward) puis `Plot_Problem_Forward_fn` (export PNG).
* Tient un **checkpoint** persistant `BatchLogs/completed.tsv` (lignes `patient_id \t z`) et un **log** daté dans `BatchLogs/*.log`.

**Comment le lancer**

```matlab
% Tout le dataset (plots obligatoires pour considérer “done”)
Batch('RequirePlots', true);

% Restreindre à un patient et à un sous-ensemble de z (debug)
Batch('OnlyPatients','s0011','OnlyZSlices',280:300);

% Accepter "done" dès que le forward est présent, même si PNG manquent
Batch('RequirePlots', false);
```

**Paramètres utiles**

* `RequirePlots` (bool, def `false`) : si `true`, une slice est “done” seulement si les PNG clés existent.
* `OnlyPatients` (`char` ou `cellstr`) : liste d’IDs patients à traiter.
* `OnlyZSlices` (vector num) : indices `z` à limiter.

**Ce que ça génère**

* `BatchLogs/yyyymmdd_HHMMSS_BatchTorsoForward.log` (trace)
* `BatchLogs/completed.tsv` (checkpoint `patient\tz`)
* Pour chaque `(patient,z)`:

  ```
  Outputs/<patient>/slice_<ZZZ>/
  ├─ mesh/mesh_slice<ZZZ>.mat      % g, H, triGroup, shapes, domain, contour
  ├─ E_electrodes.mat              % E (cell d’arêtes), el_centers
  ├─ cond_values.mat               % struct cond
  ├─ viz_config.mat                % displayMode, sigma_clim, cmap
  ├─ Imeas.mat, U_cell.mat, Vmat.mat, Ipat.mat, Efield_tri.mat
  ├─ eit_pack.mat                  % pack forward complet (pour plot/inverse)
  ├─ data_sources.mat              % chemins CT/seg + (patient,z)
  └─ plots_forward/                % PNG haute def (voir plus bas)
  ```

---

## Problem_Forward_fn.m — **forward** EIT (fonction)

**Ce que fait le script**

* Prépare maillage 2D (via `MeshBuilder.from_polys`), place `Ne` électrodes d’arc `ElWidthMM` sur le contour lissé, fixe les **conductivités** par groupe, et appelle `EITSim(...).simulate_forward_full()`.
* Sauvegarde toutes les **sorties forward** (+ config de viz) dans `Outputs/<patient>/slice_<ZZZ>/`.

**Lancer (exemple)**

```matlab
Problem_Forward_fn('s0011', 320, ...
  'Ne', 16, 'ElWidthMM', 33, 'Iamp', 0.35e-3, ...
  'Zcontact', 0.05, 'NoiseRel', 1e-3, 'DisplayMode','neurological');
```

**Sorties**

* Idem qu’en batch (voir arborescence ci-dessus).

---

## Plot_Problem_Forward_fn.m — **export** de toutes les figures forward (sans affichage)

**Ce que fait le script**

* Charge `mesh_slice*.mat`, `eit_pack.mat`, `viz_config.mat`, `data_sources.mat`, `E_electrodes.mat`.
* Produit **tous** les PNG (sans ouvrir de fenêtres par défaut).

**Lancer**

```matlab
% Exporte en PNG (300 DPI), n’affiche pas
Plot_Problem_Forward_fn('s0011', 320, 'SavePNG', true, 'DPI', 300, 'ShowFigures', false);

% Pour voir les figures à l’écran
Plot_Problem_Forward_fn('s0011', 320, 'ShowFigures', true);
```

**PNG générés** (dans `plots_forward/`, clés pour `RequirePlots=true`)

* `ct_brut.png`
* `slice_segmentee.png`
* `contours_avant_mesh.png`
* `mesh_only.png`
* `mesh_electrodes_numbered.png`
* `mesh_electrodes_sigma_forward.png`
* `electrode_voltages_heatmap.png`, `electrode_voltages_profile.png`
* `potentials_montage_4x4.png`
* `Efield_mag_tri.png`, `Efield_quiver.png`
* `currents_heatmap.png`, `currents_profile.png`
* `reciprocity_heatmap.png`

---

## Problem_Forward.m — **exemple script** forward (batch-only, pas d’affichage)

**Ce que fait le script**

* Version “driver” monofile du forward (équivalent à `Problem_Forward_fn`) pour un `(patient,z)` codé en dur.
* N’affiche **aucune** figure (destiné au batch/démo).

**Lancer**

```matlab
edit Problem_Forward.m   % adapte patient/z et chemins OOEIT si besoin
Problem_Forward
```

**Sorties**

* Identiques à `Problem_Forward_fn` (voir arborescence).

---

## Plot_Problem_Forward.m — **visualisation** complète du forward (avec affichage)

**Ce que fait le script**

* Charge les mêmes fichiers que la version `_fn` mais **affiche** les figures (et exporte en PNG).
* Ajoute des plots diagnostiques (histogrammes aires/angles, histogramme σ, SNR, profils par motif).

**Lancer**

```matlab
edit Plot_Problem_Forward.m   % adapte patient/z si besoin
Plot_Problem_Forward
```

**Sorties**

* `plots_forward/` (mêmes PNG qu’au-dessus + diagnostics : `mesh_area_hist.png`, `mesh_minangle_hist.png`, `sigma_hist.png`, `noise_hist.png`, `electrode_profile_motifXX.png`…)

---

## Problem_Inverse.m — **reconstruction** (init **morphologique** depuis un slice “donneur”)

**Ce que fait le script**

* Charge le **pack forward** de la cible `(patient_this, slice_this)`.
* Récupère les **organes** d’un `(patient_prev, slice_prev)` et les **transplante** (warp TPS) dans le contour cible pour construire une **init discrète structurée** (pas homogène).
* Construit un **maillage inverse** (contour seul, paramétrisation nodale), définit le planning NOSER+GN et lance `simulate_inverse_slice`.
* Sauvegarde la reco + **snapshots** d’itérations.

**Lancer**

```matlab
edit Problem_Inverse.m  % adapte patients/slices et hyperparams
Problem_Inverse
```

**Sorties**

```
Outputs/<patient_this>/slice_<ZZZ>/inverse/prev_<patient_prev>_<ZZZprev>/
├─ reconstruction_inverse.mat   % sigma0, sigma_rec(_tri), misfit, Umeas, Ufwd, scaleU, initInfo, g,H...
└─ plots_inverse/
   └─ iterations/               % iter_*.mat (snapshots pour plots)
```

---

## Plot_Problem_Inverse.m — **plots** de la reco (init morpho)

**Ce que fait le script**

* Affiche **Référence** (si disponible), **Init** (discrète), **Reco**, **intermédiaires** regroupés par stage, **signaux collés** (meas/init/reco + intermédiaires), **erreurs**, et **trajectoires** (objective/L2) si dispo.
* Utilise l’échelle couleur de `viz_config.mat` si possible.

**Lancer**

```matlab
edit Plot_Problem_Inverse.m
Plot_Problem_Inverse
```

**PNG générés** (dans `plots_inverse/`)

* `sigma_ref_*.png` (si GT)
* `sigma_init_*.png`
* `sigma_reco_*.png`
* `montage_ref_init_reco_*.png`
* `intermediate_recos_stage*.png`
* `signals_*.png`, `signal_errors_*.png`
* `l2_datafit_*.png`, `objective_vs_iter_*.png` (si traqueur)

---

## Problem_Inverse_ConstInit.m — **reconstruction** avec **init homogène** (σ constant)

**Ce que fait le script**

* Ignore les organes “donneurs”, **désactive** la transplantation et **impose** une init homogène `sigma_0 = INIT_CONST`.
* Calcule un clip dynamique depuis les conds du pack forward.
* Lance NOSER(+GN) et sauve tout comme ci-dessus, mais dans un dossier **spécifique à la valeur** d’init.

**Lancer**

```matlab
edit Problem_Inverse_ConstInit.m   % adapte patient/z, hyperparams
Problem_Inverse_ConstInit
```

**Sorties**

```
Outputs/<patient>/slice_<ZZZ>/inverse/initConst_<valFormatee>/
├─ reconstruction_inverse.mat
└─ plots_inverse/ (+ iterations/)
```

---

## Plot_Problem_Inverse_ConstInit.m — **plots** de la reco (init homogène)

**Ce que fait le script**

* Retrouve automatiquement le **dernier run initConst_*** si la valeur initiale ne correspond pas.
* Affiche Référence (si dispo), **Init discrète** (constante), Reco, Intermédiaires, Signaux, Erreurs, Trajectoires — mêmes figures que la version “morpho” mais taguées `initConst`.

**Lancer**

```matlab
edit Plot_Problem_Inverse_ConstInit.m
Plot_Problem_Inverse_ConstInit
```

**PNG générés**

* `sigma_ref_sliceXXX_initConst_<val>.png`
* `sigma_init_sliceXXX_initConst_<val>.png`
* `sigma_reco_sliceXXX_initConst_<val>.png`
* `montage_ref_init_reco_sliceXXX_initConst_<val>.png`
* `intermediate_recos_stage*_sliceXXX_initConst_<val>.png`
* `signals*_sliceXXX_initConst_<val>.png`, `signal_errors*_sliceXXX_initConst_<val>.png`
* `l2_datafit_sliceXXX_initConst_<val>.png` (si dispo)

---

## Récap : quoi lancer pour quoi obtenir ?

| Besoin                                       | Script à lancer                                    | Affichage          | Dossier/ fichiers produits                                                       |
| -------------------------------------------- | -------------------------------------------------- | ------------------ | -------------------------------------------------------------------------------- |
| Forward **à grande échelle** (auto-skip)     | `Batch(...)`                                       | Non (PNG exportés) | `Outputs/<patient>/slice_<ZZZ>/...` + `plots_forward/` + `BatchLogs/*`           |
| Forward **unitaire** (fonction)              | `Problem_Forward_fn(patient,z,...)`                | Non                | `Outputs/<patient>/slice_<ZZZ>/...`                                              |
| Plots forward **unitaires** (export/affiche) | `Plot_Problem_Forward_fn` / `Plot_Problem_Forward` | Oui/Non            | `plots_forward/*.png`                                                            |
| Inverse **morphologique**                    | `Problem_Inverse`                                  | Non                | `inverse/prev_.../reconstruction_inverse.mat` + `plots_inverse/iterations/*.mat` |
| Plots inverse **morphologique**              | `Plot_Problem_Inverse`                             | Oui                | `plots_inverse/*.png`                                                            |
| Inverse **init homogène**                    | `Problem_Inverse_ConstInit`                        | Non                | `inverse/initConst_*/reconstruction_inverse.mat`                                 |
| Plots inverse **init homogène**              | `Plot_Problem_Inverse_ConstInit`                   | Oui                | `plots_inverse/*.png`                                                            |

---

### Notes pratiques

* **Chemins OOEIT** : vérifie la ligne `addpath(genpath('/Users/anis/.../OOEIT-main'))` dans les scripts `Problem_*`/`Plot_*` et adapte-la à ta machine (ou place OOEIT à la racine du repo et utilise la détection déjà présente dans `*_fn.m`).
* **DisplayMode** : `'neurological'` (X normal) vs `'radiological'` (X inversé, convention DICOM).
* **Idempotence batch** : pour forcer la régénération d’une slice, **supprime** son dossier `Outputs/<patient>/slice_<ZZZ>/plots_forward` (si `RequirePlots=true`) et/ou les fichiers `mesh/eit_pack/E_electrodes` selon le besoin — ou **retire** la ligne correspondante dans `BatchLogs/completed.tsv`.



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

