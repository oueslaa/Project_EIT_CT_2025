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
from scipy.sparse import csr_matrix, bmat,hstack
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
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from Extract_skin_mask import *
from help_functions import *

from pyeit.mesh.utils import edge_list
from pyeit.mesh import PyEITMesh






class CEMFEM:
    def __init__(
        self,
        mesh: PyEITMesh,
        protocol: PyEITProtocol,
        ref_node: int = 0,
    ):
        """
        CEM‐FEM forward solver with integrated PyEITProtocol.

        mesh     : PyEITMesh
        protocol : PyEITProtocol (ex_mat, meas_mat, keep_ba)
        ref_node :index of reference electrode(Dirichlet comme dans pyeit,)
        """
        self.mesh     = mesh
        self.protocol = protocol
        self.N        = mesh.n_nodes
        self.L        = mesh.n_el
        self.ref      = ref_node

        # calcul des stiffness locaux(ce terme est indéendant des electrodes)
        self.ke = calculate_ke(mesh.node, mesh.element)
        # Rrécuper  des arêtes de bord du domaine
        self.boundary_edges = edge_list(mesh.element)
        # les regrouper par éléctrode
        self.electrode_edges = [
            [e for e in self.boundary_edges if mesh.el_pos[l] in e]
            for l in range(self.L)
        ]
        

    def _edge_integrals_1D(self, edge):
        """∫_edge φiφj dS et ∫_edge φi dS on a segment
           OUTPUT: - first and secend integral
        """
        i, j = edge
        xi, yi = self.mesh.node[i,0], self.mesh.node[i,1]
        xj, yj = self.mesh.node[j,0], self.mesh.node[j,1]
        Le = np.hypot(xi - xj, yi - yj)
        intS = {(i,i):Le/3, (j,j):Le/3, (i,j):Le/6, (j,i):Le/6}
        intM = {i:Le/2, j:Le/2}
        return Le, intS, intM
    
    def assemble_CEM(self, perm, zeta, mode: str,ex_line):
        """
        Assemble the CEM system A u = b for ONE excitation (self.protocol.ex_mat[0]),
        in either 'current' or 'potential' mode.

        Parameters
        ----------
        perm : array_like
        Conductivity values per element.
        zeta : array_like of length L
        Contact impedances.
        mode : str, either "current" or "potential"
        - "current":  inject I, solve for U_electrodes
        - "potential": inject U, solve for I_electrodes

        Returns
        -------
        A : csr_matrix, shape ((N+L-1),(N+L-1))
        Complete electrode‐FEM matrix.
        b : ndarray, shape (N+L-1,)
        Right‐hand side for the chosen mode.

        Notes : -- Reference for mode "current" : https://www.researchgate.net/publication/12807310_Three-dimensional_electrical_impedance_tomography_based_on_the_complete_electrode_model
                -- For mode "potential" I just implemented the same code as OOEIT
        """
        # 0) Sizes
        N = self.mesh.n_nodes    # nombre de nœuds
        L = self.mesh.n_el       # nombre d'électrodes
        Lr = L - 1               # hors référence

        #  Global stiffness B (indepndant des electrodes)
        B = assemble(
            self.ke,
            self.mesh.element,
            self.mesh.get_valid_perm_array(perm),
            N, ref=-1
        )
        
        # Surface‐integrals pour M et S_sum
        M     = np.zeros((L, N), dtype=float)
        S_sum = lil_matrix((N, N), dtype=float)
        intB  = np.zeros(L, dtype=float)

        for l, edges in enumerate(self.electrode_edges):
            # on parcourt TOUTES les arêtes de l'électrode l
            for edge in edges:
                i, j = map(int, edge)                   # deux nœuds du segment
                Le, intS_e, intM_e = self._edge_integrals_1D((i, j))
                intB[l] += Le                           # somme des longueurs/aires

                # contribution à M pour CHAQUE segment
                for node_idx, val in intM_e.items():
                    M[l, node_idx] += val

                # contribution à S_sum si besoin
                for (ni, nj), val in intS_e.items():
                    S_sum[ni, nj] += val / zeta[l]

        # Blocs A11, A12 and D
        A11 = (B + S_sum.tocsr()).tocsr()
      
        # A12[:, j] = –(1/z₁ M[0,:] – 1/z_{j+1} M[j+1,:])
        A12_arr = np.zeros((N, Lr), dtype=float)
        for j in range(Lr):
            A12_arr[:, j] = -(
                (1.0 / zeta[0]) * M[0, :]
                - (1.0 / zeta[j+1]) * M[j+1, :]
            )
        A12 = csr_matrix(A12_arr)

        # D[i,i] = |e₁|/z₁ + |e_{i+1}|/z_{i+1},  D[i≠j] = |e₁|/z₁
        n_ref = intB[0]
        Dmat = np.full((Lr, Lr), fill_value=(n_ref / zeta[0]), dtype=float)
        for i in range(Lr):
            
            Dmat[i, i] = (n_ref / zeta[0]
                        + intB[i+1] / zeta[i+1])
        D = csr_matrix(Dmat)

        # final CEM matrix
        A = bmat([
            [A11,     A12  ],
            [A12.T,   D    ]
        ], format="csr")
        

        ### construction of the b vector
        a, b_idx = ex_line
        if mode == "current":
            # Inject current  +1/–1 → get U
            I = np.zeros(L, float)
            I[a] = +1.0
            I[b_idx] = -1.0

            b_int = np.zeros(N, dtype=float)

            #    b_el[j] = I[0] - I[j+1], for j=0..L-2 (according to the ref cited in the notes up)
            b_el = np.array([I[0] - I[j+1] for j in range(L-1)], dtype=float)
            

           
            b = np.concatenate([b_int, b_el])
            
        elif mode == "potential":
             
            A_12=csr_matrix((self.N, Lr), dtype=float)
      
            C_sparse = vstack([csr_matrix(np.ones((1, self.L-1))),-eye(self.L-1, format='csr')], format='csr')
            CtC_sparse = C_sparse.T.dot(C_sparse)
            D_inv_sp = diags(1.0/zeta)     
            intM_sp     = csr_matrix(M)            
            intM_scaled = D_inv_sp.dot(intM_sp)     

            S21_sp = - C_sparse.T.dot(intM_scaled)   # sparse (L-1 × N)

            # Assembler la matrice complète
            A = bmat([
                [A11,     A12  ],
                [S21_sp,   CtC_sparse    ]
            ], format="csr")

            #  Construire le vecteur b
            #    – on impose U[a]=+1, U[b_idx]=0
            U = np.zeros(L, dtype=float)
            U[a] = +1.0              #excitation with 1/0 V values like OOEIT
            U[b_idx] = 0

            L, N = M.shape
            n_r = L - 1

            
            intM = M.T    

            
          
            #    b_dom[i] = sumₗ intM[i,ℓ] * (U[ℓ]/zeta[ℓ])
            b_dom = intM.dot((1.0/zeta) * U)   # shape (N,)

         
            tmp  = (intB / zeta) * U           # shape (L,)
            b_el = - C_sparse.T.dot(tmp)              # shape (L-1,)

            
            b = np.concatenate([b_dom, b_el])  # shape (N + L - 1,)



            

        else:
            raise ValueError("mode must be 'current' or 'potential'")
        
        # print("le vecteur b est")
        # print(b)
        # print("la taille de b est ")
        # print(b.shape)
        # print("le nbre de noueds est ",self.N)
        # print("les eectrodes sont")
        # print(listee)

        return A, b
    
    def solve_CEM_ex_line(self, perm, zeta, mode:str,ex_line):
        """
          Solves for one excitation line
          returns: Eletrodes current/potential
        """
        A,b=self.assemble_CEM(perm, zeta, mode=mode,ex_line=ex_line)
        sol = spsolve(A, b)
        if mode=="current":
            return  np.concatenate(([-(sol[self.N:].sum())],sol[self.N:]))
        return  np.concatenate(([-sol[self.N:].sum()],sol[self.N:]))                    ##   np.concatenate(([sol[self.N:].sum()],sol[self.N:]))
    
    def solve_CEM(self,perm,zeta,mode:str):
        """
          Solves for the excitation defined in self.protocol
          returns: Eletrodes current/potential
          """
        ex_mat = self.protocol.ex_mat
        n_exc = ex_mat.shape[0]
        Lr = self.L - 1

        # pre-allocate result array
        Uall = np.array([])  
        # loop over each excitation line
        for k, ex_line in enumerate(ex_mat):
            # solve for this single excitation
            # print(ex_line)
            sol_el = self.solve_CEM_ex_line(perm, zeta, mode, ex_line)

            # sol_el is of length L-1
            Uall=np.concatenate((Uall,sol_el))  

        return Uall
    
        


def scale_pyeit_mesh_to_mm(mesh, case_id, slice_index, mask2d):
    """
    Prend un maillage PyEIT 2D (dont les noeuds sont en [-1,1],car pyeit renvoie un mesh normalisé sur [-1,1]),
    et le remet dans le repère réel (mm) de ton image CT.
    
    - mesh.node : (Nnodes,2) en coordonnées normalisées [-1,1]
    - case_id    : pour charger ct.nii.gz et son affine
    - slice_index: indice z de la coupe
    - mask2d     : le masque 2D original, shape (h, w)
    """
    # Récupère l'affine du CT
    ct = nib.load(f"Data_set/{case_id}/ct.nii.gz")
    affine = ct.affine  # 4×4
    
    h, w = mask2d.shape

    #  On renverse la normalisation [-1,1] → pixel indices [0..w-1], [0..h-1]
    #    x_norm = mesh.node[:,0] ∈ [-1,1] correspond à col
    #    y_norm = mesh.node[:,1] ∈ [-1,1] correspond à row (inversion d'axe Y)
    xn = mesh.node[:, 0]
    yn = mesh.node[:, 1]
    # colonnes pixels
    px = (xn + 1) * 0.5 * (w - 1)
    # lignes pixels (on a fait y = 1 - row/(h-1)*2 dans Extract_contour)
    py = (1 - yn) * 0.5 * (h - 1)

    #  Passe en coordonnées homogènes (i, j, z, 1) puis affine → (x_mm, y_mm, z_mm)
    ones = np.ones_like(px)
    zs   = np.full_like(px, slice_index)
    vox4 = np.vstack([px, py, zs, ones]).T         # (Nnodes,4)
    xyz  = vox4.dot(affine.T)[:, :3]               # (Nnodes,3)

    #  Remplace node par (x_mm, y_mm)
    mesh.node[:, 0] = xyz[:, 0]
    mesh.node[:, 1] = xyz[:, 1]

    return mesh



# def view_EIT_BP(mask_eit,mesh_obj,perm0,perm):
#         """ 
#         Visualize the EIT BP results.
#         :param mesh_obj: Mesh object created from the mask
#         :param ds_bp: BP reconstruction result
#         :param perm0: Initial conductivity distribution
#         :param v0: Initial potential
#         :param v1: Potential after applying the BP method
#         :param eit_bp: EIT BP object
#         :param protocol_obj: Protocol object used for EIT
#         :return: None
#         """

#         # delta_perm = perm - perm0
#         # print(delta_perm)
#         # pts      = mesh_obj.node      
#         # tri      = mesh_obj.element  
#         # tri_centers = mesh_obj.elem_centers 

#         # mesh_obj.perm = delta_perm  # Update the mesh object with the reconstructed conductivity

      

#         # Définis l'étendue du domaine en coordonnées réelles (en mm)
#         # extent = [xmin, xmax, ymin, ymax]

        

#         # # Affiche avec l'extent pour avoir les axes en mm
#         # im = axes[0].imshow(
#         #     mask_eit,
#         #     cmap='viridis',
#         #     extent=extent,
#         #     origin='lower',   # pour que (0,0) soit en bas à gauche
#         # )
#         # axes[0].set_title("Masque (0=fond,1=jaune,2=vert)")
#         # axes[0].set_xlabel("X (mm)")
#         # axes[0].set_ylabel("Y (mm)")
#         # axes[0].set_aspect('equal')#
#         fig,ax=plt.subplots(figsize=(6,6))

#         xmin, xmax = mesh_obj.node[:,0].min(), mesh_obj.node[:,0].max()
#         ymin, ymax = mesh_obj.node[:,1].min(), mesh_obj.node[:,1].max()

#         # Choisis ton pas en mm
#         step = 10.0
#         # Gère la gamme des ticks
#         ticks_x = np.arange(np.floor(xmin/step)*step, np.ceil(xmax/step)*step+step, step)
#         ticks_y = np.arange(np.floor(ymin/step)*step, np.ceil(ymax/step)*step+step, step)
  
#         # plots the mesh and the electrodes
#         create_mesh_plot(ax, mesh_obj,
#                         electrodes=mesh_obj.el_pos,
#                         coordinate_labels="radiological",
#                         marker_text_kwargs={ "color": "red", "fontsize": 6 })
#         ax.set_title("mesh + électrodes")
#         ax.axis("off")
       
     
#         ax.tick_params(axis='both', which='both', bottom=True, left=True, labelbottom=True, labelleft=True)
#         ax.set_xlabel("X (mm)")
#         ax.set_ylabel("Y (mm)")

#         # 4) Optionnel : ajouter des grilles à chaque 10 mm
#         ax.xaxis.set_major_locator(MultipleLocator(step))
#         ax.yaxis.set_major_locator(MultipleLocator(step))
#         ax.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.3)

#         plt.tight_layout()



#         # ax1.axis("equal")
#         # ax1.set_title(r"Input $\Delta$ Conductivities")
#         # im = ax1.tripcolor(pts[:, 0], pts[:, 1], tri,delta_perm, shading="flat")



#         # #plot the reconstructed BP conductivity
#         # ax2=axes[2]
#         # im = ax2.tripcolor(pts[:,0], pts[:,1], tri,ds,shading='flat', cmap='jet')
#         # ax2.set_title("BP reconstruit Δσ")
#         # ax2.axis("off")
#         # ax2.set_aspect("equal")
#         # fig.colorbar(im, ax=axes.ravel().tolist())
#         # # fig.savefig('../doc/images/demo_bp.png', dpi=96)
#         plt.show()


# def display_v0_v1(mesh_obj, protocol_obj, perm0, perm_true):
#     """
#     Calcule et affiche les matrices v0 et v1 (tensions aux électrodes)
#     pour perm0 (milieu homogène) et perm_true (milieu réel).

#     Paramètres
#     ----------
#     mesh_obj : PyEITMesh
#         Maillage PyEIT déjà construit.
#     protocol_obj : PyEITProtocol
#         Protocole EIT déjà créé (ex_mat, meas_mat, etc.).
#     perm0 : array_like, shape (n_elem,)
#         Vecteur de conductivité de référence (milieu homogène).
#     perm_true : array_like, shape (n_elem,)
#         Vecteur de conductivité “réel”.

#     Comportement
#     ------------
#     - Calcule v0 = fwd.solve_eit(perm=perm0) et v1 = fwd.solve_eit(perm=perm_true).
#     - Détermine n_exc = nombre d’électrodes, n_meas = len(v0) // n_exc.
#     - Remet en forme v0_mat, v1_mat de shape (n_exc, n_meas).
#     - Imprime v0_mat, v1_mat et leurs différences.
#     - Trace v0_mat[i], v1_mat[i], et Δv_mat[i] pour i=0 (pattern 0) en deux sous‐plots.

#     """
#     from pyeit.eit.fem import EITForward

#     # 1) Calcul des tensions v0 et v1
#     fwd = EITForward(mesh_obj, protocol_obj)

#     v0 = fwd.solve_eit(perm=perm0)
#     mesh_obj.perm = perm_true
#     v1 = fwd.solve_eit(perm=perm_true)

#     n_exc = mesh_obj.el_pos.size
#     total_meas = v0.shape[0]
#     n_meas = total_meas // n_exc

#     v0_mat = v0.reshape(n_exc, n_meas)
#     v1_mat = v1.reshape(n_exc, n_meas)
#     dv_mat = v1_mat - v0_mat

#     np.set_printoptions(precision=4, suppress=True)
#     # print(f"v0_mat shape = {v0_mat.shape}")
#     # print(v0_mat, "\n")
#     # print(f"v1_mat shape = {v1_mat.shape}")
#     # print(v1_mat, "\n")
#     # print("Δv_mat = v1_mat - v0_mat\n", dv_mat)

#     # 5) Tracé : pattern 0 (i=0)
#     i = 4   
#     fig, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True)

#     # (a) v0 vs v1 pour pattern 0
#     axes[0].plot(v0_mat[i, :], '-o', label="v0 (réf)")
#     axes[0].plot(v1_mat[i, :], '-s', label="v1 (réel)")
#     axes[0].set_title("Tensions v0 et v1 (pattern 0)")
#     axes[0].set_ylabel("Tension")
#     axes[0].legend(loc="upper right", fontsize="small")
#     axes[0].grid(True)

#     # (b) Δv pour pattern 0
#     axes[1].plot(dv_mat[i, :], '-^', label="Δv = v1-v0", color="red")
#     axes[1].set_title("Δv (pattern 0)")
#     axes[1].set_xlabel("Indice de mesure")
#     axes[1].set_ylabel("ΔTension")
#     axes[1].legend(loc="upper right", fontsize="small")
#     axes[1].grid(True)

#     plt.tight_layout()
#     plt.show()
