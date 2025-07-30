import os
import pandas as pd
import numpy as np
import nibabel as nib
from skimage import measure
import trimesh
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
from nibabel.affines import apply_affine



def make_mesh_MC(case_id,organ_name,level):

    #loads: case_id: subject, organ_name: organ name(with .nii.gz), level: iso-surface in wich we slice through.##Needs automatic level value juste like with alpha, this might be difficul i'll come back later.

    #Output :mesh: surface mesh to visualise the volume
    #        verts_mm: coordonnées des points formant les triangles
    #        faces_mm : les faces des triangles
    #        volume_cm3,surface_cm2: volume et surface
    base_dir=r'Data_set'
    subject_path = os.path.join(base_dir, case_id)
    path=os.path.join(subject_path, "segmentations", organ_name)
    

    img=nib.load(path)
    data=img.get_fdata().astype(np.float32)
    affine=img.affine  # matrix that passes from the coordinates base to actual base (in mm)

    if not np.any(data > 0):
        # return None or raise a controlled exception
        raise ValueError("Empty segmentation (no organ present)")

    verts,faces,_,_ =measure.marching_cubes(data,level=level) #verts are coordinates
    verts_mm = apply_affine(affine, verts)  
    
    mesh=trimesh.Trimesh(vertices=verts_mm,faces=faces)
    volume_cm3=abs(mesh.volume)/1000         
    surface_cm2=mesh.area/100

    
    return mesh,verts_mm,faces,volume_cm3,surface_cm2

def make_multi_part_mesh_MC(case_id,part_list,organe_name,level):
    #Same as previous, with multipart organ
    mshs=[]
    verts=[]
    faces=[]
    volume_cm3=0
    surface_cm2=0
    for part in part_list:
        try:
            ms,v,f,vo,su=make_mesh_MC(case_id,part,level)
            mshs.append(ms)
            verts.append(v)
            faces.append(f)
            volume_cm3+=vo
            surface_cm2+=su
        except ValueError as e:
            msg = str(e)
            if "Empty segmentation" in msg:
              print(f"{case_id} skipped: no organ present")
            else:
                print(f"Erreur pour {case_id}: {e}")
    
    return mshs,verts,faces,volume_cm3,surface_cm2

#plots mesh surface,whatever method used
def view_mesh_3D(verts_mm,faces,organ_name):


    x,y,z=verts_mm.T
    i,j,k=faces.T
    edges=[]
    for tri in faces:
        pts=[tri[0],tri[1],tri[2],tri[0]]    #a triagnle ABC is the three segments A->B, B->C,C->A
        for p1, p2 in zip(pts,pts[1:]):
            edges.append(((verts_mm[p1],verts_mm[p2]))) 


    #consruct the 3D representation: link each edges one to the next
    xe,ye,ze=[],[],[]
    for (p1, p2) in edges:
        xe+=[p1[0],p2[0],None]
        ye+=[p1[1],p2[1],None]
        ze+=[p1[2],p2[2],None]

    fig = go.Figure()

    #plot the 3D mesh
    fig.add_trace(go.Mesh3d(
        x=x,y=y,z=z,
        i=i,j=j,k=k,
        color='magenta',
        opacity=0.6,
        name=organ_name,
        lighting=dict(ambient=0.4,diffuse=0.8),
        lightposition=dict(x=100,y=100,z=100)))
    # opptionel : show the edges of each triangle
    #fig.add_trace(go.Scatter3d(x=xe,y=ye,z=ze,mode='lines',line=dict(color='black', width=1),name='Edges'))

    fig.update_layout(
        title=f"3D Surface with Wireframe-{organ_name}",
        scene=dict(
            xaxis_title='X(mm)',
            yaxis_title='Y(mm)',
            zaxis_title='Z(mm)',
            aspectmode='data'
        ),
        width=900,
        height=750
    )


    fig.write_html(f"{organ_name}_surface.html")

def view_multi_part_mesh_3D(case_id, part_list, organ_name,level):
    #Same as previous, for  multipart organs,whatever method used
 
    _,verts,faces,_,_=make_multi_part_mesh_MC(case_id,part_list,organ_name,level)


    fig = go.Figure()

    for part,vert,face in zip(part_list,verts,faces):
        x,y,z=vert.T
        i,j,k=face.T

        # edges
        edges=[]
        for tri in face:
            pts=[tri[0],tri[1],tri[2],tri[0]]
            for p1, p2 in zip(pts, pts[1:]):
                edges.append((vert[p1],vert[p2]))

        xe,ye,ze=[],[],[]
        for (p1, p2) in edges:
            xe+=[p1[0],p2[0],None]
            ye+=[p1[1],p2[1],None]
            ze+=[p1[2],p2[2],None]

        fig.add_trace(go.Mesh3d(
            x=x,y=y,z=z,
            i=i,j=j,k=k,
            color='magenta',
            opacity=0.5,
            name=part,
            lighting=dict(ambient=0.3,diffuse=0.6),
            lightposition=dict(x=100,y=100,z=100)
        ))
        #fig.add_trace(go.Scatter3d(x=xe,y=ye,z=ze,mode='lines',line=dict(color='black', width=1),name=part))
    

    fig.update_layout(
        title=f"3D Multi-Part Organ Surface-{organ_name}",
        scene=dict(
            xaxis_title='X(mm)',
            yaxis_title='Y(mm)',
            zaxis_title='Z(mm)',
            aspectmode='data'
        ),
        width=1000,
        height=800
    )

    fig.write_html(f"{organ_name}_surface.html")

