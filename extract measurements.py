import oiffile
import matplotlib.pyplot as plt
import numpy as np
import os as os
from glob import glob
import meshlib.mrmeshpy as mr
import meshlib.mrmeshnumpy as mrn
import skimage
from math import isnan
from skimage.measure import label, regionprops
import pandas as pd


def surface_area_3d(binary_volume, voxel_size):
    """
    Calculate the surface area of a 3D binary object using marching cubes.

    Parameters:
        binary_volume (3D np.array): binary array, 1 = object, 0 = background
        voxel_size (tuple): (z, y, x) voxel spacing in real-world units (e.g., microns)

    Returns:
        float: surface area in voxel_size units² (e.g., µm²)
    """
    verts, faces, _, _ = marching_cubes(binary_volume.astype(float), level=0.5, spacing=voxel_size)

    # Compute triangle areas
    v0 = verts[faces[:, 0]]
    v1 = verts[faces[:, 1]]
    v2 = verts[faces[:, 2]]

    # Area of each triangle = 0.5 * ||(v1-v0) × (v2-v0)||
    tri_areas = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    surface_area = tri_areas.sum()
    return surface_area




def mainscript(folder_name):
    os.chdir(folder_name)
    oib_files= glob('*.oib')
    dataframe=pd.DataFrame(columns = ['name','num','Actin vol','nucl vol','SA','sphere','major_axis length','orientation','touching boundary','image vol'])
    for f_name in oib_files:
        #import original oib image and voxel size
        with oiffile.OifFile(f_name) as oif:
            xstep=(oif.mainfile['Axis 0 Parameters Common']['EndPosition']-oif.mainfile['Axis 0 Parameters Common']['StartPosition'])/oif.mainfile['Axis 0 Parameters Common']['GUI MaxSize']
            ystep=(oif.mainfile['Axis 1 Parameters Common']['EndPosition']-oif.mainfile['Axis 1 Parameters Common']['StartPosition'])/oif.mainfile['Axis 1 Parameters Common']['GUI MaxSize']
            zstep=np.abs((oif.mainfile['Axis 3 Parameters Common']['EndPosition']-oif.mainfile['Axis 3 Parameters Common']['StartPosition'])/(oif.mainfile['Axis 3 Parameters Common']['GUI MaxSize']*1000))
            imag_org=(oif.asarray()/16)[:,:,:,:]
 
        #import cell outline
        cell_outline=np.load(f_name+'.npy')
        actin,nucl,num=cell_outline['array1'],cell_outline['array2']
        
        shap=np.shape(actin)


        props=measure.regionprops(actin)
        num,sa,act_vol,nu_vol,MAJ=[],[],[],[],[]


        for i in range(0,np.max(actin)-2):
            num.append(i)

            bounds=props[i]['bbox']
            act_mask=(actin==(i+2))[bounds[0]:bounds[3],bounds[1]:bounds[4],bounds[2]:bounds[5]]
            nucl_mask(nucl==(i+2))[bounds[0]:bounds[3],bounds[1]:bounds[4],bounds[2]:bounds[5]]
            actp=imag_org[int(2),bounds[0]:bounds[3],bounds[1]:bounds[4],bounds[2]:bounds[5]]*act_mask
            smap=imag_org[int(1),bounds[0]:bounds[3],bounds[1]:bounds[4],bounds[2]:bounds[5]]*act_mask

            # volume calcualtion of action and nucleus
            act_vol.append(np.sum(mask)*xstep*ystep*zstep)
            nu_vol.append(np.sum(nucl_mask)*xstep*ystep*zstep)

            #image volume for density calculation
            volume.append(xstep*ystep*zstep*shap[1]*shap[2]*shap[3])

            #calculate major axis length
            MAJ.append(props[i]['axis_major_length'])
            
            #orientation of cell based on major axis direction
            ori=props[i].orientation

            #check if cell is touching the boundary True = touching boundary False= not touching boundary
            if bounds[0] == 0 or bounds[1] == 0 or bounds[2] == 0 or bounds[3]== shap[0] or bounds[4]== shap[1] or bounds[5]== shap[2]:
                touch.append(True)
            else:
                touch.append(False)

            #calculate surface area and spheroid approximation
            sa.append(surface_area_3d(actp!=0, [zstep,xstep,ystep]))
            sphere.append((np.pi**(1/3)*(6*acts[i])**(2/3))/sa[i])

         
        df1=pd.DataFrame(zip(act_vol,num,act_vol,nu_vol,sa,MAJ,ori,touch,volume), columns = ['name','num','Actin vol','nucl vol','SA','sphere','major_axis length','orientation','touching boundary','image vol'])
        df1=df1.assign(name=f_name)
        #append data frames from different images together
        dataframe=pd.concat([df,df1], ignore_index = True)

    with pd.ExcelWriter('outputalot.xlsx') as writer:
        dataframe.to_excel(writer)   



mainscript(r"directory...")


