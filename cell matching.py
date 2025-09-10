import oiffile
import numpy as np
import os as os
from glob import glob
from skimage import data, filters, measure, morphology
from scipy import ndimage
from skimage import measure, color, io, util, measure, exposure, img_as_ubyte
from skimage.segmentation import clear_border, find_boundaries
import napari
import json
import tifffile
from scipy import ndimage as ndi
import pandas as pd
from skimage.filters import threshold_otsu, gaussian

def split_cells(cells,numb,min_volume):
    #splits integer array of segemented image into a binary array where each cell have its own ZxYxX array
    cellsplit=[]
    for i in range(0,numb):
        single=(cells==(i+1)).astype(bool)
        #filter out objects that are two small and to large based on pixel measurements
        if (np.sum(single))>min_volume and np.sum(single)<1000000:
            cellsplit.append(single)
    return cellsplit

def seperate(actin,nucli):
    #separte Actin objects that is associated with multiple DAPI objects, assigning each actin pixel to each DAPI object based on shortest euclidian distance.
    # loop through all actin objects 
    for cell in range(2,np.max(actin)):
        imgc=nucli==(cell)
        label=measure.label(imgc)
        #seperate cells only if they have multiple DAPI objects
        if (np.max(label))>1:
             img=actin==(cell)
             imgc=nucli==(cell)
             max=int(np.max(actin))
             nucli=nucli+(label>0)*(-cell)+(label==1)*(cell)+(label>1)*(max+label-1)
             props=measure.regionprops(label)
             label1=measure.label(img)
             props1=measure.regionprops(label1)
             bbox=props1[0]['bbox']
             bounding_box_min = np.array([bbox[0], bbox[1], bbox[2]])
             bounding_box_max = np.array([bbox[3], bbox[4], bbox[5]])

             pos=[]
             points=[]
             X,Y,Z = np.indices((bbox[3]-bbox[0], bbox[4]-bbox[1], bbox[5]-bbox[2]))
             Z= Z + bbox[2]
             Y = Y + bbox[1]
             X = X + bbox[0] 
             #assign pixels to each DAPI objects in the bounding box based of the shortest euclidean distance    
             for nucl in range(0,np.max(label)):
                 pos=np.array([props[nucl]['centroid'][0],props[nucl]['centroid'][1],props[nucl]['centroid'][2]])
                 points.append((X-pos[0])**2 + (Y-pos[1])**2 + (Z-pos[2])**2)
                     
             closest_indices=img*0
             closest_indices[bbox[0]:bbox[3],bbox[1]:bbox[4],bbox[2]:bbox[5]] = (np.argmin(points, axis=0)+1)*(img[bbox[0]:bbox[3],bbox[1]:bbox[4],bbox[2]:bbox[5]]==1)
             actin=actin+(max+closest_indices-cell-1)*(closest_indices>1)*(img==1)
    return([actin,nucli])



def DAPI_Actin_match(Actin,DAPI):
    #matches Actin and DAPI objects together discarding any objects that don't have a match, create a new Aactin and DAPI where matched objects in eacha array
    #have the same integer value, can assign multiple DAPI objects to single Actin
  
    match_Actin=np.zeros([shap[1],shap[2],shap[3]]).astype(bool)
    match_DAPI=np.zeros([shap[1],shap[2],shap[3]]).astype(bool)
    count=1

    #cycle through all actin objects
    for a in range(0,len(actin)):
        #calculate bounding box of actin object
        label=measure.label(actin[a])
        props=measure.regionprops(label)  
        bbox=props[0]['bbox']
        
        nucl_index=[]
        if_nucl=False
        nucl=np.zeros([shap[1],shap[2],shap[3]]).astype(bool)

        #cycle through all DAPI objects matching Actin and DAPI objects if manders overlap is greater then 0.6 based based on DAPI object as the mask 
        for n in range(0,len(nucl)):
            if (np.sum(Actin[a][bbox[0]:bbox[3],bbox[1]:bbox[4],bbox[2]:bbox[5]]*DAPI[n][bbox[0]:bbox[3],bbox[1]:bbox[4],bbox[2]:bbox[5]])/np.sum(DAPI[n]))>0.6:
                 if_nucl=True
                 nucl=nucl+DAPI[n]
                 nucl_index.append(n)
        
        if if_nucl:
            #assign DAPI and Actin objects to array giving them the same integer
            count=count+1
            match_Actin=match_Actin+Actin[a]*count
            nuclc=match_DAPI+nucl*count

        #delete DAPI objects assigned from the list
        nucl=[i for j, i in enumerate(DAPI) if j not in np.array(nucl_index)]
        
    return match_Actin,match_DAPI



def threshold(actin,props,araw): 
    newact=np.zeros_like(actin)
    #loop through all actin objects
    for i in range(0,np.max(actin)-1):
        bounds=props[i]['bbox']
        mask=araw[bounds[0]:bounds[3],bounds[1]:bounds[4],bounds[2]:bounds[5]]
        #otsu threshold based on bounding box around actin object
        thresh = threshold_otsu(mask)
        mask = mask > thresh
        mask=mask*(actin==(i+2))[bounds[0]:bounds[3],bounds[1]:bounds[4],bounds[2]:bounds[5]]
        #fill in holes in 2D on each plane of the object
        for plane in range(0,bounds[3]-bounds[0]):
            mask[plane,:,:]=ndimage.binary_fill_holes(mask[plane,:,:])
        newact[bounds[0]:bounds[3],bounds[1]:bounds[4],bounds[2]:bounds[5]]=newact[bounds[0]:bounds[3],bounds[1]:bounds[4],bounds[2]:bounds[5]]+mask*(i+2)  
    return newact  
    



os.chdir(r"Folder directory...")
oib_files=glob('*.oib')
#loop through all oib files in folder
for f_name in oib_files:
    # read in orginal oib file 
    with oiffile.OifFile(f_name) as oif:
        xstep=(oif.mainfile['Axis 0 Parameters Common']['EndPosition']-oif.mainfile['Axis 0 Parameters Common']['StartPosition'])/oif.mainfile['Axis 0 Parameters Common']['GUI MaxSize']
        ystep=(oif.mainfile['Axis 1 Parameters Common']['EndPosition']-oif.mainfile['Axis 1 Parameters Common']['StartPosition'])/oif.mainfile['Axis 1 Parameters Common']['GUI MaxSize']
        zstep=np.abs((oif.mainfile['Axis 3 Parameters Common']['EndPosition']-oif.mainfile['Axis 3 Parameters Common']['StartPosition'])/(oif.mainfile['Axis 3 Parameters Common']['GUI MaxSize']*1000))
        imag1=(oif.asarray()/16)

    # read in binary TIFF file of thresholded image
    imag=tifffile.imread(f_name+'.tif')>0
    shap=np.shape(imag)
    
    #segment DAPI channels into distinct objects 
    DAPI,numb=measure.label(imag[0,:,:,:],return_num=True)
    #split segmented image into ZxXxYxcell binary array 
    DAPI=split_cells(DAPI,numb,200)

    #segment Actin channels into distinct objects 
    Actin,numb=measure.label(imag[2,:,:,:],return_num=True)
    #split segmented image into ZxXxYxcell binary array 
    Actin=split_cells(Actin,numb,200)

    #match Actin and DAPI objects are discard objects that do not have a match. 
    matched_actin,matched_nucl=DAPI_Actin_match(Actin,DAPI)

    matched_actin,matched_nucl=seperate(actin,nuclc)

    #gaussian blur image sigma=2
    guas_raw=gaussian(imag1[2,:,:,:],sigma=2)
    #rethreshold actin objects based on Otsu threshold of bounding box of each Actin Object
    actin=threshold(matched_actin,measure.regionprops(matched_actin) ,guas_raw)


    with open(f_name+'.npy', 'wb') as f:
        #save file as numpy array for matched actin skeletons and DAPI objects.
        np.savez(f, array1=np.array(matched_actin), array2=np.array(matched_nucl))


