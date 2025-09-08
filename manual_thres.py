import oiffile, os, napari,tifffile
import numpy as np
from skimage.filters import gaussian
from glob import glob
from scipy import ndimage as ndi
import napari

#sequntially Reads in all oib files from a folder, assuming data is arranged as CZXY and channel order is DAPI, alpha-SMA and Actin, and returns initial manual thersholding using napari as TIFF file. 

if __name__ == '__main__':
    os.chdir(r"folder_directory")
    oib_files = glob('*.oib')
    for f_name in oib_files:
        #read in OIB file
        with oiffile.OifFile(f_name) as oif:
            xstep=(oif.mainfile['Axis 0 Parameters Common']['EndPosition']-oif.mainfile['Axis 0 Parameters Common']['StartPosition'])/oif.mainfile['Axis 0 Parameters Common']['GUI MaxSize']
            ystep=(oif.mainfile['Axis 1 Parameters Common']['EndPosition']-oif.mainfile['Axis 1 Parameters Common']['StartPosition'])/oif.mainfile['Axis 1 Parameters Common']['GUI MaxSize']
            zstep=np.abs((oif.mainfile['Axis 3 Parameters Common']['EndPosition']-oif.mainfile['Axis 3 Parameters Common']['StartPosition'])/(oif.mainfile['Axis 3 Parameters Common']['GUI MaxSize']*1000))
            imag=np.float16((oif.asarray()/16)[:,:,:,:])
        shap=np.shape(imag)
        #gaussian filter sigma =2
        imag[:,:,:,:]=gaussian(imag[:,:,:,:],sigma=2,channel_axis=0)

        #get thresholding values that can be set manually using Napari
        viewer = napari.Viewer()
        image_layer = viewer.add_image(imag, channel_axis=0, name=['DAPI','SMA','Actin'],colormap=['blue','green','red'])
        viewer.layers['DAPI'].scale=[zstep,xstep,ystep]
        viewer.layers['SMA'].scale=[zstep,xstep,ystep]
        viewer.layers['Actin'].scale=[zstep,xstep,ystep]
        napari.run()

        #apply threshold and fill in holes on all planes
        for i in range(0,3):
            img=imag[i,:,:,:]
            img=img>viewer.layers[i].contrast_limits[0]
            for plane in range(0,np.shape(imag)[1]):
                img[plane,:,:] = ndi.binary_fill_holes(img[plane,:,:])
            imag[i,:,:,:]= img

        #save file as binary TIFF file
        tifffile.imwrite(f_name+'.tif',imag)




