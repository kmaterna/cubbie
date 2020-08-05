"""
# A set of functions that read and write vrt gdal grid files
# Compatible with ISCE. 
    read_complex_data(GDALfilename)
    read_scalar_data(GDALfilename)
    read_phase_data(GDALfilename)
    write_isce_data(data, nx, ny, dtype, filename)
    plot_scalar_data(GDALfilename,band=1,title="",colormap='gray',aspect=1....) 
    plot_complex_data(GDALfilename,title="",aspect=1,....)
"""

import numpy as np 
import matplotlib.pyplot as plt 
import struct


# ----------- READING FUNCTIONS ------------- # 

def read_complex_data(GDALfilename):
    # Reads data into a 2D array where each element is a complex number. 
    from osgeo import gdal            ## GDAL support for reading virtual files
    print("Reading file %s " % GDALfilename);
    ds = gdal.Open(GDALfilename, gdal.GA_ReadOnly)
    slc = ds.GetRasterBand(1).ReadAsArray()
    transform = ds.GetGeoTransform()
    ds = None
    
    # getting the min max of the axes
    firstx = transform[0]
    firsty = transform[3]
    deltay = transform[5]
    deltax = transform[1]
    lastx = firstx+slc.shape[1]*deltax
    lasty = firsty+slc.shape[0]*deltay
    ymin = np.min([lasty,firsty])
    ymax = np.max([lasty,firsty])
    xmin = np.min([lastx,firstx])
    xmax = np.max([lastx,firstx])

    # put all zero values to nan
    try:
        slc[slc==0]=np.nan
    except:
        pass

    return slc;

def read_scalar_data(GDALfilename, band=1,flush_zeros=True):
    # band = 1;  # this seems right for most applications 
    # For unwrapped files, band = 2
    from osgeo import gdal            ## GDAL support for reading virtual files
    print("Reading file %s " % GDALfilename);
    if ".unw" in GDALfilename and ".unw." not in GDALfilename and band==1:
        print("WARNING: We usually read band=2 for snaphu unwrapped files. Are you sure you want band 1 ????");
    ds = gdal.Open(GDALfilename, gdal.GA_ReadOnly)
    data = ds.GetRasterBand(band).ReadAsArray()
    transform = ds.GetGeoTransform()
    ds = None
    
    # getting the min max of the axes
    firstx = transform[0]
    firsty = transform[3]
    deltay = transform[5]
    deltax = transform[1]
    lastx = firstx+data.shape[1]*deltax
    lasty = firsty+data.shape[0]*deltay
    ymin = np.min([lasty,firsty])
    ymax = np.max([lasty,firsty])
    xmin = np.min([lastx,firstx])
    xmax = np.max([lastx,firstx])

    # put all zero values to nan
    if flush_zeros:
        try:
            data[data==0]=np.nan
        except:
            pass

    return data;

def read_phase_data(GDALfilename):
    # Start with a complex quantity, and return only the phase of that quantity. 
    slc = read_complex_data(GDALfilename);
    phasearray = np.angle(slc);
    return phasearray;

def read_scalar_data_no_isce(filename, nx, ny):
    # Takes float32 numbers from binary file into 2d array
    final_shape=(ny,nx);
    num_data = nx*ny;
    print("Reading file %s into %d x %d array" % (filename, ny, nx));
    f = open(filename,'rb')
    rawnum = f.read();
    f.close();
    floats = np.array(struct.unpack('f'*num_data, rawnum));
    scalar_field = floats.reshape(final_shape);
    return scalar_field;

def read_phase_data_no_isce(filename, nx, ny):
    final_shape=(ny,nx);
    num_data = nx*ny*2;
    print("Reading file %s into %d x %d array" % (filename, ny, nx));
    f = open(filename,'rb')
    rawnum = f.read();    
    floats = np.array(struct.unpack('f'*num_data, rawnum))
    f.close();
    real = floats[::2];
    imag = floats[1::2];
    phase = np.arctan2(imag,real);
    phase = phase.reshape(final_shape);
    return phase;

# If you don't have ISCE
def read_unw_data_from_xml(filename):
    # For the unw format produced by isce (one minor change from that format)
    # (BIL scheme assumed)
    # There must be a matching xml in this directory
    # --------------------

    # Parse xml in a slightly manual fashion, looking for length and width
    xml_file = filename+".xml";
    tree = ET.parse(xml_file);
    root = tree.getroot();  # open xml file
    for element in root:  # we can index through the root
        if element.attrib['name']=="coordinate1":
            for subE in element:
                if len(subE)>0:
                    if subE.attrib['name']=='size':
                        ncols = int(subE[0].text);
        if element.attrib['name']=="coordinate2":
            for subE in element:
                if len(subE)>0:
                    if subE.attrib['name']=='size':
                        nrows = int(subE[0].text);

    # Open the binary file
    f = open(filename,'rb');
    final_shape=(nrows,ncols*2);  # unw has two bands with BIL scheme
    num_data = final_shape[0]*final_shape[1];
    rawnum = f.read();
    floats = np.array(struct.unpack('f'*num_data, rawnum))
    data = floats.reshape(final_shape); 
    f.close();
    return data;

# ----------- WRITING FUNCTIONS ------------- # 

def write_isce_data(data, nx, ny, dtype, filename, 
    firstLat=None, firstLon=None, deltaLon=None, deltaLat=None,Xmin=None, Xmax=None):
    # This function writes ISCE data into a single-band file with given filename
    # Plus creating an associated .vrt and .xml file
    # If DTYPE=="FLOAT": you're writing scalar data (float32)
    # IF DTYPE=="CFLOAT": you're writing complex data (float32 + j*float32)
    import isce
    import isceobj    
    from osgeo import gdal
    print("Writing data as file %s " % filename);
    out = isceobj.createImage()
    out.setFilename(filename)
    out.setWidth(nx)
    out.setLength(ny)
    out.setInterleavedScheme('BIP') #'BIP'/ 'BIL' / ‘BSQ’
    out.setAccessMode('READ')
    out.setDataType(dtype)
    if firstLon is not None:  # Special options that aren't usually used. 
        out.setFirstLongitude(firstLon);
    if firstLat is not None:
        out.setFirstLatitude(firstLat);
    if deltaLon is not None:
        out.setDeltaLongitude(deltaLon);
    if deltaLat is not None:
        out.setDeltaLatitude(deltaLat);
    if Xmin is not None:
        out.setXmin(Xmin);
    if Xmax is not None:
        out.setXmax(Xmax);
    out.renderHdr()
    data.tofile(filename) # write file out
    return

def write_isce_unw(data1, data2, nx, ny, dtype, filename,
    firstLat=None, firstLon=None, deltaLon=None, deltaLat=None,Xmin=None, Xmax=None):
    # ISCE uses band=2 for the unwrapped phase of .unw files
    # Writes to float32
    import isce
    import isceobj
    from osgeo import gdal    
    print("Writing data as file %s " % filename);
    out = isceobj.Image.createUnwImage()     
    out.setFilename(filename)
    out.setWidth(nx)
    out.setLength(ny)
    out.imageType = 'unw'
    out.bands = 2
    out.scheme="BIL"
    out.setAccessMode('read')
    if firstLon is not None:  # Special options that aren't usually used. 
        out.setFirstLongitude(firstLon);
    if firstLat is not None:
        out.setFirstLatitude(firstLat);
    if deltaLon is not None:
        out.setDeltaLongitude(deltaLon);
    if deltaLat is not None:
        out.setDeltaLatitude(deltaLat);
    if Xmin is not None:
        out.setXmin(Xmin);
    if Xmax is not None:
        out.setXmax(Xmax);    
    out.setDataType(dtype)
    out.renderHdr()
    data_to_file_2_bands(data1, data2, filename); #  dump the data into a binary file
    return;

def data_to_file_2_bands(data1, data2, filename):
    data1=np.float32(data1);  # we should be consistent about float types here. 
    data2=np.float32(data2);
    data=np.hstack((data1,data2));  # establishing two bands
    data.tofile(filename)        
    return;

def data_to_file_1_bands(data1, filename):
    data1=np.float32(data1);  # we should be consistent about float types here. 
    data1.tofile(filename)        
    return;

def plot_scalar_data(GDALfilename,band=1,title="",colormap='gray',aspect=1, 
    datamin=None, datamax=None,draw_colorbar=True,colorbar_orientation="horizontal",background=None, outname=None):
    from osgeo import gdal    
    ds = gdal.Open(GDALfilename, gdal.GA_ReadOnly)
    data = ds.GetRasterBand(band).ReadAsArray()
    transform = ds.GetGeoTransform()
    ds = None
    
    # getting the min max of the axes
    # Note: this assumes that the transform is north-up
    # There are transform[2] and transform[4] for other projections. 
    # In general, I might need a more general visualizer. 
    firstx = transform[0]
    firsty = transform[3]
    deltay = transform[5]
    deltax = transform[1]
    lastx = firstx+data.shape[1]*deltax
    lasty = firsty+data.shape[0]*deltay
    ymin = np.min([lasty,firsty])
    ymax = np.max([lasty,firsty])
    xmin = np.min([lastx,firstx])
    xmax = np.max([lastx,firstx])

    # put all zero values to nan and do not plot nan
    if background is None:
        try:
            data[data==0]=np.nan
        except:
            pass
    
    fig = plt.figure(figsize=(18, 16))
    ax = fig.add_subplot(111)
    cax = ax.imshow(data, vmin = datamin, vmax=datamax, cmap=colormap,extent=[xmin,xmax,ymin,ymax])
    ax.set_title(title)
    if draw_colorbar is not None:
        cbar = fig.colorbar(cax,orientation=colorbar_orientation)
    ax.set_aspect(aspect)    
    if outname is None:
        plt.show()
    else:
        fig.savefig(outname);
    return;


def plot_complex_data(GDALfilename,title="",aspect=1, band=1, colormap='rainbow',
    datamin=None, datamax=None,draw_colorbar=None,colorbar_orientation="horizontal", outname=None):
    from osgeo import gdal   
    ds = gdal.Open(GDALfilename, gdal.GA_ReadOnly)
    slc = ds.GetRasterBand(band).ReadAsArray()
    transform = ds.GetGeoTransform()
    ds = None
    
    # getting the min max of the axes
    firstx = transform[0]
    firsty = transform[3]
    deltay = transform[5]
    deltax = transform[1]
    lastx = firstx+slc.shape[1]*deltax
    lasty = firsty+slc.shape[0]*deltay
    ymin = np.min([lasty,firsty])
    ymax = np.max([lasty,firsty])
    xmin = np.min([lastx,firstx])
    xmax = np.max([lastx,firstx])

    # put all zero values to nan and do not plot nan
    try:
        slc[slc==0]=np.nan
    except:
        pass

    
    fig = plt.figure(figsize=(18, 16))
    ax = fig.add_subplot(1,2,1)
    cax1=ax.imshow(np.abs(slc),vmin = datamin, vmax=datamax, cmap='gray',extent=[xmin,xmax,ymin,ymax])
    ax.set_title(title + " (amplitude)")
    if draw_colorbar is not None:
        cbar1 = fig.colorbar(cax1,orientation=colorbar_orientation)
    ax.set_aspect(aspect)

    ax = fig.add_subplot(1,2,2)
    cax2 =ax.imshow(np.angle(slc),cmap=colormap,extent=[xmin,xmax,ymin,ymax])
    ax.set_title(title + " (phase [rad])")
    if draw_colorbar is not None:
        cbar2 = fig.colorbar(cax2,orientation=colorbar_orientation)
    ax.set_aspect(aspect)
    if outname is None:
        plt.show()
    else:
        fig.savefig(outname);

    return;
