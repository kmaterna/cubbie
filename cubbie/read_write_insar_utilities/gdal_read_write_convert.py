
import numpy as np
from osgeo import gdal, osr
import matplotlib.pyplot as plt
from . import isce_read_write


def read_geotiff(filename):
    """
    :param filename: string, name of geotiff
    :return: xarray, yarray, zarray raster
    """
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    rb = ds.GetRasterBand(1)
    img_array = rb.ReadAsArray()
    transform = ds.GetGeoTransform()
    ylen, xlen = np.shape(img_array)
    xs, ys = isce_read_write.get_xarray_yarray_from_transform(transform, xlen, ylen)
    return xs, ys, img_array


def write_as_geotiff(x, y, zdata, out_tif):
    """
    :param x: 1d array associated with longitude values, in degrees
    :param y: 1d array associated with latitude values, in degrees
    :param zdata: 2d array of raster data
    :param out_tif: string, filename
    :return:
    """
    gdal.UseExceptions()  # opt-in now (preferred)

    # 2D array
    array = np.flipud(zdata)  # shape: (rows, cols)

    # GeoTIFF metadata
    x_min = np.min(x)  # top-left corner longitude
    y_max = np.max(y)  # top-left corner latitude
    xpixel_size, ypixel_size = x[1]-x[0], y[1]-y[0]  # resolution (pixel size)
    no_data_value = -9999
    epsg_code = 4326  # WGS-84

    # Create GeoTIFF
    rows, cols = array.shape
    driver = gdal.GetDriverByName('GTiff')  # :contentReference[oaicite:0]{index=0}
    ds = driver.Create(
        out_tif, cols, rows, 1, gdal.GDT_Float32,  # :contentReference[oaicite:1]{index=1}
        options=['COMPRESS=LZW', 'TILED=YES']
    )

    # Set geotransform (top left x, w-e pixel res, 0, top left y, 0, n-s pixel res)
    geotransform = (x_min, xpixel_size, 0, y_max, 0, -ypixel_size)
    ds.SetGeoTransform(geotransform)

    # Set spatial reference (optional, here WGS84)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg_code)
    ds.SetProjection(srs.ExportToWkt())

    # Write data to band 1
    band = ds.GetRasterBand(1)
    band.WriteArray(array)
    band.SetNoDataValue(no_data_value)
    band.FlushCache()

    print("Writing file %s " % out_tif)

    # Clean up
    ds = None  # Closes and saves the file
    srs = None

    return


def write_colored_geotiff(x, y, zdata, out_tif, cmap='viridis'):
    """
    Write raster data into a colored geoTiff. In the process, the data will be truncated to exactly 256 values.

    :param x: 1d array associated with longitude values, in degrees
    :param y: 1d array associated with latitude values, in degrees
    :param zdata: 2d array of raster data
    :param out_tif: string, filename
    :param cmap: name of a valid Python color map
    :return:
    """
    gdal.UseExceptions()  # opt-in now (preferred)

    # Sample 2D array
    array = np.flipud(zdata)  # shape: (rows, cols)

    # GeoTIFF metadata
    x_min = np.min(x)  # top-left corner longitude
    y_max = np.max(y)  # top-left corner latitude
    xpixel_size, ypixel_size = x[1]-x[0], y[1]-y[0]  # resolution (pixel size)
    epsg_code = 4326  # WGS-84

    # Create GeoTIFF
    rows, cols = array.shape

    # Set geotransform (top left x, w-e pixel res, 0, top left y, 0, n-s pixel res)
    geotransform = (x_min, xpixel_size, 0, y_max, 0, -ypixel_size)

    # Set spatial reference (optional, here WGS84)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg_code)

    # Write geotiff with colors: Normalize data to 0–255
    data_norm = (255 * (array - array.min()) / (array.max() - array.min())).astype(np.uint8)

    # Apply a matplotlib colormap (e.g., viridis)
    cmap = plt.get_cmap(cmap)
    rgba_img = (cmap(data_norm / 255.0)[:, :, :3] * 255).astype(np.uint8)

    # Save as 3-band RGB GeoTIFF
    driver = gdal.GetDriverByName('GTiff')
    rgb_ds = driver.Create(out_tif, cols, rows, 3, gdal.GDT_Byte)

    for i in range(3):  # R, G, B
        rgb_ds.GetRasterBand(i + 1).WriteArray(rgba_img[:, :, i])

    # Optional: set geotransform and projection if known
    rgb_ds.SetGeoTransform(geotransform)
    rgb_ds.SetProjection(srs.ExportToWkt())

    print("Writing file %s " % out_tif)

    # Clean up
    rgb_ds.FlushCache()
    rgb_ds = None
    srs = None
    return


def convert_tiff_to_kml():
    """
    Convert a geotiff file into a KML. Not yet written.
    :return:
    """
    return


def convert_tiff_to_kmz(tif_name, kmz_name):
    """
    Convert a geotiff file into a KMZ

    :param tif_name: string, name of geoTiff file
    :param kmz_name: string, name of KMZ file
    :return:
    """
    gdal.Translate(
        kmz_name, tif_name,
        format='KMLSUPEROVERLAY',
        creationOptions=[
            'FORMAT=PNG', 'KMZ=YES'
        ]
    )
    return
