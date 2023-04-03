import os
from statistics import median
from sys import builtin_module_names
import numpy as np
import rasterio as ra
from rasterio.plot import show
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import earthpy as et
import earthpy.spatial as es
import earthpy.plot as ep
from skimage import exposure
from osgeo import gdal
from PIL import Image
from earthpy.io import path_to_example
from mmap import PAGESIZE
from matplotlib.pyplot import draw
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Frame, PageBreak
from reportlab.platypus import Image as img
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.lib.styles import ParagraphStyle

my_path= r'C:\Users\Admin420\Desktop\NDVI\BandAnalysis.pdf'
NDVImap = r'C:\Users\Admin420\Desktop\NDVI\NDVImap.jpg'
NDVIhistogram = r'C:\Users\Admin420\Desktop\NDVI\NDVIhistogram.jpg'
FCCimg = r'C:\Users\Admin420\Desktop\NDVI\FCCimg.jpg'
RGBimg = r'C:\Users\Admin420\Desktop\NDVI\RGBimg.jpg'

output_NDVI = r'C:\Users\Admin420\Desktop\NDVI\ndvi.tif'
output_RGB = r'C:\Users\Admin420\Desktop\NDVI\rgb.tif'
output_FCC = r'C:\Users\Admin420\Desktop\NDVI\fcc.tif'

BLUE = "C:/Users/Admin420/Desktop/NDVI/S2A_MSIL2A_20220701T030531_N0400_R075_T50UPA_20220701T064612.SAFE/GRANULE/L2A_T50UPA_A036680_20220701T031319/IMG_DATA/R10m/T50UPA_20220701T030531_B02_10m.jp2"
GREEN = "C:/Users/Admin420/Desktop/NDVI/S2A_MSIL2A_20220701T030531_N0400_R075_T50UPA_20220701T064612.SAFE/GRANULE/L2A_T50UPA_A036680_20220701T031319/IMG_DATA/R10m/T50UPA_20220701T030531_B03_10m.jp2"
RED = "C:/Users/Admin420/Desktop/NDVI/S2A_MSIL2A_20220701T030531_N0400_R075_T50UPA_20220701T064612.SAFE/GRANULE/L2A_T50UPA_A036680_20220701T031319/IMG_DATA/R10m/T50UPA_20220701T030531_B04_10m.jp2"
NIR = "C:/Users/Admin420/Desktop/NDVI/S2A_MSIL2A_20220701T030531_N0400_R075_T50UPA_20220701T064612.SAFE/GRANULE/L2A_T50UPA_A036680_20220701T031319/IMG_DATA/R10m/T50UPA_20220701T030531_B08_10m.jp2"

with ra.open(BLUE) as src:              
    band_blue = src.read()

with ra.open(GREEN) as src:
    band_green = src.read()

with ra.open(RED) as src:                    
    band_red = src.read()

with ra.open(NIR) as src:
    band_nir = src.read()

band2=ra.open(BLUE)
band3=ra.open(GREEN)
band4=ra.open(RED)
band8=ra.open(NIR)

np.seterr(divide='ignore', invalid='ignore')       

ndvi = (band_nir.astype(float) - band_red.astype(float)) / (band_nir + band_red) 

max_NDVI = np.nanmax(ndvi)      
min_NDVI = np.nanmin(ndvi)

mean_NDVI = np.nanmean(ndvi)      
median_NDVI = np.nanmedian(ndvi)

profile = src.meta                               
profile.update(driver='GTiff')
profile.update(dtype=ra.float32)

with ra.open(output_NDVI, 'w', **profile) as dst:     
    dst.write(ndvi.astype(ra.float32))

rgbprofile= band2.profile
rgbprofile.update(driver="GTiff", count=3, photometric="RGB", compress="LZW", bigtiff="IF_NEEDED", tiled=True, blockxsize=256, blockysize=256)

with ra.open(output_RGB, "w", **rgbprofile) as rgb:
    rgb.write(band4.read(1), 1)
    rgb.write(band3.read(1), 2)
    rgb.write(band2.read(1), 3)

with ra.open(output_FCC, "w", **rgbprofile) as rgb:
    rgb.write(band8.read(1), 1)
    rgb.write(band4.read(1), 2)
    rgb.write(band3.read(1), 3)

with ra.open(output_NDVI) as src:
    dem = src.read()
    fig, ax = plt.subplots(figsize = (10, 5))

im = ax.imshow(dem.squeeze())
ep.colorbar(im)
ax.set(title="NDVI MAP")
ax.set_axis_off()
plt.savefig("NDVImap.jpg")

ep.hist(ndvi,figsize=(12, 6), title=[""]) # CREATING A HISTOGRAM OF ndvi VALUE spread to better understand the data
plt.savefig("NDVIhistogram.jpg")

rgb_img = ra.open(output_RGB)
image = np.array([rgb_img.read(1), rgb_img.read(2), rgb_img.read(3)]).transpose(1,2,0)
p2, p98 = np.percentile(image, (2,98))
rgb_image = exposure.rescale_intensity(image, in_range=(p2, p98)) / 100000

im = Image.fromarray((rgb_image * 255).astype(np.uint8))
im.save("RGBimg.jpg")

fcc_img = ra.open(output_FCC)
image = np.array([fcc_img.read(1), fcc_img.read(2), fcc_img.read(3)]).transpose(1,2,0)
p2, p98 = np.percentile(image, (2,98))
r_channel = exposure.rescale_intensity(image[:,:,0], in_range=(p2, p98))/100000
g_channel = exposure.rescale_intensity(image[:,:,1], in_range=(p2, p98))/20000
b_channel = exposure.rescale_intensity(image[:,:,2], in_range=(p2, p98))/100000
auto = np.stack((r_channel,g_channel,b_channel),axis=2)

im = Image.fromarray((auto* 255).astype(np.uint8))
im.save("FCCimg.jpg")

styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='Justify', alignment=TA_CENTER, fontName='Times-Bold', fontSize=12))
styles.add(ParagraphStyle(name='large', alignment=TA_CENTER, fontName='Times', fontSize=20))

c = canvas.Canvas('Band_Analysis.pdf')
f = Frame(inch, inch, 7*inch, 10*inch, showBoundary=0)
f2= Frame(inch, inch, 7*inch, 10*inch, showBoundary=0)

c.line(50, 760, 560, 760)

Story=[]

ptext = 'Band Analysis'
Story.append(Paragraph(ptext, styles['large']))
Story.append(Spacer(1, 25))

im = img(NDVImap, 7*inch, 3.5*inch)
Story.append(im)
Story.append(Spacer(1, 1))
ptext = 'Normalized Difference Vegetation Index (NDVI)'
Story.append(Paragraph(ptext, styles["Justify"]))
Story.append(Spacer(1, 12))

im2 = img(NDVIhistogram, 6*inch, 3*inch )
Story.append(im2)
Story.append(Spacer(1, 1))
ptext = 'NDVI Histogram (Distribution of Pixels)'
Story.append(Paragraph(ptext, styles["Justify"]))
Story.append(Spacer(1, 12))

Story2=[]

im3= img(RGBimg, 4*inch, 4*inch )
Story2.append(im3)
ptext = 'Red, Green, Blue (RGB) MAP'
Story2.append(Paragraph(ptext, styles["Justify"]))
Story2.append(Spacer(1, 40))

im4 = img(FCCimg, 4*inch, 4*inch )
Story2.append(im4)
ptext = 'False Color Composite(FCC) MAP'
Story2.append(Paragraph(ptext, styles["Justify"]))

f.addFromList(Story,c)
c.showPage()
f2.addFromList(Story2,c)
c.save()