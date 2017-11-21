import tomopy
import dxchange
import numpy as np
import matplotlib.pyplot as plt




# prj = np.random.rand(14, 15, 16).astype(np.float32)
# ang = tomopy.angles(prj.shape[0])

# # t = time.time()
# rec1, rec2, rec3 = tomopy.vector(prj, prj, ang, ang)
# dxchange.write_tiff(rec1)

mx = dxchange.read_tiff('/Users/dgursoy/Documents/Data/M4R1_mx.tif')
my = dxchange.read_tiff('/Users/dgursoy/Documents/Data/M4R1_my.tif')
mz = dxchange.read_tiff('/Users/dgursoy/Documents/Data/M4R1_mz.tif')

ang = tomopy.angles(31)
prj1 = tomopy.project2(mx, my, mz, ang, axis=1)
prj2 = tomopy.project2(mx, my, mz, ang, axis=2)
print (mx.shape, prj1.shape, prj2.shape)

rec1, rec2, rec3 = tomopy.vector(prj1, prj2, ang, ang)
dxchange.write_tiff(rec1)

