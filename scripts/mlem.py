# Data Pathway
path = '/home/raptor/repos/'

# Libraries
import numpy as np
import tables
import time

starttime = time.time()

# Output Filename
imageName = 'image'

# Input Parameters
eps = 0.05 # regularization weighting factor
nIter = 50 # number of iterations for MLEM

# Input Files
detectorImages = np.array( [ path + 'gammaimage/data/spiral0deg.h5',
                            path + 'gammaimage/data/spiral45deg.h5',
                            path + 'gammaimage/data/spiral90deg.h5',
                            path + 'gammaimage/data/spiral135deg.h5',
                            path + 'gammaimage/data/spiral180deg.h5',
                            path + 'gammaimage/data/spiral225deg.h5',
                            path + 'gammaimage/data/spiral270deg.h5',
                            path + 'gammaimage/data/spiral315deg.h5'] )

systemMatrices = np.array( [ path + 'gammaimage/sysresp/sysresp122_0.h5',
                            path + 'gammaimage/sysresp/sysresp122_45.h5',
                            path + 'gammaimage/sysresp/sysresp122_90.h5',
                            path + 'gammaimage/sysresp/sysresp122_135.h5',
                            path + 'gammaimage/sysresp/sysresp122_180.h5',
                            path + 'gammaimage/sysresp/sysresp122_225.h5',
                            path + 'gammaimage/sysresp/sysresp122_270.h5',
                            path + 'gammaimage/sysresp/sysresp122_315.h5'] )

# Defining Source Voxel Space
sourceX, sourceY, sourceZ = np.mgrid[ -49:51:2, -69:71:2, -49:51:2 ]
sourcePixels = np.array( [ sourceX.flatten(), sourceY.flatten(), sourceZ.flatten() ] ).T

# Defining Detector Pixel Space
detectorX, detectorY = np.mgrid[ -36.5:37.5:1, -36.5:37.5:1 ]
detectorPixels = np.array( [ detectorX.flatten(), detectorY.flatten() ] ).T

#  Loading Detector Data
detectorData = np.zeros( ( len( detectorImages ), len( detectorPixels ) ) )
for i in np.arange( len( detectorImages ) ):
    detectorData[ i ] = tables.open_file( detectorImages[ i ], 'r' ).root.counts.read()

# Loading System Response and Computing Sensitivity
sensitivity = np.zeros( len( sourcePixels ), dtype = 'float32' )
systemMatrix = { str( 0 ): np.zeros( ( len( sourcePixels ), len( detectorPixels ) ), dtype = 'float32' ) }
for i in np.arange( len( detectorImages ) ):
    matrix = tables.open_file( systemMatrices[ i ], 'r').root.matrix.read()
    sensitivity = sensitivity + matrix.sum( axis = 1 )
    systemMatrix[ str( i ) ] = matrix

# Performing MLEM Image Reconstruction
lamb = np.ones( len( sourcePixels ) )
# for loop over iterations
for i in np.arange( nIter ):
    print( 'iteration: %i of %i' %( i + 1, nIter ) )
    outerSum = np.zeros( ( len( sourcePixels )  ), dtype = 'float32' )
    # for loop over projections
    for p in np.arange( len( detectorImages ) ):
        # for loop over detector pixels
        for d in np.arange( len( detectorPixels ) ):
            sysMat = systemMatrix[ str( p ) ][ :, d ]
            data = detectorData[ p, d ]
            innerSum = np.dot( sysMat, lamb )
            numerator = data * sysMat
            fraction = np.divide(numerator, innerSum, out = np.zeros_like( numerator ), where = innerSum != 0 )
            outerSum = outerSum + fraction
    lamb = lamb * ( sensitivity / ( sensitivity ** 2 + max( sensitivity ) ** 2 * eps ** 2 ) ) * outerSum

f = tables.open_file( path + 'gammaimage/images/' + imageName + '.h5', 'w')
f.create_array('/', 'image', lamb)
f.close()

print( 'Computation Time:', time.time() - starttime )
