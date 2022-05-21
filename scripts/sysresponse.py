# Data Pathway
path = '/home/raptor/'

# Libraries
import numpy as np
import time
import tables
import sys

name = 'gammaimage/sysresp/sysresp_'

# System Response Adjustable Parameters
poses = np.array( [ 0, 45, 90, 135, 180, 225, 270, 315 ] ) * np.pi / 180 # Rotation angles (radians), i.e. detector poses
T = 2.4 # mask thickness (mm)
D = 115 # focal distance (mm), i.e. distance between detector and mask
mu = 2.672 * 19.3 * 0.1 # attenuation factor in tungsten (mm-1), 2.672 at 122-keV, 0.6422 at 218-keV, 0.9321 at 186-keV
c = np.array( [ 1, 0, 95] ) # center of body rotation (mm), # 97 for all Ac-225, 117 for Co-57 flexible
n = 1  # detector pixel size

# Initializing Source Voxels
sourceX, sourceY, sourceZ = np.mgrid[ -49:51:2, -69:71:2, -49:51:2 ]
sourcePixels = np.array( [ sourceX.flatten(), sourceY.flatten(), sourceZ.flatten() ] ).T

# Initializing Detector Pixels
detectorX, detectorY = np.mgrid[ -36.5:37.5:n, -36.5:37.5:n ]
detectorPixels = np.array( [ detectorX.flatten(), detectorY.flatten() ] ).T

# Loading and Oversampling Efficiency Data
f1 = tables.open_file( path + 'gammaimage/efficiency.h5', 'r' ).root.eff.read()
f2 = f1.repeat( int( 2 ), axis = 0 )
f3 = f2.reshape( int( np.sqrt( len( f1 ) ) ), int( 2 * np.sqrt( len( f1 ) ) ) )
eff = f3.repeat( int( 2 ), axis = 0 )

# Loading Far-field Lookup Table
lookTab = tables.open_file( path + 'gammaimage/farfieldmap.h5', 'r' ).root.Matrix.read()

# Initilizing Original Cartesian Far-Field Grid from Lookup Table
xi = np.arange( -1, 1, 0.01 )
yi = np.arange( -1, 1, 0.01 )
gridX, gridY = np.mgrid[-1:1:0.01, -1:1:0.01]

# Initilizing Original Mask Plane from Lookup Table
maskX, maskY = np.mgrid[ -69.75:70.25:0.5, -69.75:70.25:0.5 ]    # mask plane, 0.5-mm pixels
maskPixels = np.array( [ maskX.flatten(), maskY.flatten() ] ).T

lookTab = lookTab.reshape( ( gridX.shape[0], gridX.shape[1], maskX.shape[0], maskX.shape[1] ) )

# Building System Response Matrix
tme = time.time()
for p in np.arange( len( poses ) ):
    print( 'detector pose: %i of %i' %( p + 1, len( poses ) ) )
    matrix = np.zeros( ( len( sourcePixels ), len( detectorPixels ) ), dtype = 'float32' )
    ang = poses[ p ]
    R = np.array( [  [ np.cos( ang ), 0, np.sin( ang ) ], [ 0, 1, 0 ], [ -np.sin( ang ), 0, np.cos( ang ) ] ] ).T
    B = np.array( [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ] )
    K = np.array( np.dot( R, B ) )
    sourcePixelsNew = np.dot( sourcePixels, K ) + c # Rotate and Translate Source Image Space

    for d in np.arange( len( detectorPixels ) ):
        print( 'detector pixel: %i of %i' %( d + 1, len( detectorPixels ) ) )
        weights = np.zeros( len( sourcePixels ) )
        r = np.hstack( ( sourcePixelsNew[ :, :2 ] - detectorPixels[ d ],  ( D + T + sourcePixelsNew[ :, 2 ] )[ :, np.newaxis ] ) )
        rHat = r / np.linalg.norm( r, axis = 1 )[ :, None ]
        const = ( n ** 2 * rHat[ :, 2 ] ) / ( 4 * np.pi * np.linalg.norm( r, axis = 1 ) ** 2 + 2 * n ** 2 * rHat[ :, 2 ] )
        mP = ( ( D * sourcePixelsNew[ :, :2 ] ) + ( ( T + sourcePixelsNew[ :, 2 ] )[ :, np.newaxis ] * detectorPixels[ d ] ) ) / ( D + T + sourcePixelsNew[ :, 2 ] )[ :, np.newaxis ]
        phi = np.arctan2( np.sqrt( rHat[:, 1] ** 2 +  rHat[:, 0] ** 2 ), rHat[:, 2] )
        theta = np.arctan2( rHat[:, 1],  rHat[:, 0] )

        # Finding index of the (x, y) mask plane coordinate in the far-field lookup table
        mP = np.round( ( mP + 0.25 ) * 2 ) / 2 - 0.25
        i1 = ( ( mP[ :, 0 ] - maskPixels[ 0, 0 ] ) / ( maskPixels[ 1, 1 ] - maskPixels[ 0, 1 ] ) ).astype(int)
        i2 = ( ( mP[ :, 1 ] - maskPixels[ 0, 1 ] ) / ( maskPixels[ 1, 1 ] - maskPixels[ 0, 1 ] ) ).astype(int)
        # photon is modulated by the mask
        mask = ( i1 >= 0 ) & ( i2 >= 0 ) & ( i1 < maskX.shape[ 0 ] ) & ( i2 < maskY.shape[ 0 ] )

        # Finding index of the Cartesian (theta, phi) coordinate in the far-field lookup table
        xx = np.round( np.sin( phi ) * np.cos( theta ), 2 )
        yy = np.round( np.sin( phi ) * np.sin( theta ), 2 )
        i3 = np.round( ( xx - xi[ 0 ] ) / ( xi[ 1 ] - xi[ 0 ] ) ).astype( int )
        i4 = np.round( ( yy - yi[ 0 ] ) / ( yi[ 1 ] - yi[ 0 ] ) ).astype( int )

        # Getting system response weights
        L = lookTab[ i3[mask], i4[mask], i1[mask], i2[mask] ]
        weights[ mask ] = const[ mask ] * np.exp( -mu * L )

        # Zeroing out non-modulated photons ( opposite: const[ ~mask ] )
        weights[ ~mask ] = 0
        matrix[ :, d ] = weights

    # Applying Detector Pixel Efficiency
    matrix = matrix * eff

    filename = path + name + str( int( poses[ p ] * 180 / np.pi ) ) + '.h5'
    f = tables.open_file( filename, 'w' )
    f.create_array( '/', 'matrix', matrix )
    f.close()
print( time.time() - tme )
