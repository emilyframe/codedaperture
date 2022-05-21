# Data Pathway
path = '/home/raptor/repos/'

import numpy as np
import tables
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def make2DMesh( x, y, vals, vmin, vmax ):
    """generates 2D mesh plot
    """
    fig, ax = plt.subplots()
    ax.set_xlabel( 'X (mm)', fontsize = 15 )
    ax.set_ylabel( 'Y (mm)', fontsize = 15 )
    ax.tick_params( labelsize = 15 )
    if ( x is False ) and ( y is False ):
        im = ax.pcolormesh( vals.T, vmin = vmin, vmax = vmax, shading = 'gouraud' )
    else:
        im = ax.pcolormesh( x, y, vals.T, vmin = vmin, vmax = vmax, shading = 'gouraud' )
    cbar = fig.colorbar( im, format = '%.0e' )
    cbar.set_label( label = 'Intensity', rotation = 270, fontsize = 15, labelpad = 25 )
    cbar.ax.tick_params( labelsize = 15 )
    return fig, ax, im

# Opening Reconstructed Image Data File
image = path + 'gammaimage/images/image.h5'
#
X, Y, Z = np.mgrid[ -49:51:2, -69:71:2, -49:51:2 ]
x, y = np.mgrid[ -49:51:2, -69:71:2 ]
z = np.arange( -49, 51, 2 )[:]

f = tables.open_file( image, 'r' )
vals = f.root.image.read().reshape( X.shape )
f.close()

fig, ax, im = make2DMesh( x, y, vals[ :, :, 0].T, vmin = min( vals.flatten() ), vmax = max( vals.flatten() ) )
title = plt.title( 'depth = %i mm' %z[ 0 ], fontsize = 20 )

def animate( i ):
    grid = vals[ :, :, i ]
    grid = grid.reshape( x.shape[ 0 ], x.shape[ 1 ] )
    im.set_array( grid.ravel() )
    title.set_text( 'depth = %i mm' %( z[i] ) )
    return im

anim = animation.FuncAnimation( fig, animate, frames = len( z ), interval = 400 )
plt.show()
