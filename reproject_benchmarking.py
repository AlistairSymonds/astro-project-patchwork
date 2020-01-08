from reproject import reproject_adaptive, reproject_interp, reproject_exact
from reproject.utils import reproject_blocked
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import multiprocessing

def main():
    processes = [False,1,2,4]
    reproj_funcs = [reproject_adaptive, reproject_interp, reproject_interp]
    block_sizes = [500, 1000, 5000]

    mra=180.0
    mdec=45.0
    cellsize=1.0/3600
    ctype=('RA---SIN','DEC--SIN')

    # pick sizes here
    rxsize=rysize=3000
    xsize=ysize=1024

    rwcs=WCS(naxis=2)
    rwcs.wcs.ctype=ctype
    rwcs.wcs.cdelt=(-cellsize,cellsize)
    rwcs.wcs.crval=[mra,mdec]
    rwcs.wcs.crpix=[1,1]

    rheader=rwcs.to_header()
    rheader['NAXIS']=2
    rheader['NAXIS1']=rxsize
    rheader['NAXIS2']=rysize

    header=rwcs.to_header()
    header['NAXIS']=2
    header['NAXIS1']=xsize
    header['NAXIS2']=ysize
    orig_wcs = WCS(header)

    data=np.random.rand(ysize,xsize)
    for func in reproj_funcs:
        for num_p in processes:
            for bsize in block_sizes:
                print("Func is:" +str(func) + ", Processes "+ str(num_p) + ", block size " + str(bsize))
                time_start = timer()
                r, footprint = reproject_blocked(func,input_data=(data, header), output_projection=rheader,
                                                 parallel=num_p, block_size=(bsize, bsize))
                time_taken = timer() - time_start
                print("Time taken was: " + str(time_taken))


    #ax1 = plt.subplot(1,2,1, projection=WCS(rheader))
    #ax1.imshow(data, origin='lower')
    #ax1.coords['ra'].set_axislabel('Right Ascension')
    #ax1.coords['dec'].set_axislabel('Declination')
    #ax1.set_title('Input image')
    #
    #
    #ax2 = plt.subplot(1,2,2, projection=WCS(rheader))
    #ax2.imshow(r, origin='lower')
    #ax2.coords['ra'].set_axislabel('Right Ascension')
    #ax2.coords['dec'].set_axislabel('Declination')
    #ax2.coords['dec'].set_axislabel_position('r')
    #ax2.coords['dec'].set_ticklabel_position('r')
    #ax2.set_title('Reprojected image')
    #
    #plt.show()
if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()