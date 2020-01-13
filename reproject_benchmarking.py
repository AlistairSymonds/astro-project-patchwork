from reproject import reproject_adaptive, reproject_interp, reproject_exact
from reproject.utils import reproject_blocked
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import multiprocessing
import json

def main():
    processes = [False,1,2,4, True]
    reproj_funcs = [reproject_interp, reproject_adaptive] #, reproject_exact]
    block_sizes = [500, 1000, 2500]
    results = []

    mra=180.0
    mdec=45.0
    cellsize=1.0/3600
    ctype=('RA---SIN','DEC--SIN')

    # pick sizes here
    rxsize=rysize=5000
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
        reference_reprojection, reference_footprint = func(input_data=(data, header), output_projection=rheader)
        for num_p in processes:
            for bsize in block_sizes:
                time_start = timer()
                r, footprint = reproject_blocked(func,input_data=(data, header), output_projection=rheader,
                                                 parallel=num_p, block_size=(bsize, bsize), return_footprint=True)
                time_taken = timer() - time_start

                array_closenesss = np.isclose(reference_reprojection, r, equal_nan=True)
                if np.any(array_closenesss == False):
                    print("arrays not close!")
                    plt.figure()
                    plt.imshow(array_closenesss)
                    plt.show()

                print("### Time taken was: " + str(time_taken) + " for func is:" +str(func) + ", Processes "+ str(num_p) + ", block size " + str(bsize))
                results.append({'time': time_taken, 'func': str(func), 'processes': num_p, 'block_size':bsize})

    for result in results:
        print(result)

    results_file = open("benchmark_results.json", 'w+')
    json.dump(results, results_file)

    #insert analysis code here

if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()