WIP for /r/astrophotography's #project-patchwork 

requires a python environment with:
* astroquery
* astropy
* reproject

Example conda environment:
    
    conda create -n astro -c astropy astropy astroquery reproject matplotlib


Example usage using my Horsehead and Orion images
    
    main.py --imgs_csv alistair_images.csv --astrometry_net_api_key "YOU_API_KEY_HERE"