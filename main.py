import argparse
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
import csv
from pathlib import Path
import urllib.request
import urllib.parse
import astropy.io.fits.header
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from reproject import reproject_interp
from astropy.wcs import WCS

from astroquery.astrometry_net import AstrometryNet

def prepare_image(img_url: str, img_name: str, img_author: str, imgs_path: Path, astrometry_net_key):

    img_file_name = urllib.parse.quote_plus(img_url)
    print(img_file_name)
    img_prep_folder = imgs_path / Path(urllib.parse.quote_plus(img_name) + "_"+ urllib.parse.quote_plus(img_author)+"_" +urllib.parse.quote_plus(img_url))
    img_prep_folder.mkdir(exist_ok=True, parents=True)

    downloaded_img_path = img_prep_folder / urllib.parse.quote_plus(img_url)
    if downloaded_img_path.exists():
        print("Already downloaded file from: " + img_url + " - skipping download!")
    else:
        print("Preparing image from: " + img_url + " and downloading to: " + str(downloaded_img_path))
        urllib.request.urlretrieve(img_url, downloaded_img_path)

    wcs_file_name = str(downloaded_img_path.stem) + ".wcs"
    wcs_file_path = img_prep_folder / wcs_file_name
    if wcs_file_path.exists():
        print("Found WCS data at: " + str(wcs_file_path))
        fp = open(wcs_file_path)
        wcs_header = astropy.io.fits.header.Header()
        wcs_header = wcs_header.fromtextfile(fp)
        fp.close()
    else:
        print("Couldn't find WCS data at: " + str(wcs_file_path) + " - beginning Astrometry.net platesolving")
        ast = AstrometryNet()
        ast.api_key = astrometry_net_key

        try_again = True
        submission_id = None

        while try_again:
            try:
                if not submission_id:
                    wcs_header = ast.solve_from_image(str(downloaded_img_path),
                                                      submission_id=submission_id,
                                                      solve_timeout=600)
                else:
                    wcs_header = ast.monitor_submission(submission_id,
                                                        solve_timeout=600)
            except TimeoutError as e:
                submission_id = e.args[1]
            else:
                # got a result, so terminate
                try_again = False

        fp = open(wcs_file_path, "wb")
        wcs_header.totextfile(fp)
        fp.close()

    #combien into one fits file

    img_data = np.array(Image.open(downloaded_img_path))
    img_data = np.rollaxis(img_data,2,0)
    fits = astropy.io.fits.PrimaryHDU(data=img_data, header=wcs_header)

    fits_path = downloaded_img_path.with_suffix(".fits")
    print(fits)
    fits.writeto(fits_path, overwrite=True)



    return


def load_all_wcs_meta(img_dir_path):
    wcs_glob = img_dir_path.glob("**/*.wcs")

    img_paths = []
    wcs_headers = []
    for w in wcs_glob:
        fp = open(w)
        wcs_h = astropy.io.fits.header.Header()
        wcs_h = wcs_h.fromtextfile(fp)
        wcs_headers.append(wcs_h)
        img_paths.append(w.with_suffix(".fits"))
        fp.close()


    #wcs_headers.sort(key=)
    #find comment keyword

    sky_areas = []
    for idx, h  in enumerate(wcs_headers):
        comments = h['comment']
        for c in comments:
            if 'scale' in str(c) and 'Field' not in str(c):
                print(c)

                l = []
                for t in str(c).split():
                    try:
                        l.append(float(t))
                    except ValueError:
                        pass
                print(l)
                if len(l) != 1:
                    assert AttributeError("Multiple numbers detected in scale")


                sky_area = h['IMAGEH'] * h['IMAGEW'] * l[0]
                print(sky_area)
                patchwork_meta = {}
                patchwork_meta['path'] = img_paths[idx]
                patchwork_meta['sky_area'] = sky_area
                patchwork_meta['scale'] = l[0]
                sky_areas.append(patchwork_meta)


    print(sky_areas)

    sky_areas_sorted = sorted(sky_areas,key=lambda x: x['sky_area'], reverse=True)
    print(sky_areas_sorted)
    return sky_areas_sorted




def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--imgs_csv", help="CSV file containing links to patch work images", default="patchwork.csv")
    ap.add_argument("--arc_sec_per_px", help="An override for final arc sec/px, if no value is given the smallest "
                                             "values in all images will be used")
    ap.add_argument("--imgs_dir", help="Directory for images to be downloaded to and prepared in, defaults to the same name as CSV file used")
    ap.add_argument("--astronometry_net_api_key", "-ast_key")

    args = ap.parse_args()

    #check args

    imgs_csv_path = Path(args.imgs_csv)


    img_dir_path = None
    if args.imgs_dir is not None:
        img_dir_path = Path(args.imgs_dir)
        if img_dir_path.exists() == False or img_dir_path.exists() == False:
            assert NotADirectoryError("Couldn't find directory: "+ str(args.imgs_dir))
    else:
        img_dir_path = (imgs_csv_path.resolve().parent / imgs_csv_path.stem)

    img_dir_path.mkdir(exist_ok=True, parents=True)


    csv_rows = []
    with open(imgs_csv_path, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in csvreader:
            print(row)
            csv_rows.append(row)



    prep_pool = ProcessPoolExecutor()
    for row in csv_rows[1:]: #
        url = row[0]
        name = row[1]
        author = row[2]
        #prep_pool.submit(prepare_image, url, name, author, img_dir_path, args.astronometry_net_api_key)
        prepare_image(url, name, author, img_dir_path, args.astronometry_net_api_key)

    print("all images submitted for prep, waiting for all to be done...")
    prep_pool.shutdown()

    imgs = load_all_wcs_meta(img_dir_path)

    smallest_scale = args.arc_sec_per_px
    if smallest_scale is None:
        smallest_scale = imgs[0]['scale']
        for i in imgs:
            if i['scale'] < smallest_scale:
                smallest_scale = i['scale']

    print("Using image scale of: " + str(smallest_scale) + " arc sec per pixel")


    ## project and align
    base_file_path = imgs[0]['path']
    base_hdu = astropy.io.fits.open(base_file_path)[0]
    for patch_tuple in imgs[1:]:
        patch = astropy.io.fits.open(str(patch_tuple['path']))
        print(base_hdu.header)
        wcs = WCS(base_hdu.header, naxis=2) #this is the line that has issues
        print(wcs)
        array, footprint = reproject_interp(patch, base_hdu.header)


        ax1 = plt.subplot(1, 2, 1, projection=WCS(base_hdu.header, naxis=0))
        ax1.imshow(array, origin='lower', vmin=-2.e-4, vmax=5.e-4)
        ax1.coords.grid(color='white')
        ax1.coords['ra'].set_axislabel('Right Ascension')
        ax1.coords['dec'].set_axislabel('Declination')
        ax1.set_title('Reprojected MSX band E image')

        ax2 = plt.subplot(1, 2, 2, projection=WCS(base_hdu.header))
        ax2.imshow(footprint, origin='lower', vmin=0, vmax=1.5)
        ax2.coords.grid(color='white')
        ax1.coords['ra'].set_axislabel('Right Ascension')
        ax1.coords['dec'].set_axislabel('Declination')
        ax2.coords['dec'].set_axislabel_position('r')
        ax2.coords['dec'].set_ticklabel_position('r')
        ax2.set_title('MSX band E image footprint')

        plt.show()

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()