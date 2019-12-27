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
from reproject import reproject_exact, mosaicking, reproject_interp
from astropy.wcs import WCS
from astropy import units as u

from astroquery.astrometry_net import AstrometryNet


def calc_sky_area(header):

    sky_area = None
    comments = header['comment']
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


            sky_area = header['IMAGEH'] * header['IMAGEW'] * l[0]
            print(sky_area)
    return sky_area



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
                                                      solve_timeout=1200)
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

    channel_names = ['r','g','b']
    for i in range(0, len(channel_names)):
        channel_fits = astropy.io.fits.PrimaryHDU(data=img_data[:,:,i], header=wcs_header)

        channel_fits_path_parent = downloaded_img_path.parent
        channel_fits_path_filename = downloaded_img_path.stem
        channel_fits_path = Path(str(channel_fits_path_parent / channel_fits_path_filename) +"_"+channel_names[i]+".fits")
        print(channel_fits)
        channel_fits.writeto(channel_fits_path, overwrite=True)



    return


def load_all_wcs_meta(img_dir_path):
    wcs_glob = img_dir_path.glob("**/*.wcs")

    #img_paths = []
    #wcs_headers = []
    data = []
    for w in wcs_glob:
        img_dict = {}
        fp = open(w)
        wcs_h = astropy.io.fits.header.Header()
        wcs_h = wcs_h.fromtextfile(fp)
        img_dict['wcs'] = wcs_h

        red_path = Path(str(w.parent / w.stem) + "_r.fits")
        green_path = Path(str(w.parent / w.stem) + "_g.fits")
        blue_path = Path(str(w.parent / w.stem) + "_b.fits")
        img_channels = [red_path, green_path, blue_path]

        img_dict['paths'] = (img_channels)
        img_dict['area'] = calc_sky_area(wcs_h)
        data.append(img_dict)
        fp.close()
    return data




def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--imgs_csv", help="CSV file containing links to patch work images", default="patchwork.csv")
    ap.add_argument("--arc_sec_per_px", help="An override for final arc sec/px, if no value is given the smallest "
                                             "values in all images will be used")
    ap.add_argument("--imgs_dir", help="Directory for images to be downloaded to and prepared in, defaults to the same name as CSV file used")
    ap.add_argument("--astronometry_net_api_key", "-ast_key")

    args = ap.parse_args()

    px_scale = None
    if args.arc_sec_per_px is not None:
        px_scale = float(args.arc_sec_per_px) * u.arcsec

    imgs_csv_path = Path(args.imgs_csv)

    # check args
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
    for row in csv_rows[1:]:
        url = row[0]
        name = row[1]
        author = row[2]
        prep_pool.submit(prepare_image, url, name, author, img_dir_path, args.astronometry_net_api_key)

    print("all images submitted for prep, waiting for all to be done...")
    prep_pool.shutdown()

    all_data = load_all_wcs_meta(img_dir_path)
    #sort based upon sky area,
    all_data = sorted(all_data,key=lambda x: x['area'], reverse=True)
    #open all the red fits and hope we have enough RAM I guess....
    red_hdus = []
    for img in all_data:
        red_hdus.append(astropy.io.fits.open(img['paths'][0]))
    final_wcs, final_shape = mosaicking.find_optimal_celestial_wcs(red_hdus, resolution=px_scale)
    print(final_wcs.to_header())
    print("Final pixel resolution will be: " + str(final_shape))


    ## project and align
    final_image = [] #channels will be added to this array
    channel_names = ['r', 'g', 'b']
    for c in range(0, len(channel_names)): #for RGB, once per channel

        print("#####")
        print("Reprojecting channel:" + str(c))
        print("#####")

        channel_canvas_array = np.zeros(final_shape)
        for img in all_data:
            curent_img_path = str(img['paths'][c])
            print(curent_img_path)
            patch = astropy.io.fits.open(curent_img_path)
            patch.info()


            patch[0].data = patch[0].data / 255


            array, footprint = reproject_interp(patch[0], final_wcs, shape_out=final_shape)

            channel_canvas_array = channel_canvas_array * (footprint-1)*-1
            channel_canvas_array = channel_canvas_array + np.nan_to_num(array)
            #plt.imshow(channel_canvas_array)
            #plt.show()

        final_image.append(channel_canvas_array)
        #plt.imshow(channel_canvas_array)
        #plt.show()

    patch_werked = np.array(final_image)
    patch_werked = np.moveaxis(patch_werked, 0, -1)

    plt.imshow(patch_werked)
    plt.show()
    plt.imsave('final_image.png', np.clip(patch_werked, a_min=0, a_max=1))


if __name__ == "__main__":
    multiprocessing.freeze_support()
    Image.MAX_IMAGE_PIXELS = None #keep PIL happy when opening stupidly 3x drizzled files :P
    main()