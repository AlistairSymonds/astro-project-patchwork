import argparse
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
import csv
from pathlib import Path
import urllib.request
import urllib.parse
import astropy.io.fits.header

from astroquery.astrometry_net import AstrometryNet

def prepare_image(img_url: str, img_name: str, img_author: str, imgs_path: Path, astrometry_net_key):

    img_file_name = urllib.parse.quote_plus(img_url)
    print(img_file_name)
    img_prep_folder = imgs_path / Path(urllib.parse.quote_plus(img_name) + "_"+ urllib.parse.quote_plus(img_author)+"_" +urllib.parse.quote_plus(img_url))
    img_prep_folder.mkdir(exist_ok=True,parents=True)

    downloaded_file_path = img_prep_folder / urllib.parse.quote_plus(img_url)
    if downloaded_file_path.exists():
        print("Already downloaded file from: " + img_url + " - skipping download!")
    else:
        print("Preparing image from: " + img_url + " and downloading to: " + str(downloaded_file_path))
        urllib.request.urlretrieve(img_url, downloaded_file_path)




    wcs_file_name = str(downloaded_file_path.stem) + ".wcs"
    wcs_file_path = img_prep_folder / wcs_file_name
    if wcs_file_path.exists():
        print("Found WCS data at: " + str(wcs_file_path))
        fp = open(wcs_file_path)
        wcs_header = astropy.io.fits.header.Header()
        wcs_header.fromtextfile(fp)
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
                    wcs_header = ast.solve_from_image(str(downloaded_file_path),
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










    return "xd"





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
        #prep_pool.submit(prepare_image, url, name, author, img_dir_path)
        prepare_image(url, name, author, img_dir_path, args.astronometry_net_api_key)

    print("all images submitted for prep, waiting for all to be done...")
    prep_pool.shutdown()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()