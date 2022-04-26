from pipeline import Utilities, Pipeline, Config, PeriodFinder

from datetime import datetime
import os


def main():
    """
    Main entrypoint for the pipeline process.

    Here the end-user can write config objects for each specific field and 
    then run the pipeline.

    Use one `Config` object per field.

    Things you WILL need to change per field config:

    - `Config.raw_image_dir` - set this to the directory where the raw image files are
    located for a the field.
    - `Config.image_dir` - Set this to the directory where you want the reduced image
    files to be stored. You will need a lot of space (~3GB per field).
    - `Config.image_prefix` - name of the field, which is also the prefix for each image.


    Provided is a loop or two which allows for easy running of all fields in consecutive
    order.

    Fields taken with Nebulosity and SX10 have different fits headers, so they need to 
    be processed in separate loops.

    The `field_details` variable is an array containing fundamental details of 
    each field and things that are different across each field.
    Things like the image prefix, the number of sets the field has, 
    what the flat or bias prefix is and the base location of the raw images.

    The only things you'll need to change (to run on your computer) are where the raw
    images are stored on your drive.

    Tip: `~/` is a shortcut for your home folder, even on Windows.

    """

    start_time = datetime.now()
    print("[MAIN] Started everything at {}".format(start_time.strftime("%H:%M:%S")))


    ## Change these to where you copied the JGT images to and where you want the 
    ## reduced images to be stored
    raw_base_dir = os.path.expanduser("~/mnt/jgt/")
    img_base_dir = os.path.expanduser("~/mnt/data/tmp/jgt_images")



    ## ===================================================
    ## Nebulosity fields
    
    ## Details for fields taken with Nebulosity
    ## SX10 configs are too different to be able to be put into the same array
    field_details = [
            # image_dir             image_prefix    n_sets  flat_prefix bias_prefix
            #[raw_base_dir+"/2019/1028", "l135",         6,  "flat",     "bias"],
            #[raw_base_dir+"/2021/1210Trius", "l135_5",  4,  "dflat7",   "bias2"],
            #[raw_base_dir+"/2022/0105", "l136",         9,  "dflat",    "bias2"],
            #[raw_base_dir+"/2022/0117", "l136_5",       7,  "dflat",    "bias4"],
            #[raw_base_dir+"/2022/0121", "l137_0",       7,  "dflat",    "bias_shutter"],
            #[raw_base_dir+"/2022/0124", "l137_5",       7,  "dflat",    "bias"],
            #[raw_base_dir+"/2019/0218", "l140_0",       7,  "domeflat", "BiasLast"],
            #[raw_base_dir+"/2019/0221", "l140_5",       7,  "domeflat", "Bias_end_"],
            [raw_base_dir+"/2020/0206Trius", "l141",    7,  "Flat",     "bias2"],
            #[raw_base_dir+"/2020/0212", "l141_5",       8,  "Flat",     "bias2"],
            #[raw_base_dir+"/2019/0226", "l196",         7,  "domeflat", "bias_end"],
            #[raw_base_dir+"/2019/0225", "l196_5",       6,  "domeflat", "bias_end"],
            #[raw_base_dir+"/2019/0204", "l197",         7,  "domeflat", "Bias"],
            #[raw_base_dir+"/2019/0131", "l197.5",       9,  "flat",     "biasend"],
            #[raw_base_dir+"/2019/0128", "l198",         7,  "flat",     "bias_2"],
            #[raw_base_dir+"/2019/0129", "l198_5",       7,  "flat",     "bias-1"],
            #[raw_base_dir+"/2019/0207", "l199",         5,  "domeflat", "Bias"],
    ]

    ## Run for all Nebulosity fields
    for raw_image_dir, image_prefix, n_sets, flat_prefix, bias_prefix in field_details:

        config = Config(
            raw_image_dir = os.path.expanduser(raw_image_dir),
            image_dir = os.path.expanduser(img_base_dir),
            image_prefix = image_prefix,
            n_sets = n_sets,
            flat_prefix = flat_prefix,
            bias_prefix = bias_prefix,
        )

        Pipeline.run(config)
        #Pipeline.run_analysis(config, assume_already_adjusted=True)



    ## ===================================================
    ## SX10 fields

    field_details = [
            # image_dir                 image_prefix    n_sets  flat_prefix bias_prefix
            [raw_base_dir+"/2022/0301", "l138_0",       9,      "dflat",    "bias_end"],
    ]

    for raw_image_dir, image_prefix, n_sets, flat_prefix, bias_prefix in field_details:

        config = Config(
            raw_image_dir = os.path.expanduser(raw_image_dir),
            image_dir = os.path.expanduser(img_base_dir),
            image_prefix = image_prefix,
            n_sets = n_sets,
            flat_prefix = flat_prefix,
            bias_prefix = bias_prefix,
            fits_extension = ".fits",
            fits_date_format = "%Y.%m.%dT%H:%M:%S.%f",
            has_filter_in_header = False,
        )

        Pipeline.run(config)
        #Pipeline.run_analysis(config, assume_already_adjusted=True)

    _ = Utilities.finished_job("all fields", start_time)



## TODO: Auto-copy flats for ones that need it?
def rename():
    """
    Rename certain problematic images/fields.
    Should only need to run once per machine (unless you delete your
    JGT images folder).

    Assumes they are copied straight from the JGT backup.
    Change the base path for each field and they should rename to a format
    that will play nice with the pipeline.
    Also if you make any fields that need renaming, please add them here.
    It will make future observers' lives easier.

    Ideal format:
    lxxx_y_s_iii.fit
    xxx: longitude integer part
    y: longitude decimal part, 0 or 5
    s: set number, 1..n_sets
    iii: image_number, 001..050 

    """

    print("[Rename] Starting rename")

    ## l141 first set
    base_path = os.path.expanduser("F:/2020/0206Trius")
    if os.path.exists(os.path.join(base_path, "l141_001.fit")):
        for i in range(1, 51):
            os.rename(
                    os.path.join(base_path, "l141_{:03d}.fit".format(i)),
                    os.path.join(base_path, "l141_1_{:03d}.fit".format(i))
            )

    ## l140.5 -> l140_5
    bsase_path = os.path.expanduser("F:/2019/0221")
    if os.path.exists(os.path.join(base_path, "l140.5_1_001.fit")):
        for s in range(1, 8):
            for i in range(1, 51):
                os.rename(
                        os.path.join(base_path, "l140.5_{:1d}_{:03d}.fit".format(s, i)),
                        os.path.join(base_path, "l140_5_{:1d}_{:03d}.fit".format(s, i))
                )

    ## L140 -> l140_0
    base_path = os.path.expanduser("F:/2019/0218")
    if os.path.exists(os.path.join(base_path, "L140_1_001.fit")):
        for s in range(1, 8):
            for i in range(1, 51):
                os.rename(
                        os.path.join(base_path, "L140_{:1d}_{:03d}.fit".format(s, i)),
                        os.path.join(base_path, "l140_0_{:1d}_{:03d}.fit".format(s, i))
                )

    ## Make l138_0 image numbers range from 1-50
    base_path = os.path.expanduser("F:/2022/0301")
    if not os.path.exists(os.path.join(base_path, "l138_0_1_001.fits")):
        for i in range(40, 40+50):
            os.rename(
                    os.path.join(base_path, "l138_0_1_{:03d}.fits".format(i)),
                    os.path.join(base_path, "l138_0_1_{:03d}.fits".format(i-39)),
            )

        for i in range(51, 51+50):
            os.rename(
                    os.path.join(base_path, "l138_0_3_{:03d}.fits".format(i)),
                    os.path.join(base_path, "l138_0_3_{:03d}.fits".format(i-50)),
            )

    ## l196
    base_path = os.path.expanduser("F:/2019/0226")
    if os.path.exists(os.path.join(base_path, "l196_1-1_001.fit")):
        os.rename(os.path.join(base_path, "l196_1_001.fit"), os.path.join(base_path, "test_l196_1_001.fit"))
        for i in range(1, 51):
            file_from = os.path.join(base_path, "l196_1-1_{:03d}.fit".format(i))
            file_to   = os.path.join(base_path, "l196_1_{:03d}.fit".format(i))
            os.rename(
                    file_from,
                    file_to
            )
            print("[RENAME] '{}' -> '{}'".format(file_from, file_to))

    ## l196.5
    base_path = os.path.expanduser("F:/2019/0225")
    if os.path.exists(os.path.join(base_path, "l196_5_001.fit")):
        for i in range(1, 51):
            file_from = os.path.join(base_path, "l196_5_{:03d}.fit".format(i))
            file_to   = os.path.join(base_path, "l196_5_1_{:03d}.fit".format(i))
            os.rename(
                    file_from,
                    file_to
            )
            print("[RENAME] '{}' -> '{}'".format(file_from, file_to))


    print("[Rename] Finished renaming")

if __name__ == "__main__":
    #rename()
    main()

