import os


## TODO: Append image_prefix to times and shifts files
class Config:
    """
    Config object for each run.
    
    Set parameters here which won't change for a single run.
    Pass this object around to the other classes for a nicer
    way of having 'global' variables.

    Note: Anything with "_dir" suffix is an absolute path to a directory.
          Anything with "_subdir" is a subdirectory of the workspace.
          "_fname" is a file name (without extension).
          "_path" is an absolute path to a file (not directory).

    Don't change the things in here, these are just defaults.
    If you want to change something for a run or user, do it 
    when making a new config object like:

    ```
    config = Config(
                pipeline_root = "/where/i/am/working",
                image_dir = "/path/to/reduced/images",
                has_sets = True,
                n_sets = 8,
                set_size = 50,
             )
    ```

    Any directories listed here will be created if they don't already exist.

    Hopefully most things are self-explanatory.
    Further explanations are below:

    Parameters
    ----------

    """
    
    ## Defaults for the config.
    ## These shouldn't need changing here
    ## Set these when creating the constructor in Main.py
    def __init__(self,
            pipeline_root       = ".",
            workspace_subdir    = "workspace",
            image_subdir        = "../images",
            image_dir           = "",
            raw_image_subdir    = "../raw_images",
            raw_image_dir       = "",
            flux_subdir         = "fluxes",
            light_curve_subdir  = "light_curves",
            adjusted_curve_subdir = "adjusted_light_curves",
            output_subdir       = "results",
            streak_subdir       = "streaks",
            moving_obj_subdir   = "moving_objects",
            rejects_subdir      = "rejects",
            periods_subdir      = "periods",
            testing_subdir      = "testing",
            hand_picked_subdir  = "hand_picked",

            id_map_fname        = "id_mapper",
            moving_obj_fname    = "moving_obj",
            astrometry_job_fname = "astrometryjob",
            avg_curve_fname     = "avg",

            reduced_prefix      = "r_",
            catalogue_prefix    = "catalogue_",
            time_prefix         = "times_",
            shift_prefix        = "shift_",
            flux_prefix         = "flux_",
            bias_prefix         = "bias",
            flat_prefix         = "dflat",

            fits_extension      = ".fit",
            standard_file_extension = ".txt",
            plot_file_extension = ".jpg",

            line_ending         = "\n",
            identifier          = "id",
            table_format        = "ascii",
            fits_date_format    = "%Y-%m-%dT%H:%M:%S",

            avg_exclude_threshold   = 0.4,
            variability_threshold   = 1.2,
            variability_max         = 5,
            amplitude_score_threshold = 0.85,
            check_radius            = 5,
            cosmic_threshold        = 5,
            min_flux_cutoff         = 1.5,
            max_flux_cutoff         = 100,
            edge_limit              = 50,
            inner_radius            = 8,
            outer_radius            = 13,
            moving_obj_check_image  = 30,
            n_reference_stars       = 10,
            min_signal_to_noise     = 2,
            n_clip_iterations       = 3,
            counts_clip_threshold   = 3,

            period_chi2_range       = 3,
            period_max_iterations   = 20,
            period_width_adjustment = 1.5,
            n_sample_periods        = 500,
            plot_fit_chi2           = False,
            plot_fit_comparison     = False,

            image_width         = 2432,
            image_height        = 1616,
            image_prefix        = "l137_0",
            has_sets            = True,
            set_size            = 50,
            n_sets              = 7,
            has_filter_in_header= True,
            catalogue_set_number= 1,
            catalogue_image_number=1,

            astrometry_timeout  = 1200,
            api_key             = "",
            api_key_file        = "astrometry_api_key.txt"
            ):

        ## Set all the values to the object

        self.pipeline_root          = pipeline_root
        self.workspace_dir          = os.path.join(pipeline_root, workspace_subdir)

        ## Only overwrite image_dir if the user didn't give us an explicit path
        if len(image_dir) > 0:
            self.image_dir = image_dir
        else:
            self.image_dir = os.path.join(self.workspace_dir, image_subdir)

        if len(raw_image_dir) > 0:
            self.raw_image_dir = raw_image_dir
        else:
            self.raw_image_dir = os.path.join(self.workspace_dir, raw_image_subdir)

        self.flux_dir               = os.path.join(self.workspace_dir, flux_subdir)
        self.light_curve_dir        = os.path.join(self.workspace_dir, light_curve_subdir)
        self.adjusted_curve_dir     = os.path.join(self.workspace_dir, adjusted_curve_subdir)
        self.output_dir             = os.path.join(self.workspace_dir, output_subdir)
        self.streak_dir             = os.path.join(self.workspace_dir, streak_subdir)
        self.moving_obj_dir         = os.path.join(self.workspace_dir, moving_obj_subdir)
        self.rejects_dir            = os.path.join(self.workspace_dir, rejects_subdir)
        self.periods_dir            = os.path.join(self.workspace_dir, periods_subdir)
        self.testing_dir            = os.path.join(self.workspace_dir, testing_subdir)
        self.hand_picked_dir        = os.path.join(self.workspace_dir, hand_picked_subdir)

        fname = "{}{}{}".format(time_prefix, image_prefix, standard_file_extension)
        self.time_path              = os.path.join(self.workspace_dir, fname)
        fname = "{}{}{}".format(shift_prefix, image_prefix, standard_file_extension)
        self.shift_path             = os.path.join(self.workspace_dir, fname)
        fname = "{}{}".format(id_map_fname, standard_file_extension)
        self.id_map_path            = os.path.join(self.workspace_dir, fname)
        fname = "{}{}".format(moving_obj_fname, standard_file_extension)
        self.moving_obj_path        = os.path.join(self.workspace_dir, fname)
        fname = "{}{}".format(astrometry_job_fname, standard_file_extension)
        self.astrometry_job_path    = os.path.join(self.workspace_dir, fname)
        fname = "{}_{}{}".format(image_prefix, avg_curve_fname, standard_file_extension)
        self.avg_curve_path         = os.path.join(self.workspace_dir, fname)

        self.reduced_prefix         = reduced_prefix
        self.catalogue_prefix       = catalogue_prefix
        self.flux_prefix            = flux_prefix
        self.bias_prefix            = bias_prefix
        self.flat_prefix            = flat_prefix

        self.fits_extension             = fits_extension
        self.standard_file_extension    = standard_file_extension
        self.plot_file_extension        = plot_file_extension

        self.line_ending        = line_ending
        self.identifier         = identifier
        self.table_format       = table_format
        self.fits_date_format   = fits_date_format

        self.avg_exclude_threshold  = avg_exclude_threshold
        self.variability_threshold  = variability_threshold
        self.variability_max        = variability_max
        self.amplitude_score_threshold = amplitude_score_threshold
        self.check_radius           = check_radius
        self.cosmic_threshold       = cosmic_threshold
        self.min_flux_cutoff        = min_flux_cutoff
        self.max_flux_cutoff        = max_flux_cutoff
        self.edge_limit             = edge_limit
        self.inner_radius           = inner_radius
        self.outer_radius           = outer_radius
        self.moving_obj_check_image = moving_obj_check_image
        self.n_reference_stars      = n_reference_stars
        self.min_signal_to_noise    = min_signal_to_noise
        self.n_clip_iterations      = n_clip_iterations
        self.counts_clip_threshold  = counts_clip_threshold

        self.plot_fit_chi2           = plot_fit_chi2
        self.plot_fit_comparison     = plot_fit_comparison
        self.period_chi2_range       = period_chi2_range
        self.period_max_iterations   = period_max_iterations  
        self.period_width_adjustment = period_width_adjustment
        self.n_sample_periods        = n_sample_periods

        self.image_width        = image_width
        self.image_height       = image_height
        self.image_prefix       = image_prefix
        self.has_sets           = has_sets
        self.set_size           = set_size
        self.n_sets             = n_sets
        self.has_filter_in_header= has_filter_in_header
        self.catalogue_set_number = catalogue_set_number
        self.catalogue_image_number = catalogue_image_number

        self.astrometry_timeout = astrometry_timeout
        self.api_key_file       = api_key_file
        self.api_key            = api_key


        ## ===============
        ## Setting up other stuff

        ## numpy data type for light curve table
        self.light_curve_dtype = [
            ('time', 'float64'),
            ('counts', 'float64'),
            ('counts_err', 'float64')]

        ## Read in astrometry api key from file
        with open(os.path.join(pipeline_root, api_key_file)) as f:
            self.api_key = f.read().replace(line_ending, "")

        ## Create catalogue path
        fname = "{}{}{}".format(catalogue_prefix, image_prefix, standard_file_extension)
        self.catalogue_path = os.path.join(self.workspace_dir, fname)

        ## Format strings used to identify image file names and 
        ## light curve file names
        if self.has_sets:
            self.image_format_str = "{}{}_{}_{}{}".format(
                    reduced_prefix, image_prefix, "{:1}", "{:03}", fits_extension)
            self.raw_image_format_str = "{}_{}_{}{}".format(
                    image_prefix, "{:1}", "{:03}", fits_extension)
        else:
            self.image_format_str = "{}{}_{}{}".format(
                    reduced_prefix, image_prefix, "{:04}", fits_extension)
            self.raw_image_format_str = "{}_{}{}".format(
                    image_prefix, "{:04}", fits_extension)

        self.source_format_str      = "{}_{}{}{}".format(
                image_prefix, identifier, "{:04}", standard_file_extension)


        ## ===============
        ## Create directories if they don't exist

        if not os.path.exists(self.workspace_dir):
            os.mkdir(self.workspace_dir)

        if not os.path.exists(self.image_dir):
            os.mkdir(self.image_dir)

        if not os.path.exists(self.raw_image_dir):
            print("[Config] WARNING: Raw image dir doesn't exist: {}".format(self.raw_image_dir))
            print("........ are you sure you set the right variable?")
            os.mkdir(self.raw_image_dir)

        if not os.path.exists(self.flux_dir):
            os.mkdir(self.flux_dir)

        if not os.path.exists(self.light_curve_dir):
            os.mkdir(self.light_curve_dir)

        if not os.path.exists(self.adjusted_curve_dir):
            os.mkdir(self.adjusted_curve_dir)

        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

        if not os.path.exists(self.streak_dir):
            os.mkdir(self.streak_dir)

        if not os.path.exists(self.moving_obj_dir):
            os.mkdir(self.moving_obj_dir)

        if not os.path.exists(self.rejects_dir):
            os.mkdir(self.rejects_dir)

        if not os.path.exists(self.periods_dir):
            os.mkdir(self.periods_dir)

        if not os.path.exists(self.testing_dir):
            os.mkdir(self.testing_dir)

        if not os.path.exists(self.hand_picked_dir):
            os.mkdir(self.hand_picked_dir)

