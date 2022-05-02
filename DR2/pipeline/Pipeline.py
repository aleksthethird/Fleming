from pipeline import Config, Reducer, Cataloguer, ShiftFinder, FluxFinder
from pipeline import DataAnalyser, Utilities, VariableDetector

from datetime import datetime
import os
import numpy as np

HOME = os.path.expanduser("~")

def run(config, show_plots=False, show_errors=False, solve_astrometry=True, skip_existing_images=True):
    """
    Run the whole pipeline for a single field.

    Parameters
    ----------

    config: Config
        Config objecet for field

    show_plots: bool, optional
        Show plots to the output (only useful in ipython notebooks).

    show_errors: bool, optional
        Plot error bars on the light curves

    solve_astrometry: bool, optional
        Flag to skip the astrometry plate solve process. Setting to false will not give RA/DEC
        coordinates of the stars.

    skip_existing_images: bool, optional
        Flag to skip reducing images which have already been reduced.

    """

    start_time = datetime.now()
    print("[JOB] Started {} at {}".format(config.image_prefix, start_time.strftime("%H:%M:%S")))
    
    n_sources, time = images_to_light_curves(config, start_time,
            skip_existing_images, solve_astrometry=solve_astrometry)

    _results_table, output_time = run_existing(config, n_sources, time,
            show_plots=show_plots, show_errors=show_errors)

    
    _ = Utilities.finished_job("everything for {}".format(config.image_prefix), start_time)


def run_analysis(config, show_plots=False, show_errors=False, assume_already_adjusted=True):
    """
    Only run the variable star selection part of the pipeline.

    Parameters
    ----------

    config: Config
        Config objecet for field

    show_plots: bool, optional
        Show plots to the output (only useful in ipython notebooks).

    show_errors: bool, optional
        Plot error bars on the light curves

    assume_already_adjusted: bool
        Skip adjusting the raw light curves if we already have done in previous runs.


    Returns
    -------

    results_table: numpy array
        Table of all of the stars deemed variable (also written to disk).

    end_time: datetime
        datetime object of when the function finished.
    
    """

    start_time = datetime.now()

    c = Cataloguer(config)

    catalogue_set_number = 1
    catalogue_image_number = 1

    catalogue_image_path = os.path.join(config.image_dir,
            config.image_format_str
            .format(catalogue_set_number, catalogue_image_number))

    print("[Pipeline] Cataloguing image {}".format(catalogue_image_path))
    n_sources = c.generate_catalogue(catalogue_image_path, solve=False, write_file=False)

    cataloguer_time = Utilities.finished_job("cataloguing stars", start_time)

    run_existing(config, n_sources, cataloguer_time,
            show_plots=show_plots, show_errors=show_errors,
            assume_already_adjusted=assume_already_adjusted)

    _ = Utilities.finished_job("everything for {}".format(config.image_prefix), start_time)


## =============================================================================================
## The above two functions are all you need if you just want results.
## The next functions give finer control over which parts of the pipeline to run, but are not
## necessary

def images_to_light_curves(config, start_time, skip_existing_images=True, solve_astrometry=True):
    """
    Takes raw images and creates un-adjusted light curves.

    Produces catalogue, shifts and times of the images.

    Parameters
    ----------

    config: Config
        Config object of the field to analyse.

    start_time: datetime 
        datetime object of the start time of the pipeline (or field).

    skip_existing_images: bool, optional
        Flag to skip reducing images which have already been reduced.

    solve_astrometry: bool, optional
        Flag to skip the astrometry plate solve process. Setting to false will not give RA/DEC
        coordinates of the stars.


    Returns
    -------

    
    n_sources: int
        Number of sources in the catalogue

    light_curve_time: datetime
        datetime object of when the process finished.

    """

    ## Reducer
    ## Takes raw images, subtracts bias and divides by flat field

    r = Reducer(config, "No filter")        ## Only "No filter" for Trius
    r.reduce(skip_existing=skip_existing_images)    ## Skip images that have already been reduced

    reducer_time = Utilities.finished_job("reducing images", start_time)
    
    ## Cataloguer
    ## Creates a catalogue of stars found in the given image
    c = Cataloguer(config)

    catalogue_image_path = os.path.join(config.image_dir,
            config.image_format_str
            .format(config.catalogue_set_number, config.catalogue_image_number))

    print("[Pipeline] Cataloguing image {}".format(catalogue_image_path))
    n_sources = c.generate_catalogue(catalogue_image_path, solve=solve_astrometry)
    c.generate_image_times()

    cataloguer_time = Utilities.finished_job("cataloguing stars", reducer_time)
    

    ## Moving object finder
    #mof = MovingObjectFinder.MovingObjectFinder(Constants.folder)
    #mof.find_moving_objects()
     
     
    ## ShiftFinder
    ## Gets the shift of each star for each image in the series
    print("[Pipeline] Finding shifts in each image")
    sf = ShiftFinder(config, n_sources)
    sf.generate_shifts()
    reference_ids = sf.get_reference_ids()
     
    shift_finder_time = Utilities.finished_job("finding shifts", cataloguer_time)
     

    ## FluxFinder
    ff = FluxFinder(config, n_sources)
    print("[Pipeline] Making initial light curves")

    ## Find the flux of each star in each image then create a light curve
    ## Write the light curves to file
    ff.make_light_curves()

    light_curve_time = Utilities.finished_job("making light curves", shift_finder_time)

    return n_sources, light_curve_time


def run_existing(config, n_sources, start_time,
        show_plots=False, show_errors=False, assume_already_adjusted=True):
    """
    Run for a field and assume the light curves are already present.
    Very useful for fast debugging the later stages of the pipeline like the variable finder etc.

    Parameters
    ----------

    config: Config
        Config objecet for field

    n_sources: int
        Number of stars in the field

    start_time: datetime
        datetime object of the beginning of the function call

    show_plots: bool, optional
        Show plots to the output (only useful in ipython notebooks).

    show_errors: bool, optional
        Plot error bars on the light curves

    assume_already_adjusted: bool
        Skip adjusting the raw light curves if we already have done in previous runs.


    Returns
    -------

    results_table: numpy array
        Table of all of the stars deemed variable (also written to disk).

    end_time: datetime
        datetime object of when the function finished.
    
    """

    ## FluxFinder
    ff = FluxFinder(config, n_sources)
    
    ## Assume already have unadjusted light curves

    ## Sometimes we want to remake adjusted light curves
    if not assume_already_adjusted:

        print("[Pipeline] Creating average light curve")
        da = DataAnalyser(config, adjusted=False)
        mean, std, med, n_positive = da.get_means_and_stds()
        source_ids = da.get_source_ids()
        da.plot_means_and_stds()
        
        vd = VariableDetector(config, source_ids, mean, std, med, n_positive, adjusted=False)
        exclude_ids = vd.std_dev_search(config.avg_exclude_threshold, 1e4, 0, 0)
        avg_ids = da.get_ids_for_avg(exclude_ids)

        da.make_avg_curve(avg_ids)
        make_avg_curve_time = Utilities.finished_job("making average curve", start_time)

        ## 'adjusts' light curves by dividing by average
        print("[Pipeline] Adjusting")
        ff.create_adjusted_light_curves(source_ids, std)

        adjustment_time = Utilities.finished_job("adjusting light curves", make_avg_curve_time)
    else:
        adjustment_time = start_time



    ## Now we have adjusted light curves
    ## New DataAnalyser for non-adjusted
    da = DataAnalyser(config, adjusted=True)
    means_adj, stds_adj, medians_adj, n_positive_adj = da.get_means_and_stds()
    source_ids = da.get_source_ids()
    da.plot_means_and_stds()

    print("[Pipeline] Getting variables post-adjustment")
    vd = VariableDetector(config, source_ids, means_adj, stds_adj,
            medians_adj, n_positive_adj, adjusted=True)

    n_measures = config.n_sets * config.set_size
    variable_ids_s = vd.std_dev_search(config.variability_threshold,
            config.variability_max, config.min_signal_to_noise, n_measures)
    variable_ids_a = vd.amplitude_search(config.amplitude_score_threshold)
    variable_ids_a = vd.filter_variables(variable_ids_a)

    variable_ids = np.unique(np.concatenate((variable_ids_a, variable_ids_s)))
    #variable_ids = vd.filter_variables(variable_ids)

    post_time = Utilities.finished_job("post-adjustment", adjustment_time)

    print("[Pipeline] Plotting variable curves")
    ff.plot_avg_light_curve(config.avg_curve_path, show=show_plots, show_errors=show_errors)
    ff.plot_given_light_curves(variable_ids, adjusted=True, show=show_plots, show_errors=show_errors)
    #ff.plot_adjusted_comparison(variable_ids, plot_dir=config.testing_dir,
    #        show=show_plots, show_errors=show_errors)

    print("[Pipeline] Outputting results")
    results_table = da.output_results(variable_ids, vd)
    #results_table = da.output_results(source_ids, vd)
    ff.create_thumbnails(results_table)
    print("[Pipeline] Variables found: total {}; std {}; amp {}"
            .format(len(variable_ids), len(variable_ids_s), len(variable_ids_a)))

    output_time = Utilities.finished_job("outputting results", post_time)
    
    return results_table, output_time


