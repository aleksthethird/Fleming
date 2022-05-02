from astropy.table import Table
from astropy.visualization import (ZScaleInterval, LinearStretch, ImageNormalize)
from PIL import Image
from photutils import aperture_photometry
from photutils import CircularAperture
from photutils import CircularAnnulus
from astropy.io import fits

import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

from . import Utilities, Config

class FluxFinder:
    """
    Flux Finder class.

    Anything to do with creating light curves and getting 
    statistics from the images.

    Only one `FluxFinder` object per field is needed.

    Attributes
    ----------

    config: Config
        Config object for the field

    x_shifts, y_shifts: numpy array
        Pixel shifts of each star since first image

    catalogue: numpy array
        DAO star catalogue containing IDs, fluxes etc of each source.

    n_sources: int
        Number of sources the object knows about
    
    n_measures: int
        Number of total images for the field. (n_sets * set_size)

    """

    config  = None
    
    #list of all total shifts
    x_shifts = None
    y_shifts = None
    
    #catalogue of image 1
    catalogue = None
    n_sources = None
    n_measures = None
    
    
    def __init__(self, config, n_sources):
        """
        Initialises class and reads shifts and catalogue from files.

        Parameters
        ----------
        
        config: Config
            Config object for the field

        n_sources: int
            Number of sources the object knows about

        """

        self.config = config
        self.n_sources = n_sources
        self.n_measures = self.config.set_size*self.config.n_sets
        
        self.catalogue = Utilities.read_catalogue(self.config)
        self.x_shifts, self.y_shifts = np.genfromtxt(self.config.shift_path).transpose()

        ## TODO: Sanity check if catalogue has same number of sources
        ## as arg passed



    def get_total_shift(self, set_number=0, image_number=0):
        """
        Get the total shift between the specified image and the first image.
        Requires `has_sets=True`

        Parameters
        ----------

        set_number: int
            Set of the image

        image_number: int
            Number of image in set

        """
    
        #file is just one long list of shifts, with a shift associated with 
        #each image - below is the expression required
        #to get the shift for the image belong to set number 'set' and image 
        #number 'image_number'
        index = (set_number-1)*self.config.set_size + image_number-1
        
        total_x = self.x_shifts[index]
        total_y = self.y_shifts[index]

        return total_x, total_y
    
    

    def find_fluxes(self, set_number=0, image_number=0):
        """
        Find the fluxes of all stars in a given image.

        Sums up the counts in a circle around a center.
        Subtracts background via an annulus.
        Backgrounds can be overestimated if another star is too close.

        Estimates uncertainties using Poisson noise models from
        2nd year observational techniques.
        Uncertainties aren't concrete, but it's something.

        Parameters
        ----------

        set_number: int
            Set of the image

        image_number: int
            Number of image in set

        Returns
        -------

        phot_table2: Table
            Table containing the fluxes and uncertainties of all sources.

        """

        print("[FluxFinder] Finding fluxes in image: set {:1}; image: {:03}"
                .format(set_number, image_number))

        #get total shift from first image to this one
        x_shift, y_shift = self.get_total_shift(set_number=set_number, image_number=image_number)
        
        ## TODO: Assumes catalogue is first image here
        #add the shift onto the positions of the stars in the first image
        #to find their positions in this image
        x = self.catalogue['xcentroid'] + x_shift
        y = self.catalogue['ycentroid'] + y_shift
        

        positions = np.array([x, y]).transpose()
        
        image_data = Utilities.get_image(self.config, set_number, image_number)[0].data
        x_max = self.config.image_width
        y_max = self.config.image_height

        ## Local background subtraction

        ## Define size of background aperture
        background_annuli = CircularAnnulus(
                positions,
                r_in=self.config.inner_radius,
                r_out=self.config.outer_radius)

        ## Star apertures
        ## TODO: Variable aperture sizes?
        star_apertures = CircularAperture(positions, r=self.config.inner_radius-1)

        all_apertures = [star_apertures, background_annuli]

        # find counts in each aperture and annulus
        phot_table2 = aperture_photometry(image_data, all_apertures)

        ## Probably not needed, will almost always be the same
        n_sources_phot = len(phot_table2['id'])

        if n_sources_phot != self.n_sources:
            print("[FluxFinder] Warn: Aperture photometry found {} sources, self has {}"
                    .format(n_sources_phot, self.n_sources))

        for col in phot_table2.colnames:
            phot_table2[col].info.format = '%.8g'  # for consistent table output
    

        ## Set new columns to some unphysical value
        phot_table2['mean'] = phot_table2['aperture_sum_1']*0 - 1e9
        phot_table2['median'] = phot_table2['aperture_sum_1']*0 - 1e9
        phot_table2['residual_aperture_sum_mean'] = phot_table2['aperture_sum_1']*0 - 1e9
        phot_table2['residual_aperture_sum_med'] = phot_table2['aperture_sum_1']*0 - 1e9
        phot_table2['counts_err'] = phot_table2['aperture_sum_1']*0 + 1e9

        for i in range(n_sources_phot):
            
            x, y = positions[i]

            ## Check that largest aperture does not exceed the boundaries of the image
            if Utilities.is_within_boundaries(x, y, x_max, y_max, self.config.outer_radius):

                # Calc the mean background in the second aperture ring
                bkg_mean = phot_table2['aperture_sum_1'][i] / background_annuli[i].area
                phot_table2['mean'][i] = bkg_mean
                
                # Calc background level in each main aperture and subtract
                bkg_sum = bkg_mean * star_apertures[i].area
                final_sum = phot_table2['aperture_sum_0'][i] - bkg_sum
                phot_table2['residual_aperture_sum_mean'][i] = final_sum

                # Calc background median
                ann_mask = background_annuli[i].to_mask(method='center')
                weighted_data = ann_mask.multiply(image_data)
                phot_table2['median'][i] = np.median(weighted_data[weighted_data != 0])

                # Calc 
                bkg_med = phot_table2['median'][i]
                bkg_sum = bkg_med * star_apertures[i].area
                final_sum = phot_table2['aperture_sum_0'][i] - bkg_sum
                phot_table2['residual_aperture_sum_med'][i] = final_sum

                # Calc errors (using method from 2nd year ObsTech lectures)
                star_var = phot_table2['aperture_sum_0'] 
                bkg_var  = (star_apertures[i].area/background_annuli[i].area)**2 * phot_table2['aperture_sum_0'][i]
                counts_err = np.sqrt(star_var + bkg_var)
                phot_table2['counts_err'] = counts_err

            ## If aperture doesn't fit, don't bother with the star.
            ## If it's too close to the edge of the chip, it's probably not the best data anyway

            #else:
            #    print("[FluxFinder] Source {} is not within boundary {},{}; {},{}"
            #            .format(i, x_max, y_max, x, y))

        phot_table2['residual_aperture_sum_mean'].info.format = '%.8g'  # for consistent table output
        phot_table2['residual_aperture_sum_med'].info.format = '%.8g'  # for consistent table output
                
        return phot_table2
        


    ## TODO: Breaks if no sets
    def make_light_curves(self):
        """
        Turns a series of images into a series of light curves.
        Writes light curve for each source to a file.

        """

        if self.config.has_sets == False:
            print("[FluxFinder] Error: Cannot make light curves, not implemented for has_sets=False")
            exit()
        
        times = [line.rstrip(self.config.line_ending) for line in open(self.config.time_path)]
        n_times = len(times)

        ## Datetime object of the start of the observation
        dt_obs_start = datetime.strptime(times[0], self.config.fits_date_format)

        ## Light curves
        ## Dim 0: each source
        ## Dim 1: time/counts/error
        ## Dim 2: measurement number
        light_curves = np.zeros((self.n_sources, 3, n_times))
          
        ## iterate over each image
        ## image_index iterates over each image/time/measurement
        ## j iterates over source index
        image_index = 0
        for _, s, i in Utilities.loop_images(self.config):
            #print("[FluxFinder] Making light curve for image: set {:1}; image: {:03}".format(s, i))
            
            source_table = self.find_fluxes(s, i)

            ## TODO: Decide if using mean or median
            #counts = source_table['residual_aperture_sum_median'][j]
            counts = source_table['residual_aperture_sum_mean']
                            
            ## Iterate through each source in the table
            for j in range(self.n_sources):
                    
                #get time matching the image
                date_and_time = times[image_index]

                ## Use datetime objects to calc time difference
                dt_image = datetime.strptime(date_and_time, self.config.fits_date_format)
                dt_elapsed = dt_image - dt_obs_start
                seconds_elapsed = dt_elapsed.seconds
                


                ## If counts has reasonable value
                if counts[j] > 0: 
                    light_curves[j][1][image_index] = float(counts[j])
                    light_curves[j][2][image_index] = float(source_table['counts_err'][j])
                else:
                    ## If not (probably if we failed to find fluxes, 
                    ## or if we ran off the edge of the chip)
                    ## Set to some unphysical error values.
                    light_curves[j][1][image_index] = -1.0
                    light_curves[j][2][image_index] = 1e9

                ## Include time, regardless of value
                light_curves[j][0][image_index] = seconds_elapsed

            image_index += 1

        #loop through all light curves, writing them out to a file
        for j in range(self.n_sources):
            #build light curve path
            fname = self.config.source_format_str.format(self.catalogue['id'][j])
            path = os.path.join(self.config.light_curve_dir, fname)

            ## Write light curve
            ## Transpose to have columns
            ## Only one light curve per file
            np.savetxt(path, light_curves[j].transpose())
    
            
        
        
        
    def plot_light_curve(self, source_id=None, curve_path=None, plot_dir=None,
            adjusted=False, show=False, close=True, show_errors=True):
        """
        Plot light curve of star with the given ID from catalogue

        Parameters
        ----------

        source_id: int
            ID of the source to plot light curve for.
            Plots average light curve if `None`

        curve_path: string, optional
            Path to save the image to

        plot_dir: string, optional
            Directory to save the finished plots to

        adjusted: bool, optional
            Has the light curve been divided by the average flux?

        show: bool, optional
            Should we show the plot?

        close: bool, optional
            Should we close the plot? (useful for overplotting)

        show_errors: bool, optional
            Should we plot error bars?


        Returns
        -------

        times: numpy array
            x-axis of the plot

        """

        if source_id == None:
            print("[FluxFinder] Plotting average light curve")
        else:
            print("[FluxFinder] Plotting light curve for source {:04} (adjusted={})"
                .format(source_id, adjusted))

        
        ## If we weren't given a path, make one
        if curve_path == None:
            fname = self.config.source_format_str.format(source_id)
            if adjusted:
                curve_path = os.path.join(self.config.adjusted_curve_dir, fname)
            else:
                curve_path = os.path.join(self.config.light_curve_dir, fname)

        ## Get the light curve data
        curve = np.genfromtxt(curve_path, dtype=self.config.light_curve_dtype).transpose()
        
        times = curve['time']/3600
        fluxes = curve['counts']
        err = curve['counts_err']
        
        ## Normalise to ~1
        #mean = np.mean(fluxes)
        #normalised_fluxes = fluxes / mean
        median = np.median(fluxes)
        normalised_fluxes = fluxes/median
        normalised_err    = err/median
        
        minimum = min(normalised_fluxes)
        maximum = max(normalised_fluxes)
    
        ## TODO: Multiple axes, seconds and minutes?
        if show_errors:
            plt.errorbar(times, normalised_fluxes, yerr=normalised_err,
                    fmt='.', elinewidth=0.7, ecolor="gray", capsize=3)
        else:
            plt.scatter(times, normalised_fluxes,
                    marker='.')

        plt.xlabel("Time [hours]")
        plt.ylabel("Relative flux [counts/mean]")
        

        ## If we don't have an ID, assume we're looking at the average
        if source_id == None:
            fname = "LC_{}_avg.jpg".format(self.config.image_prefix)
            plt.title("Average light curve of bright sources in {} (adjusted={})"
                .format(self.config.image_prefix, adjusted))

        elif source_id >= 0:
            fname = "LC_{}_{}{:04}.jpg".format(
                self.config.image_prefix, self.config.identifier, source_id)
            plt.title("Light curve for source {:04} in {} (adjusted={})"
                .format(source_id, self.config.image_prefix, adjusted))

        else:
            print("[FluxFinder] Error: Cannot plot light curve, source id '{}' invalid"
                    .format(source_id))
            plt.title("Light curve for unknown source (id {}) in {} (adjusted={})"
                .format(source_id, self.config.image_prefix, adjusted))
            
        ## Should we show the image to the user? (useful for ipython notebooks)
        if show:
            plt.show()

        ## If we don't have a place to save, default to the output directory
        if plot_dir == None:
            image_file = os.path.join(self.config.output_dir, fname)
        else:
            image_file = os.path.join(plot_dir, fname)
        plt.savefig(image_file)

        ## If we should close the plot.
        ## Usually yes, but can be used to overplot things.
        if close:
            plt.close()

        ## Return the x-axis for things to be overplot
        return times


        
    def plot_given_light_curves(self, ids, plot_dir=None, adjusted=False, show=False, show_errors=False):
        """
        Give a list of source ids and plot the light curve for all of them.

        Parameters
        ----------

        ids: numpy array
            Array of source ids to plot

        plot_dir: string, optional
            Directory to save the finished plots to

        adjusted: bool, optional
            Has the light curve been divided by the average flux?

        show: bool, optional
            Should we show the plot?

        show_errors: bool, optional
            Show error bars?

        """
        for _i, path, source_id in Utilities.loop_variables(self.config, ids, adjusted=adjusted):
            _ = self.plot_light_curve(
                    source_id=source_id, curve_path=path, plot_dir=plot_dir,
                    adjusted=adjusted, show=show, close=True, show_errors=show_errors)

    ## TODO: make cleaner, allow passing optional function/array to plot_light_curve?
    def plot_given_curves_periods(self, ids, period_stats, 
            plot_dir=None, show=False, show_errors=False, adjusted=True):
        """

        Plot a light curve, as well as its primary period.

        Parameters
        ----------

        ids: numpy array
            Array of source ids to plot

        period_stats: numpy array
            Array of period, amplitude, phase, offset and errors

        plot_dir: string, optional
            Directory to save the finished plots to

        adjusted: bool, optional
            Has the light curve been divided by the average flux?

        show: bool, optional
            Should we show the plot?

        show_errors: bool, optional
            Show error bars?

        """

        for i, path, source_id in Utilities.loop_variables(self.config, ids, adjusted=adjusted):
            ## Plot base light curve
            times = self.plot_light_curve(source_id=source_id, curve_path=path, plot_dir=plot_dir,
                    adjusted=adjusted, show=False, close=False, show_errors=show_errors)

            ## If we have a nice period, plot it
            if period_stats['period'][i] > 0:
                model = period_stats['amplitude'][i] * np.sin(2*np.pi/period_stats['period'][i] * times
                        + period_stats['phi'][i]) + period_stats['offset'][i]
            else:
            ## If not, just use a flat line
                model = np.ones(len(times))
            plt.plot(times, model, color="red")


            if show:
                plt.show()

            ## Save with different prefix
            fname = "LCP_{}_{}{:04}.jpg".format(
                self.config.image_prefix, self.config.identifier, source_id)
            plt.title("Light curve for source {:04} with period {:.4f}h in {} (adjusted={})"
                .format(source_id, period_stats['period'][i]/3600, self.config.image_prefix, adjusted))

            if plot_dir == None:
                image_file = os.path.join(self.config.output_dir, fname)
            else:
                image_file = os.path.join(plot_dir, fname)

            plt.savefig(image_file)

            plt.close()


    def plot_adjusted_comparison(self, ids, plot_dir=None, show=False, show_errors=False):
        """
        Plot light curve and its adjusted version together.
        Useful for checking if the adjustment went right, or if it was a false
        source detection (like a hot pixel)

        Parameters
        ----------

        ids: numpy array
            Array of source ids to plot

        plot_dir: string, optional
            Directory to save the finished plots to

        show: bool, optional
            Should we show the plot?

        show_errors: bool, optional
            Show error bars?

        """

        for _i, _path, source_id in Utilities.loop_variables(self.config, ids, adjusted=False):
            self.plot_light_curve(
                    source_id=source_id, curve_path=None, plot_dir=plot_dir,
                    adjusted=False, show=False, close=False, show_errors=show_errors)
            self.plot_light_curve(
                    source_id=source_id, curve_path=None, plot_dir=plot_dir,
                    adjusted=True, show=show, close=True, show_errors=show_errors)


    def plot_avg_light_curve(self, curve_path, adjusted=False, show=False, show_errors=False, plot_dir=None):
        """
        Plot the light curve used as the average of bright sources.

        Parameters
        ----------

        ids: numpy array
            Array of source ids to plot

        plot_dir: string, optional
            Directory to save the finished plots to

        adjusted: bool, optional
            Has the light curve been divided by the average flux?

        show: bool, optional
            Should we show the plot?

        show_errors: bool, optional
            Show error bars?

        """

        self.plot_light_curve(source_id=None, curve_path=curve_path,
                adjusted=adjusted, show=show, show_errors=show_errors, plot_dir=plot_dir)



    def create_adjusted_light_curves(self, source_ids, stds):
        """
        Divides all light curves by the average light curve to remove noise.
        Also removes global effects of the field that vary over time.

        Creates an 'adjusted' light curve.
        Also removes cosmic rays.

        Parameters
        ----------

        source_ids: numpy array
            IDs of the sources to 'adjust'

        stds: numpy array
            Standard deviations of each light curve

        """
        #print("[DEBUG] Calling `divide_by_average` in FluxFinder")

        avg_curve = np.genfromtxt(self.config.avg_curve_path,
                dtype=self.config.light_curve_dtype).transpose()

        #for all files
        for i, path, source_id in Utilities.loop_variables(self.config, source_ids, adjusted=False):
            curve = np.genfromtxt(path, dtype=self.config.light_curve_dtype).transpose()
            n_measures = len(curve['counts'])

            ## Might not be needed since we're sigma clipping
            self.remove_cosmics(curve, stds[i])

            curve['counts'] /= avg_curve['counts']
            med = np.median(curve['counts'])
            curve['counts'] /= med
            curve['counts_err'] /= avg_curve['counts'] * med

            ## Scuffed sigma clip
            for _i in range(self.config.n_clip_iterations):
                med = np.median(curve['counts'])
                std = np.std(curve['counts'])
                clip_idx = np.where(
                        np.abs(curve['counts']-med) > std * self.config.counts_clip_threshold
                    )[0]
                #curve['counts'][clip_idx] = med     ## TODO: What do we replace it with?
                replace_idx = clip_idx + 1

                ## TODO: Useless if neighbour is also flagged as an outlier
                broken_bounds = np.where(replace_idx >= n_measures-1)[0]
                if len(broken_bounds) > 0:
                    replace_idx[broken_bounds] = clip_idx[broken_bounds] - 1

                curve['counts'][clip_idx] = curve['counts'][replace_idx]
                curve['counts_err'][clip_idx] = 1e9

            
            fname = self.config.source_format_str.format(source_id)
            out_path = os.path.join(self.config.adjusted_curve_dir, fname)

            np.savetxt(out_path, curve)

        

    ## TODO: Use set_number and image_number
    def get_thumbnail(self, image_n, x, y, size, add_shift=True):
        """
        Get a thumbnail of a source, returns pixel values.

        Parameters
        ----------

        image_n: int
            Total image number, accumulated across sets
        size: int
            Radius of box to select around center.

        """
        
        n = image_n % self.config.set_size
        
        if n == 0:
            n = self.config.set_size
            
        set_number = int((image_n-1) / self.config.set_size) + 1
        if add_shift:
            x_shift, y_shift = self.get_total_shift(set_number=set_number, image_number=n)
        else:
            x_shift = 0
            y_shift = 0
        
        self.has_sets = True
        image_file = Utilities.get_image(self.config, set_number=set_number, image_number=n)
        
        image = image_file[0].data
        
        
        x = x + x_shift
        y = y + y_shift
        
        
        ly = int(y-size)
        ry = int(y+size)
        
        lx = int(x-size)
        rx = int(x+size)
        
        if lx < 0:
            lx = 0
        
        if rx >= self.config.image_width:
            rx = self.config.image_width-1
        
        if ly < 0:
            ly = 0
        
        if ry >= self.config.image_height:
            ry = self.config.image_height - 1
            
        thumb_data = image[ly:ry, lx:rx] # note your x,y coords need to be an int
        thumb_data = np.flip(thumb_data)
            
        return thumb_data


    def remove_cosmics(self, lc, std):
        """
        Remove cosmic rays in light curve
        I don't know the algorithm but it seems to work

        Parameters
        ----------

        lc: numpy array 2d
            Light curve of source

        std: numpy array
            Standard deviation of the curve
            
        """
        
        counts = lc['counts']
        cosmic_index = -1
        
        n_measures = len(counts)
        
        for i in range(n_measures):
            
            m = counts[i]
            
            if i == 0:
                l = counts[i+1]
            else:
                l = counts[i-1]
                
            if i == len(counts) - 1:
                r = counts[i-1]
            else:
                r = counts[i+1]
        
            if m - r > self.config.cosmic_threshold * std and m - l > self.config.cosmic_threshold * std:
                
                if cosmic_index != -1:
                    return False
              
                cosmic_index = i
            
        if cosmic_index == -1:
            return False
        
        if i == 0:
            replacement = counts[1]
        else:
            replacement = counts[i-1]
                
        lc['counts'][cosmic_index] = replacement
            
        return True


    ## TODO: Magic numbers
    def create_thumbnails(self, results_table, adjusted=False, show=False):
        """
        Creates images of the brightest and dimmest frames of each source deemed variable

        Parameters
        ----------

        ff: FluxFinder
            FluxFinder object used to get a thumbnail slice

            
        """

        print("[FluxFinder] Creating thumbnails")

        for i, path, source_id in Utilities.loop_variables(self.config, results_table['id']):
            curve = np.genfromtxt(path, dtype=self.config.light_curve_dtype).transpose()
            
            c = curve['counts']

            i_dim = np.argmin(c)
            i_bright = np.argmax(c)
            
            i_x = results_table['xcentroid'][i]
            i_y = results_table['ycentroid'][i]
            
            print("[DataAnalyser] Creating thumbnail for source id {:04}, centroid {},{}"
                    .format(source_id, i_x, i_y))
            
            ## Magic numbers
            dim = self.get_thumbnail(i_dim+1, i_x, i_y, 20, True)
            bright = self.get_thumbnail(i_bright+1, i_x, i_y, 20, True)

        
            fig = plt.figure()
            fig.add_subplot(1, 2, 1)
            plt.axis('off')

            dim_norm = ImageNormalize(dim, interval=ZScaleInterval(), stretch=LinearStretch())
            plt.imshow(dim, origin='upper', cmap='gray', norm = dim_norm)

            
            fig.add_subplot(1, 2, 2)
            plt.axis('off')

            bright_norm = ImageNormalize(bright, interval=ZScaleInterval(), stretch=LinearStretch())
            plt.imshow(bright, origin='upper', cmap='gray', norm = bright_norm)
            
            fname= "thumb_{}_{}{:04}{}".format(self.config.image_prefix, self.config.identifier,
                    source_id, self.config.plot_file_extension)
            path = os.path.join(self.config.output_dir, fname)

            if show:
                plt.show()

            plt.savefig(path)
            plt.close()


    def map_id(self, id2, cat1, cat2, shifts, set_number):
        """
        Find ID of a star from catalogue 2 in catalogue 1
        Expects everything in Table objects

        Not used as of DR2

        """
        
        x_shift = 0
        y_shift = 0
        
        #find total shift between catalogues
        for i in range(set_number-1):
            x_shift += shifts['xshifts'][i]
            y_shift += shifts['yshifts'][i]
        
        #find expected position of star in catalogue 1 with id 'id2' in
        #the second catalogue
        expected_x = cat2['xcenter'][id2] - x_shift
        expected_y = cat2['ycenter'][id2] - y_shift
        
        #get x-y coords of stars in catalogue 1
        x1s = cat1['xcentroid']
        y1s = cat1['ycentroid']
        
        #find matching star id in catalogue 1 - error of 0.01 seems very large
        #perhaps implement iterative search? (increasing error if required)
        for i in range(len(x1s)):
            if Utilities.equals(expected_x, x1s[i], 0.01) and Utilities.equals(expected_y, y1s[i], 0.01):
                return i
            
        return None

