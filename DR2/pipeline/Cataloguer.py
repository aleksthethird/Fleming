from photutils import DAOStarFinder
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.table import Table
from astroquery.astrometry_net import AstrometryNet

import os
import matplotlib.pyplot as plt

from . import Utilities, Config


class Cataloguer:
    """
    Cataloguer object.

    Anything to do with cataloguing all stars in a field, finding coordinates etc.
    Only need one `Cataloguer` object per field.

    Uses DAO star finder to get ids and basic info.
    Uses astrometry.net local solving to find coordinate system.

    Attributes
    ----------

    config: Config
        Config object for the field.

    n_sources: int
        Number of sources that the cataloguer found.

    wcs: WCS
        World Coordinate System (thing that holds ra/dec information from astrometry)
    """

    config = None           # Config object
    n_sources = 0           # Number of sources discovered
    wcs = None              ## TODO: Do we need to keep this?

    
    def __init__(self, config):
        """
        Cataloguer constructor.

        Parameters
        ----------

        config: Constants.Config
            Config object used for the run
        """

        self.config = config


    def generate_catalogue(self, image_path, solve=True, write_file=True):
        """
        Generates a catalogue of all the stars in the given image.
        Catalogue location is defined in self.config.catalogue_prefix.

        Parameters
        ----------

        image_path: string
            Absolute path to the fits image to catalogue all the stars.

        solve: bool, optional
            Should we solve astrometry.net wcs?

        Returns
        -------

        n_sources: int
            How many sources we found in the image.

        """

        ## Read in image data
        image_data = fits.getdata(image_path, ext=0)
        
        ## TODO: Magic number, threshold passed to find_stars
        ## Build a catalogue of all stars in the image
        ## sources is a Table object
        sources = self.find_stars(image_data, 5)
        self.filter_stars(sources, image_data)
        self.n_sources = len(sources['id'])

        print("[Cataloguer] Found {} suitable sources".format(self.n_sources))

        ## Write all sources to a reg file
        ## TODO: Necessary? we already have a catalogue file.
        #Utilities.make_reg_file(self.config.workspace_dir, self.config.image_prefix, sources)
        
        ## Add the RA and DEC for each star to the catalogue
        self.convert_to_ra_and_dec(image_path, sources, solve=solve)
        
        if write_file:
            ## Write the catalogue to the catalogue file
            sources.write(self.config.catalogue_path, format=self.config.table_format, overwrite=True)

        return self.n_sources

    


    def filter_stars(self, sources, image_data):
        """
        Removes stars from a source catalogue given cutoff parameters
        defined in the config object.

        Note: 'id' column in sources is now no longer a continuous
        range from 0..n_sources

        Parameters
        ----------

        sources: Table
            Information on all of the sources found.
            
        image_data: fits data
            Data of the fits image used to create the catalogue.
            
        """
        
        to_remove = []
        n_sources = len(sources)
                
        ## Loop over each source
        for i in range(n_sources):
            
            ## Conditions to remove a source
            ## Note: is_within_boundaries does not account for later shifts
            is_too_bright = sources['flux'][i] > self.config.max_flux_cutoff
            is_too_dim = sources['flux'][i] < self.config.min_flux_cutoff
            is_within_boundaries = Utilities.is_within_boundaries(
                    sources['xcentroid'][i],
                    sources['ycentroid'][i],
                    self.config.image_width,
                    self.config.image_height,
                    self.config.edge_limit)

            if is_too_bright or is_too_dim or not is_within_boundaries:
                to_remove.append(i)
        

        ## Remove rows from sources table
        for i in range(len(to_remove)):
            sources.remove_row(to_remove[i] - i)
    
            
        print("[Cataloguer] Filtered out {} objects".format(len(to_remove)))
        



    def convert_to_ra_and_dec(self, image_file, sources, solve=True):
        """
        Converts all the source positions in the image to RA and DEC.
        Uses World Coordinate System (wcs) from the fits image.
        Does this via astrometry.net api.

        Sets the self.wcs variable to the WCS header found.

        Parameters
        ----------

        image_file: string
            Absolute path to the image to use as catalogue.

        sources: Table
            Table of information on all the sources.

        solve: bool, optional
            Should we solve the header using astrometry.net?

        """
        
        print("[Cataloguer] Getting coordinate system for image '{}'".format(image_file))

        ## Find the wcs assosiated with the fits image using astropy and the header
        self.wcs, _head = self.get_wcs_header(image_file, solve=solve)
        
        ## Make two new coloums in source table, set to zero
        sources['RA'] = sources['xcentroid'] * 0
        sources['DEC'] = sources['xcentroid'] * 0
    
        ## Replace the zeroes with ra and dec
        for i in range(self.n_sources):
            ra, dec = self.wcs.all_pix2world(sources['xcentroid'][i], sources['ycentroid'][i], 0) 
            sources['RA'][i] = ra
            sources['DEC'][i] = dec
    

    
    def get_wcs_header(self, file, solve=True):
        """
        Get World Coordinate System header.
        [Astroquery docs.](https://astroquery.readthedocs.io/en/latest/astrometry_net/astrometry_net.html)

        Parameters
        ----------

        file: string
            Absolute path to file to solve with astrometry.net

        solve: bool, optional
            Should we do the solve? (used for debug as the process can be slow)

        Returns
        -------

        wcs: WCS
            FITS world coordinate system header

        """

        ## Set up astrometry.net instance
        ast = AstrometryNet()
        ast.TIMEOUT = self.config.astrometry_timeout
        ast.api_key = self.config.api_key
    
        wcs = None

        ## If we don't want to solve, return empty header
        if not solve:
            return WCS(header=wcs), wcs
    

        print("[Astrometry] Starting job")
        try:
            ## Try to solve on local machine
            wcs = ast.solve_from_image(
                    file, 
                    solve_timeout=self.config.astrometry_timeout,
                    force_image_upload=False)
        except Exception as e:
            ## If we failed, print why (likely timeout)
            print("\n[Astrometry] Error: WCS solve failed")
            print("............ Exception: {}".format(e))
            print("............ (self.config.astrometry_timeout={}s)"
                    .format(self.config.astrometry_timeout))
        
        print("\n[Astrometry] Finished job")
            
        ## TODO: Gives warning with no image data
        return WCS(header=wcs), wcs




    def generate_image_times(self):
        """
        Get a list of times of each image/measurement.
        File path specified in self.config.

        """

        print("[Cataloguer] Generating image time file")
        ## Clear time file if it exists
        if(os.path.exists(self.config.time_path)):
            open(self.config.time_path, "w").close()

        ## Open time file
        with open(self.config.time_path, "a+") as time_file:
    
            ## Loop through all images
            for image_name, _s, _i in Utilities.loop_images(self.config):
            
                ## Get the image header
                image_file = os.path.join(self.config.image_dir, image_name)
                image_header = fits.getheader(image_file)

                ## Write time to file
                ## TODO: Sometimes get warning on name
                time_file.write("{}{}".format(
                            image_header['DATE-OBS'],
                            self.config.line_ending))

        

    ## TODO: Magic numbers, make configs
    def find_stars(self, image_data, threshold):
        """
        Catalogue all sources that meet the thresholds in the image.
        Uses DAOStarFinder.
    
        Parameters
        ----------
    
        image_data: fits image data
            Image data to find stars in.
    
        threshold: float
            Number of standard deviations above the background which
            an object is classed as a star.
    
        Returns
        -------
    
        sources: Table
            Table of all the information DAO found about the star.
    
        """
        
        ## Get mean median and standard deviation of the image data
        mean, median, std = sigma_clipped_stats(image_data, sigma=3.0, maxiters=5)  
        
        ## Initiate finder object. Will find objects with a FWHM of 8 pixels
        ## and 3-sigma times the background
        daofind = DAOStarFinder(fwhm=8, threshold=threshold*std) 
        
        ## Find sources
        sources = daofind(image_data)
    
        ## Set table format
        for col in sources.colnames:    
            sources[col].info.format = '%.8g'
           
        return sources

                    

