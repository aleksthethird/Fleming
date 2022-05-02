from astropy.table import Table
from astropy.io import fits
from photutils.centroids import centroid_2dg, centroid_sources

import os
import numpy as np

from . import Utilities, Config

class ShiftFinder:

    config = None
    n_sources = None
    reference_ids = None
    
    def __init__(self, config, n_sources):
        self.config = config
        self.n_sources = n_sources

    
    ## TODO: Delete comments? Original returned a slice
    def get_reference_coordinates(self, n_reference_stars):
        """
        Finds pixel coordinates of n brightest stars in catalogue

        """

        ## Debug, skip these N brightest stars
        skip = 0
        
        #reads in catalogue
        catalogue = Utilities.read_catalogue(self.config)

        xs = catalogue["xcentroid"]
        ys = catalogue["ycentroid"]
        fluxes = catalogue["flux"]
        ids = catalogue["id"]
    
        sorted_indices = np.flip(np.argsort(fluxes))

        brightest_xs = xs[sorted_indices][skip:skip+n_reference_stars]
        brightest_ys = ys[sorted_indices][skip:skip+n_reference_stars]

        self.reference_ids = ids[sorted_indices]

        return brightest_xs, brightest_ys
    
    def get_reference_ids(self):
        return self.reference_ids


    ## TODO: Fix fit warnings
    def find_shift(self, previous_x, previous_y, image_path, size=10):
        """
        
        Parameters
        ----------

        previous_x, previous_y: int
            Previous x and y pixel coordinates of star
        image_path: string
            Path to the next image to find shifts.
        size: int, optional, default=10
            Radius of box around star. Value=10 makes 20x20 square box.

        Returns
        -------
        x_shift: float
            Shift since last image in x direction.
        y_shift: float
            Shift since last image in y direction.

        """
        
        # open the next image and load it as an array
        data_array = fits.open(image_path)[0].data

        xlow  = int(previous_x)-size
        xhigh = int(previous_x)+size

        if xlow < 0: 
            xlow = 0
        if xhigh > self.config.image_width:
            xhigh = self.config.image_width

        ylow  = int(previous_y)-size
        yhigh = int(previous_y)+size

        if ylow < 0: 
            ylow = 0
        if yhigh > self.config.image_height:
            yhigh = self.config.image_height

        ## 7 comes from centroid_2dg requirements
        ## TODO: Pick new reference star if this triggers
        ## TODO: Make better
        if (yhigh-ylow < 7) or (xhigh-xlow < 7):
            print("[ShiftFinder] Error: Reference star at {:.02},{:.02} out of frame"
                    .format(previous_x, previous_y))
            return previous_x, previous_y

        # Select just the box around the star
        data_square = data_array[ylow:yhigh,xlow:xhigh]

        # Get coordinates of the centre of the star in the data square
        # Fits a 2d Gaussian to data
        x, y = centroid_2dg(data_square)

        ## TODO: Why this weird way of calculating?
        #calculate shift between previous position and newly calculated
        #positions        
        x_shift = ((x - size) - (previous_x - int(previous_x)))
        # also correcting for the shift due to using int
        y_shift = ((y - size) - (previous_y - int(previous_y)))
        
        return x_shift, y_shift



    def generate_shifts(self):
        
        #empty shift file if it exists
        if(os.path.exists(self.config.shift_path)):
            open(self.config.shift_path, "w").close()
        
        n_images = self.config.n_sets*self.config.set_size

        star_shifts = np.empty(shape=(2, n_images))
        
        #get coordinates of N bright stars in image to use as a reference
        ref_x, ref_y = self.get_reference_coordinates(self.config.n_reference_stars)
        

        ## Initial shifts are zero
        prev_x_shift = 0
        prev_y_shift = 0


        ## Shifts of the reference stars
        ## Rewritten after each image
        ref_shifts_x = np.zeros(self.config.n_reference_stars)
        ref_shifts_y = np.zeros(self.config.n_reference_stars)

        #iterate through each image in each set
        j = 0
        for fname, s, i in Utilities.loop_images(self.config):
            image_path = os.path.join(self.config.image_dir, fname)
            print("[ShiftFinder] Finding shifts in image: set {:1}; image {:03}...".format(s, i))
                
            
            for k in range(self.config.n_reference_stars):

                #find shift between the x and y of reference star k in the 
                #previous image and the current image
                ref_shifts_x[k], ref_shifts_y[k] = self.find_shift(ref_x[k], ref_y[k], image_path)
                
                
            ## Get median shift
            med_shift_x = np.median(ref_shifts_x)
            med_shift_y = np.median(ref_shifts_y)

            ## Numpy does the loop for us
            ref_x += med_shift_x
            ref_y += med_shift_y
            
            star_shifts[0][j] = med_shift_x + prev_x_shift
            star_shifts[1][j] = med_shift_y + prev_y_shift

            ## Set the previous shift for next loop
            ## Next image's previous shift is our current shift
            prev_x_shift = star_shifts[0][j] 
            prev_y_shift = star_shifts[1][j] 

            print("[ShiftFinder] ...Shifts: {:.4f},{:.4f}".format(star_shifts[0][j], star_shifts[1][j]))
            j += 1

        ## Usually catalogue is first image, but we want
        ## all shifts relative to any general catalogue image
        catalogue_index = (self.config.catalogue_set_number-1)*self.config.set_size \
                + (self.config.catalogue_image_number-1)
        star_shifts[0] -= star_shifts[0][catalogue_index]
        star_shifts[1] -= star_shifts[1][catalogue_index]
        
        
        ## Write shifts to file
        np.savetxt(self.config.shift_path, np.transpose(star_shifts))


    

## TODO: Remove
#    #I believe this is redundant
#    def find_shift_between_all_catalogues(self, image_size):  
#        
#        previous_cat = Table.read(self.directory + Constants.working_directory + Constants.catalogue_prefix + self.image_names + "_1" + Constants.standard_file_extension, format = Constants.table_format)
#        x_shifts = []
#        y_shifts = []
#        
#        for i in range(2, self.n_sets+1):
#            cat = Table.read(self.directory + Constants.working_directory + Constants.catalogue_prefix + self.image_names + "_" + str(i) + Constants.standard_file_extension, format = Constants.table_format)
#            x_shift, y_shift = self.find_shift_between_catalogues(previous_cat, cat, image_size)
#            
#            x_shifts.append(x_shift)
#            y_shifts.append(y_shift)
#            
#            previous_cat = cat
#            
#        shift_file = self.directory + Constants.working_directory + Constants.catalogue_shift_file
#
#        table = Table([x_shifts, y_shifts], names = ('xshifts','yshifts'))
#        
#        table.write(shift_file, format = Constants.table_format, overwrite=True)
            
# =============================================================================
# #I believe this is redundant 
# def find_shift_between_catalogues(catalogue1, catalogue2):
# 
#     x1s = catalogue1['xcentroid']
#     y1s = catalogue1['ycentroid']
#     
#     x2s = catalogue2['xcentroid']
#     y2s = catalogue2['ycentroid']
#     
#     fluxes = catalogue1["flux"]
#     
#     max = 0
#     
#     for i in range(len(fluxes)):
#         if float(fluxes[i]) > fluxes[max] and x1s[i] > Constants.image_width - 200 and x1s[i] < Constants.image_width + 200 and y1s[i] > Constants.image_height - 200 and y1s[i] < Constants.image_height + 200:
#             max = i
#     
#     x = x1s[max]
#     y = y1s[max]
#     
#     distances = set()
#     
#     for i in range(len(x1s)):
#                 
#         if i != max:
#             shifts = [0, 0]
#         
#             shifts[0] = round(x1s[i] - x)
#             shifts[1] = round(y1s[i] - y)
#             
#             distances.add(str(shifts[0]) + " " +  str(shifts[1]))
#             
#     for i in range(len(x2s)):
#         matches = 0
#         for j in range(len(x2s)):
#             if i != j:
#                 shift = [round(x2s[j] - x2s[i]), round(y2s[j] - y2s[i])]
#                 
#                 for k in range(len(distances)):
#                     if distances[k][0] - 3 < shift[0] and distances[k][0] + 3 > shift[0] and distances[k][1] - 3 < shift[1] and distances[k][1] + 3 > shift[1] and not k in matched:
#                         matches += 1
#                 
#                 if j > 0.01 * len(x2s) and matches < 0.1 * j:
#                     #print(matches, j)
# 
#                     break
#         
#         print(matches)
#         if matches > 0.3*len(distances):
# 
#             return x2s[i] - x1s[max], y2s[i] - y1s[max]
#                 
# 
#     
# =============================================================================
    


        
        
    

