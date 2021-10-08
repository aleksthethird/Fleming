import os
from astropy.table import Table
from photutils import aperture_photometry
from photutils import CircularAperture
from photutils import CircularAnnulus
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import Constants 

import Utilities

class FluxFinder:
    
    #raw image folder path
    directory = None
    
    #directory path to contain flux data
    flux_dir = None
    
    #directory path to contain light curve
    light_curve_dir = None
    
    #image name regex
    image_names = None
    
    #list of all total shifts
    x_shifts = None
    y_shifts = None
    
    #catalogue of image 1
    catalogue = None
    
    #number of images in each set
    set_size = None
    
    #number of sets of images
    n_sets = None
    
    #is image data organised into sets?
    has_sets = None
    
    #time at which first image was taken
    obs_start_time = []
    
    times = []
    
    avg_fluxes = []
    
    
    def __init__(self, directory, image_names, has_sets, n_sets, set_size):
        self.directory = directory
        self.image_names = image_names
        self.n_sets = n_sets
        self.set_size = set_size
        self.has_sets = has_sets
        self.catalogue = Table.read(directory + Constants.working_directory + Constants.catalogue_prefix + image_names + Constants.standard_file_extension, format = Constants.table_format)

        self.get_shifts()
        
        #create flux directory if not exists
        self.flux_dir = directory + Constants.working_directory + Constants.flux_directory
        if not os.path.exists(self.flux_dir):
            os.mkdir(self.flux_dir)
        
        #create light curve directory if not exists
        self.light_curve_dir = directory + Constants.working_directory  + Constants.light_curve_directory
        if not os.path.exists(self.light_curve_dir):
            os.mkdir(self.light_curve_dir)
        
    #find the fluxes of each star in each image
    def find_all_fluxes(self):
        
        #iterate over each image
        for set in range(1, self.n_sets+1):
            for i in range(1, self.set_size+1):
                print(str(i))
                self.find_fluxes(i, set)
                
               
            
            
    #read the shifts file
    def get_shifts(self):
        
        #build shift file path
        shifts_path = self.directory + Constants.working_directory + Constants.shift_file
        
        #read shifts in as table
        t = Table.read(shifts_path, format = Constants.table_format)
        
        self.x_shifts = t['xshifts']
        self.y_shifts = t['yshifts']
        
    
    #get the total shift between the specified image and the first image
    def get_total_shift(self, image_number, set):
    
        #file is just one long list of shifts, with a shift associated with 
        #each image - below is the expression required
        #to get the shift for the image belong to set number 'set' and image 
        #number 'image_number'
        index = (set-1)*self.set_size + image_number-1
        
        total_x = self.x_shifts[index]
        total_y = self.y_shifts[index]

        return total_x, total_y
    
# =============================================================================
#     #find the fluxes of all stars ina given image
#     def find_fluxes(self, image_number, set):
# 
#         #get total shift from first image to this one
#         x_shift, y_shift = self.get_total_shift(image_number, set)
#         
#         #add the shift onto the positions of the stars in the first image
#         #to find their positions in this image
#         x = self.catalogue['xcentroid'] + x_shift
#         y = self.catalogue['ycentroid'] + y_shift
#         
#         positions = (x, y) 
#         
#         # Shape and size of aperture object
#         apertures = CircularAperture(positions, r=9) 
#         
#         #build path to image 
#         data_path = self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix + self.image_names
#         if not self.has_sets:
#             data_path += "000" + str(image_number)
#         else:
#             data_path += "_" + str(set) + "_" + Utilities.format_index(image_number)
#         
#         data_path += Constants.fits_extension
#             
#         #read in image data
#         image_data = fits.getdata(data_path, ext = 0)
#                 
#         # Local background subtraction
# 
#         # Define size of background aperture
#         annulus_apertures = CircularAnnulus(positions, r_in=10, r_out=15)
#         apertures = CircularAperture(positions, r=9)
#         apers = [apertures, annulus_apertures]
# 
#         # find counts in each aperture and annulus
#         phot_table2 = aperture_photometry(image_data, apers)
# 
#         for col in phot_table2.colnames:
#             phot_table2[col].info.format = '%.8g'  # for consistent table output
#     
#         # MEAN
#         # Background sub using mean
#             
#         phot_table2['residual_aperture_sum_mean'] = phot_table2['aperture_sum_1']*0 - 1
#         phot_table2['mean'] = phot_table2['aperture_sum_1']*0 - 1
# 
# 
#         for i in range(len(phot_table2)):
#             
#             x = positions[0][i]
#             y = positions[1][i]
#             
#             
#             #check that largest aperture does not exceed the boundaries of the image
#             if Utilities.is_within_boundaries(x, y, len(image_data[0]), len(image_data), 15):
#             
#                 # Calc the mean background in the second aperture ring
#                 bkg_mean = phot_table2['aperture_sum_1'][i] / annulus_apertures.area()
#                 phot_table2['mean'][i] = bkg_mean
#                 
#                 # Calc background level in each main aperture and subtract
#                 bkg_sum = bkg_mean * apertures.area()
#                 final_sum = phot_table2['aperture_sum_0'][i] - bkg_sum
#                 phot_table2['residual_aperture_sum_mean'][i] = final_sum
#         
#         phot_table2['residual_aperture_sum_mean'].info.format = '%.8g'  # for consistent table output
#                 
#         # MEDIAN
# 
#         # Background sub using median
#         phot_table2['median'] = phot_table2['aperture_sum_1']*0 - 1
#         phot_table2['residual_aperture_sum_med'] = phot_table2['aperture_sum_1']*0 - 1
# 
# 
#         #for each source
#         for q in range(0,len(phot_table2)):
#             x = positions[0][q]
#             y = positions[1][q]   
#             xypos = (x, y)
#             
#             #check that largest aperture does not exceed the boundaries of the image
#             if Utilities.is_within_boundaries(x, y, len(image_data[0]), len(image_data), 15):
#                 annulus = CircularAnnulus(xypos, r_in=10, r_out=15)
#                 ann_mask = annulus.to_mask(method='center')[0]
#                 weighted_data = ann_mask.multiply(image_data)
#                 phot_table2['median'][q] = np.median(weighted_data[weighted_data != 0])
#                 
#                 # Calc the median background in the second aperture ring
#                 bkg_med = phot_table2['median'][q]
#             
#                 # Calc background level in each main aperture and subtract
#                 bkg_sum = bkg_med * apertures.area()
#                 final_sum = phot_table2['aperture_sum_0'][q] - bkg_sum
#             
#             
#                 phot_table2['residual_aperture_sum_med'][q] = final_sum
#         
#         phot_table2['residual_aperture_sum_med'].info.format = '%.8g'  # for consistent table output
#         
#         #build path for flux table output
#         out_path = self.flux_dir + Constants.flux_prefix + self.image_names
#         
#         if not self.has_sets:
#             out_path += "000" + str(image_number)
#         else:
#             out_path += "_" + str(set) + "_" + Utilities.format_index(image_number)
#             
#         out_path += Constants.standard_file_extension
#         phot_table2.write(out_path,  format = Constants.table_format, overwrite = True)
#         
#     def make_light_curves(self):
#         
#         #get time of observation for each image
#         times_file = self.directory + Constants.working_directory + Constants.time_file
#         
#         times = [line.rstrip('\n') for line in open(times_file)]
#         
#         light_curves = []
#           
#         #iterate over each flux table
#         for set in range(1, self.n_sets+1):
#             for i in range(1, self.set_size + 1):
#                 print(i)
#                 
#                 #build flux path
#                 file = self.flux_dir + Constants.flux_prefix + self.image_names
#                 if not self.has_sets:
#                     file += "000" + str(i) 
#                 else:
#                     file += "_" + str(set) + "_" + Utilities.format_index(i)
#                 file += Constants.standard_file_extension
#                 
#                 #read flux data
#                 t = Table.read(file, format = Constants.table_format)
#                 
#                 #iterate through each source in the table
#                 for j in range(len(t['id'])):
#                     
#                     #add table objects to the light curves array when it is
#                     #empty
#                     if i == 1 and set == 1:
#                         light_curves.append(Table(names = ('time','counts')))
#                         
#                     #get time matching the image
#                     date_and_time = times[(set-1) * self.set_size + (i-1)]
#                     
#                     #just interested in time part of date and time for now
#                     time = date_and_time.split("T")[1]
#                     
#                     #split time into hours, minutes and seconds
#                     time_split = time.split(":")
#                     time_split[0] = int(time_split[0])
#                     time_split[1] = int(time_split[1])
#                     time_split[2] = int(time_split[2])
#                     
#                     #if the observation start time is empty, add the current
#                     #time in hours minutes and seconds 
#                     if len(self.obs_start_time) == 0:
#                         self.obs_start_time = time_split
#                     
#                     #get difference between time of start and time of 
#                     #observation of current image, in hours minutes and seconds
#                     time_elapsed_hms = Utilities.diff(time_split, self.obs_start_time)
#                     
#                     #if the hours value is positive this observation was taken before midnight
#                     #the time elapsed from the start would be negative as the hour value
#                     #woudl reset to 0
#                     if time_elapsed_hms[0] >= 0:
#                         #get time elapsed in seconds
#                         time_elapsed = time_elapsed_hms[0] * 3600 + time_elapsed_hms[1] * 60 + time_elapsed_hms[2]
#                     else:
#                         time_elapsed_hms_till_midnight = Utilities.diff([23, 59, 60], self.obs_start_time)
#                         time_elapsed = time_elapsed_hms_till_midnight[0] * 3600 + time_elapsed_hms_till_midnight[1] * 60 + time_elapsed_hms_till_midnight[2]
#                         time_elapsed += time_split[0] * 3600 + time_split[1] * 60 + time_split[0]
#                         
#                     
#                     #if median has reasonable value
#                     if t['median'][j] > 0: 
#                         light_curves[j].add_row((time_elapsed, str(t['residual_aperture_sum_med'][j])))
#                         
# 
#         #loop through all light curves, writing them out to a file
#         for j in range(len(light_curves)):
#             #build light curve path
#             file = self.light_curve_dir + self.image_names + Constants.identifier + str(j+1) + Constants.standard_file_extension
#            
#             table = light_curves[j]
#             table.write(file, format = Constants.table_format, overwrite=True)
#             
# =============================================================================
    

    #find the fluxes of all stars ina given image
    def find_fluxes(self, image_number, set):

        #get total shift from first image to this one
        x_shift, y_shift = self.get_total_shift(image_number, set)
        
        #add the shift onto the positions of the stars in the first image
        #to find their positions in this image
        x = self.catalogue['xcentroid'] + x_shift
        y = self.catalogue['ycentroid'] + y_shift
        
        positions = []
        
        for i in range(len(x)):
            positions.append((x[i], y[i]))
            
        
        # Shape and size of aperture object
        apertures = CircularAperture(positions, r=9) 
        
        image_data = self.get_image(image_number, set)[0].data
                
        # Local background subtraction

        # Define size of background aperture
        annulus_apertures = CircularAnnulus(positions, r_in=10, r_out=15)
        apertures = CircularAperture(positions, r=9)
        apers = [apertures, annulus_apertures]

        # find counts in each aperture and annulus
        phot_table2 = aperture_photometry(image_data, apers)

        for col in phot_table2.colnames:
            phot_table2[col].info.format = '%.8g'  # for consistent table output
    
        # MEAN
        # Background sub using mean
            
        phot_table2['residual_aperture_sum_mean'] = phot_table2['aperture_sum_1']*0 - 1
        phot_table2['mean'] = phot_table2['aperture_sum_1']*0 - 1


        for i in range(len(phot_table2)):
            
            x, y = positions[i]

            
            
            #check that largest aperture does not exceed the boundaries of the image
            if Utilities.is_within_boundaries(x, y, len(image_data[0]), len(image_data), 15):


                # Calc the mean background in the second aperture ring
                bkg_mean = phot_table2['aperture_sum_1'][i] / annulus_apertures.area
                phot_table2['mean'][i] = bkg_mean
                
                # Calc background level in each main aperture and subtract
                bkg_sum = bkg_mean * apertures.area
                final_sum = phot_table2['aperture_sum_0'][i] - bkg_sum
                phot_table2['residual_aperture_sum_mean'][i] = final_sum
        
        phot_table2['residual_aperture_sum_mean'].info.format = '%.8g'  # for consistent table output
                
        # MEDIAN

        # Background sub using median
        phot_table2['median'] = phot_table2['aperture_sum_1']*0 - 1
        phot_table2['residual_aperture_sum_med'] = phot_table2['aperture_sum_1']*0 - 1


        #for each source
        for q in range(0,len(phot_table2)):
            x, y = positions[q]
 
            xypos = (x, y)

            #check that largest aperture does not exceed the boundaries of the image
            if Utilities.is_within_boundaries(x, y, len(image_data[0]), len(image_data), 15):
                annulus = CircularAnnulus(xypos, r_in=Constants.inner_radius, r_out=Constants.outer_radius)
                ann_mask = annulus.to_mask(method='center')
                weighted_data = ann_mask.multiply(image_data)
                phot_table2['median'][q] = np.median(weighted_data[weighted_data != 0])
                
                # Calc the median background in the second aperture ring
                bkg_med = phot_table2['median'][q]
            
                # Calc background level in each main aperture and subtract
                bkg_sum = bkg_med * apertures.area
                final_sum = phot_table2['aperture_sum_0'][q] - bkg_sum
            
            
                phot_table2['residual_aperture_sum_med'][q] = final_sum
        
        phot_table2['residual_aperture_sum_med'].info.format = '%.8g'  # for consistent table output
        
        return phot_table2
        
    
    def get_image(self, image_number, set):
        
        #build path to image 
        data_path = self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix + self.image_names
        if not self.has_sets:
            data_path += "000" + str(image_number)
        else:
            data_path += "_" + str(set) + "_" + Utilities.format_index(image_number)
        
        data_path += Constants.fits_extension
            
        #read in image data
        return fits.open(data_path, ext = 0)
        
    def make_light_curves(self):
        
        #get time of observation for each image
        times_file = self.directory + Constants.working_directory + Constants.time_file
        
        times = [line.rstrip('\n') for line in open(times_file)]
        
        light_curves = []
          
        #iterate over each image
        for set in range(1, self.n_sets+1):
            for i in range(1, self.set_size+1):
                
                print("Finding fluxes in image " + str((set-1)*self.set_size + i))
                
                t = self.find_fluxes(i, set)
                                
                #iterate through each source in the table
                for j in range(len(t['id'])):
                    
                    #add table objects to the light curves array when it is
                    #empty
                    if i == 1 and set == 1:
                        light_curves.append(Table(names = ('time','counts')))
                        
                    #get time matching the image
                    date_and_time = times[(set-1) * self.set_size * 2 + (2*i-2)]
                    
                    #just interested in time part of date and time for now
                    time = date_and_time.split("T")[1]
                    
                    #split time into hours, minutes and seconds
                    time_split = time.split(":")
                    time_split[0] = int(time_split[0])
                    time_split[1] = int(time_split[1])
                    time_split[2] = int(time_split[2])
                    
                    #if the observation start time is empty, add the current
                    #time in hours minutes and seconds 
                    if len(self.obs_start_time) == 0:
                        self.obs_start_time = time_split
                    
                    #get difference between time of start and time of 
                    #observation of current image, in hours minutes and seconds
                    time_elapsed_hms = Utilities.diff(time_split, self.obs_start_time)
                    
                    #if the hours value is positive this observation was taken before midnight
                    #the time elapsed from the start would be negative as the hour value
                    #woudl reset to 0
                    if time_elapsed_hms[0] >= 0:
                        #get time elapsed in seconds
                        time_elapsed = time_elapsed_hms[0] * 3600 + time_elapsed_hms[1] * 60 + time_elapsed_hms[2]
                    else:
                        time_elapsed_hms_till_midnight = Utilities.diff([23, 59, 60], self.obs_start_time)
                        time_elapsed = time_elapsed_hms_till_midnight[0] * 3600 + time_elapsed_hms_till_midnight[1] * 60 + time_elapsed_hms_till_midnight[2]
                        time_elapsed += time_split[0] * 3600 + time_split[1] * 60 + time_split[0]
                        
                    
                    #if median has reasonable value
                    if t['median'][j] > 0: 
                        light_curves[j].add_row((time_elapsed, str(t['residual_aperture_sum_med'][j])))

        #loop through all light curves, writing them out to a file
        for j in range(len(light_curves)):
            #build light curve path
            file = self.light_curve_dir + self.image_names + Constants.identifier + str(self.catalogue['id'][j]) + Constants.standard_file_extension
            table = light_curves[j]
            table.write(file, format = Constants.table_format, overwrite=True)
            
    
    
    
    
    #find id of a star from catalogue 2 in catalogue 1
    def map_id(self, id2, cat1, cat2, shifts, set):
        
        x_shift = 0
        y_shift = 0
        
        #find total shift between catalogues
        for i in range(set-1):
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
        
        
            
       
            
        
        
        
    #plot light curve of star with the given id (from catalogue)
    def plot_light_curve(self, id, file, adjusted):
        
        if file == None:
            if not adjusted:
                file = self.light_curve_dir + self.image_names + Constants.identifier + str(id) + Constants.standard_file_extension
            else:
                file = self.directory + Constants.working_directory + Constants.adjusted_curves_directory + self.image_names + Constants.identifier + str(id) + Constants.standard_file_extension

        table = Table.read(file, format=Constants.table_format)
        
        times = table['time']
        fluxes = table['counts']
        
        mean = Utilities.mean(fluxes)
        
        for i in range(len(fluxes)):
            fluxes[i] = fluxes[i]/mean
        
        minimum = min(fluxes)
        maximum = max(fluxes)
    
        _ = plt.figure(figsize = (12, 8))
        _ = plt.scatter(times, fluxes, s=5, marker = 'x');
        _ = plt.xlabel("time (s)")
        _ = plt.ylabel("counts/mean")
        _ = plt.ylim(0.9 * minimum, maximum * 1.1)
        
        if id == None:
            id = 'avg'
     
            
        _ = plt.savefig(Constants.folder + Constants.working_directory + Constants.output_directory + self.image_names + "_id_" + str(id) + "_LC")

        
    
    #makes an average light curve with the brightest stars to subtract noise
    #from the resulting curves
    
    #use median
    def make_avg_curve(self, ids):
                
        #iterate through each star in the given list
        for i in range(len(ids)):
            file = self.light_curve_dir + self.image_names + Constants.identifier + str(ids[i]) + Constants.standard_file_extension
            
            #get the light curve for star i
            t = Table.read(file, format = Constants.table_format)
            
            fluxes = t['counts']
            
            if i == 0:
                for j in range(len(fluxes)):
                    self.avg_fluxes.append(0)
                self.times = t['time']

            #find the mean flux of star i
            mean = Utilities.mean(fluxes)
            
            #add each flux measurement in the light curve at time index k 
            #of star i to a list which will store the average flux at each point
            #in time k. Note the flux added to the list is divided by the mean
            #to reduce any bias(?). The average contribution to each k is 1.
            for k in range(len(fluxes)):
                self.avg_fluxes[k] = self.avg_fluxes[k] + (fluxes[k] / mean)
            
        #calculate average flux at each time of measurement l
        for l in range(len(self.avg_fluxes)):
            self.avg_fluxes[l] = self.avg_fluxes[l] / len(ids)
                
        #generate table of average light curve
        light_curve = Table([self.times, self.avg_fluxes], names = ('time','counts') )

        #save average light curve
        file = self.directory + Constants.working_directory + self.image_names + "_avg" + Constants.standard_file_extension
        light_curve.write(file, format = Constants.table_format, overwrite=True)

    #divide all light curves by the average light curve. Removes noise
    def divide_by_average(self):
        
        adjusted_light_curve_dir = self.directory + Constants.working_directory  + Constants.adjusted_curves_directory
        if not os.path.exists(adjusted_light_curve_dir):
            os.mkdir(adjusted_light_curve_dir)        
                    
        #for all files
        for file in os.listdir(self.light_curve_dir):
            if file[:len(self.image_names)] == self.image_names:
                
                t = Table.read(self.light_curve_dir + file, format = Constants.table_format)
                                                
                this_fluxes = t['counts']
                this_times = t['time']
                
                id = file.split("id")[1].split(".")[0]
                
                for i in range(len(this_fluxes)):
                    time = this_times[i]
                    
                    #divide each flux measurement by the average flux measurement
                    #at its respective point in time
                    for j in range(len(self.avg_fluxes)):
                        if time == self.times[j]:
                            this_fluxes[i] = this_fluxes[i] / self.avg_fluxes[j]
                            #print(this_fluxes[i])
                
                #export adjusted light curve
                light_curve = Table([this_times, this_fluxes], names = ('time','counts'))
                out_file = adjusted_light_curve_dir + self.image_names + Constants.identifier + str(id) + Constants.standard_file_extension
                light_curve.write(out_file, format = Constants.table_format, overwrite=True)

        
    def get_thumbnail(self, image_n, x, y, size, add_shift):
        
        n = image_n % self.set_size
        
        if n == 0:
            n = self.set_size
            
        set = int((image_n-1) / self.set_size) + 1
        if add_shift:
            x_shift, y_shift = self.get_total_shift(n, set)
        else:
            x_shift = 0
            y_shift = 0
        
        self.has_sets = True
        image_file = self.get_image(n, set)
        
        image = image_file[0].data
        
        
        x = x + x_shift
        y = y + y_shift
        
        
        ly = int(y-size)
        ry = int(y+size)
        
        lx = int(x-size)
        rx = int(x+size)
        
        if lx < 0:
            lx = 0
        
        if rx > Constants.image_width:
            rx = Constants.image_width-1
        
        if ly < 0:
            ly = 0
        
        if ry > Constants.image_height:
            ry = Constants.image_height - 1
            
            
        return image[ly:ry, lx:rx] # note your x,y coords need to be an int

        
        
        
        
        
        
        
        
            
            
        
            
    
    
    
    
            