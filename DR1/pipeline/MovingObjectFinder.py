from photutils import DAOStarFinder
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
import os
import Constants
import Cataloguer
import Utilities
from astropy.table import Table
import matplotlib.pyplot as plt
import FluxFinder
from astropy.table import Table
from astroquery.astrometry_net import AstrometryNet
import ShiftFinder
from photutils import centroid_2dg
import sys
import numpy as np


class MovingObjectFinder():
    
    
    directory = None

    catalogues = []
    removed_catalogues = []
    exp_positions = []

    image_width = None
    image_height = None
    catalogue = None
    x_shift = None
    y_shift = None
    first_image = None
    second_image = None
    means = []
    
    def __init__(self, dir):
        self.directory = dir
    
    def catalogue(self, file):
        
        self.first_image = fits.getdata(file, ext=0)
        
        sources = Cataloguer.find_stars(self.first_image, 3)
        

        
        return sources
    
    def get_images(self):
        self.first_image = self.get_image(1)
        self.second_image = self.get_image(Constants.moving_obj_check_image)

    
    #read the shifts file
    def get_shifts(self):
        
        #build shift file path
        shifts_path = self.directory + Constants.working_directory + Constants.shift_file
        
        #read shifts in as table
        t = Table.read(shifts_path, format = Constants.table_format)
        
        
        self.x_shift = t['xshifts'][0] - t['xshifts'][Constants.moving_obj_check_image]
        self.y_shift = t['yshifts'][0] - t['yshifts'][Constants.moving_obj_check_image]

                
        
            
        


            
    

    def find_moving_objects(self):
        
        self.get_images()        
        self.get_shifts()
        self.catalogue = Cataloguer.find_stars(self.first_image, 3)

        indexes = []
        plotx = []
        ploty = []
        
        catalogue_size = 30
        centroid_size = 9
        
        ff = FluxFinder.FluxFinder(Constants.folder, Constants.file_name, True, 7, 50)

        for i in range(len(self.catalogue['id'])):
           
            exp_x = self.catalogue['xcentroid'][i] - self.x_shift
            exp_y = self.catalogue['ycentroid'][i] - self.y_shift
            #not asure if image dimensions are the right way round
            if not Utilities.is_within_boundaries(exp_x, exp_y, len(self.second_image[0]), len(self.second_image), 200):
                #print("Not within the boundaries ", exp_x, exp_y)
                continue

            # Select just the box around the star
            catalogue_square = self.second_image[int(exp_y)-catalogue_size:int(exp_y)+catalogue_size, int(exp_x)-catalogue_size:int(exp_x)+catalogue_size] # note your x,y coords need to be an int

            #get coordinates of the centre of the star in the data square


            sources = Cataloguer.find_stars(catalogue_square, 6)

            obj_x = None
            obj_y = None

            for j in range(len(sources['id'])):
                
                x = sources['xcentroid'][j]
                y = sources['ycentroid'][j]
                
                if x > catalogue_size - 9 and x < catalogue_size + 9 and y > catalogue_size - 9 and y < catalogue_size + 9:
                    obj_x = x - catalogue_size
                    obj_y = y - catalogue_size
                    break
            
            

            


            #x = x - 10
            #y = y - 10
            
            #calculate shift between previous position and newly calculated
            #positions        
            #x_shift = ((x - size) - (exp_x - int(exp_x)))
            # also correcting for the shift due to using int
            #y_shift = ((y - size) - (exp_y - int(exp_y)))
            
            #if x_shift < self.x_shift - 10 or x_shift > self.x_shift + 10 or y_shift > self.y_shift - 10 or y_shift < self.y_shift + 10:
           
            if not obj_x:    
             
                object_square = self.second_image[int(exp_y)-centroid_size:int(exp_y)+centroid_size, int(exp_x)-centroid_size:int(exp_x)+centroid_size]

             
                    #print('expected', exp_x, exp_y)
                    #print('found', exp_x + obj_x, exp_y + obj_y)
                    
                if self.is_object(self.first_image[int(self.catalogue['ycentroid'][i])-centroid_size:int(self.catalogue['ycentroid'][i])+centroid_size, int(self.catalogue['xcentroid'][i])-centroid_size:int(self.catalogue['xcentroid'][i])+centroid_size]):
                    if not self.is_object(object_square):
                        plotx.append([x + centroid_size])
                        ploty.append([y+centroid_size])
                        
                        indexes.append(i)
                    #if len(indexes) > 10:
                     #   break
            
        output_dir = self.directory + Constants.working_directory  + Constants.output_directory + Constants.moving_obj_folder
        if not os.path.exists(output_dir):
            os.mkdir(output_dir) 
        
        moving_object_file = output_dir + Constants.moving_obj_file
        f = open(moving_object_file, "w")
        f.write("id xcentroid ycentroid")
        for i in range(len(indexes)):
            index = indexes[i]
            f.write("\r\n" + str(self.catalogue['id'][index]) + " " +  str(self.catalogue['xcentroid'][index]) + " " + str(self.catalogue['ycentroid'][index]))
            im2 = ff.get_thumbnail(Constants.moving_obj_check_image, self.catalogue['xcentroid'][indexes[i]], self.catalogue['ycentroid'][indexes[i]], 100)
            im1 = ff.get_thumbnail(1, self.catalogue['xcentroid'][indexes[i]], self.catalogue['ycentroid'][indexes[i]], 100)
            fig = plt.figure()
            fig.add_subplot(1, 2, 1)
            plt.axis('off')
            plt.imshow(im2, origin='upper', cmap=plt.cm.inferno, label = self.catalogue['id'][indexes[i]])
            #plt.scatter(plotx[i], ploty[i], s=10, c='green', marker='o')
            
            fig.add_subplot(1, 2, 2)
            plt.axis('off')

            plt.imshow(im1, origin='upper', cmap=plt.cm.inferno)
            plt.savefig(output_dir + "moving_object_id_" + str(self.catalogue['id'][indexes[i]]) + ".jpg")

        
        
        f.close()
                


                
        
        
    def get_image(self, n):
        
            set, i = Utilities.n_to_set_and_n(n)
            
            image_path = self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix + Constants.file_name
            
            file = image_path + "_" + str(set) + "_" + Utilities.format_index(i)

            file += Constants.fits_extension
        
            return fits.getdata(file, ext=0)
            
    
    def is_object(self, data):
        
        
        total_mean = np.mean(data)
        
        centre = int(len(data)/2)
        
        centre_data = data[centre-5:centre+5, centre-5:centre+5]
        
        middle_mean = np.mean(centre_data)
                
        if middle_mean > 1.2 * total_mean:
            self.means.append([total_mean, middle_mean])

            return True
        return False
        
        
        
    
        
        
        
        
                    
                    
                
                
                
                
        
        
        
            
            
            
            
    
    
    
    
    
    
    
    
    
    
    