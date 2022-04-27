# =============================================================================
# import requests
# import json
# import Constants
# from astropy.io.fits import HDUList
# from astropy.io.fits import PrimaryHDU
# from astropy.io.fits import getdata
# from astropy.io.fits import getheader
# import os
# from astropy.table import Table
# from astroquery.astrometry_net import AstrometryNet
# 
# 
# 
#     filesdir = None
#     file_name = None
#     session_key = None
#     
#     def __init__(self, sources):
#         self.filesdir = filesdir
#         self.file = file
#     
# # =============================================================================
# #     def connect(self):
# #     
# #         R = requests.post('http://nova.astrometry.net/api/login', data={'request-json': json.dumps({"apikey": Constants.api_key})})
# #         print(R.text)
# #         session_key = R.text.split(":")[4]
# #         self.session_key = session_key[2:len(session_key)-2]
# #         print(session_key)
# #     
# #     def create_file(self):
# #         
# #         form = ""
# #         NL = "\r\n"
# #                 
# #         form = form + "--===============2521702492343980833==" + NL
# #         form = form + "Content-Type: text/plain" + NL
# #         form = form + "MIME-Version: 1.0" + NL
# #         form = form + "Content-disposition: form-data; name=\"request-json\"" + NL
# #         form = form + NL
# #         form = form + "{\"publicly_visible\": \"y\", \"allow_modifications\": \"d\", \"session\": \"" + self.session_key + "\", \"allow_commercial_use\": \"d\"}"+ NL
# #         form = form + "--===============2521702492343980833=="+ NL
# #         form = form + "Content-Type: application/octet-stream"+ NL
# #         form = form + "MIME-Version: 1.0"+ NL
# #         form = form + "Content-disposition: form-data; name=\"file\"; filename=\"" + Constants.job_file_name + "\""+ NL
# #         form = form + NL
# #                         
# #         data = open(self.filesdir + "r_l141_5_1_001.fit", 'r')
# #         
# #         line = data.readline()
# #         
# #         while line:
# #             form = form + line
# #             line = data.readline()
# #         
# #         data.close
# #             
# #         form = form + NL 
# #         form = form + NL
# #         form = form + "--===============2521702492343980833==--"
# #     
# #         print(form)
# #     
# #     def getWCS(self):
# #         self.connect()
# #         self.create_file()
# # =============================================================================
#     
# 
#     def solve(self):
#         ast = AstrometryNet()
#         ast.api_key = 'XXXXXXXXXXXXXXXX'
# 
#         sources = Table.read('catalog.fits')
#         # Sort sources in ascending order
#         sources.sort('FLUX')
#         # Reverse to get descending order
#         sources.reverse()
#         
#         image_width = 3073
#         image_height = 2048
#         wcs_header = ast.solve_from_source_list(sources['X_IMAGE'], sources['Y_IMAGE'],
#                                                 image_width, image_height,
#                                                 solve_timeout=120)
#         
# 
# =============================================================================
