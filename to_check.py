import os 
import  numpy as np
import matplotlib as mpl
from matplotlib import cm



from siman import header
db = header.db

# from my_functions import *
from LFP_structures import *

LFPsPvac1 = LFPsPvac.end.write_xyz(show_around = 18, nnumber = 50, analysis = 'imp_surrounding') #to cut the exiting uc

# print(LFPsPvac1)

# from my_functions import PES_scan

# db['HPES1'] = PES_scan(cl=LFPsPvac)
# print (db['153LFP-14-5-4.0.1'].end)



