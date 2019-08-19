import os, copy
import  numpy as np
import matplotlib as mpl
from matplotlib import cm



from siman import header
db = header.db

from my_functions import PES_scan, calc_solution_energies, modify, energy2color, write_jmol_script, points_relax, separate_check
# import to_check

from LFP_structures import LFPsPvac, LFPO1min, LFPO2min, LFPO3min
from siman.geo import replic

''' # around second oxigen
LFPO1min.end.write_xyz(show_around = 12, nnumber = 50, analysis = 'imp_surrounding') #to cut the exiting uc
# print(LFPsPvac1)
# db['HPES_O2'] = PES_scan(cl=LFPO1min, submit=0, a=0.93, b=1, readfiles = 1)
dmulist, xlist, st_ref =  calc_solution_energies(db['HPES_O2'], inum_of_H=1, num_of_H=2, cl_ref = LFPO1min, shift = (0.3, 0.3, -0.3))
# print(dmulist, xlist)
print(len(dmulist))
dmulist, xlist, namelist = modify(dmulist, xlist, db['HPES_O2'], mucut=4) #
print(len(dmulist))
# print (namelist, xlist)
colors = energy2color(dmulist) #
st_ref = replic(st_ref, (2, 3, 3))
# st = replic(LFPO1min.end)
# st.nn(12, 100, from_one=0) # to print 100 number of neighbors of 12 atom
write_jmol_script(st_ref, colors, xlist, local_path = '/home/irina/LFP/xyz', num_of_H=1, shift = (0.3, 0.3, 0.3))
# points_relax (xlist, st=st_ref, inum_of_H=1, submit = 0, namelist=namelist)
# print(xlist, dmulist, colors)
'''

''' # around third oxigen
# db['HPES_O3'] = PES_scan(cl=LFPO2min, submit=0, a=0.91, b=1.1, readfiles = 1)
dmulist, xlist, st_ref =  calc_solution_energies(db['HPES_O3'], inum_of_H=2, num_of_H=3, cl_ref = LFPO2min, shift = (0.0, 0.0, 0.0)) #shift = (-0.4, 0.4, -0.1)
# print(len(dmulist))
# print(dmulist, xlist)
dmulist, xlist, namelist = modify(dmulist, xlist, db['HPES_O3'], mucut=10) 
colors = energy2color(dmulist)
# print(len(dmulist))
# print (namelist, xlist)
# print(dmulist)
st = replic(LFPO2min.end)
st_ref = replic(st_ref, (2, 3, 3)) # include_boundary = (1,5)
# st.nn(12, 100, from_one=0) # to print 100 number of neighbors of 12 atom
write_jmol_script(st_ref, colors, xlist, local_path = '/home/irina/LFP/xyz', num_of_H=2, shift = (0.3, 0.3, 0.3), shift1 = (0.0, 0.0, 0.0))
# points_relax (xlist, st=st_ref, inum_of_H=2, submit = 0, namelist=namelist)
# print(xlist, dmulist, colors)
'''

# around fourth oxigen
# db['HPES_O4'] = PES_scan(cl=LFPO3min, submit=0, a=0.91, b=1.1, readfiles = 1)
dmulist, xlist, st_ref =  calc_solution_energies(db['HPES_O4'], inum_of_H=3, num_of_H=4, cl_ref = LFPO3min, shift = (0.3, 0.3, 0.0)) #shift = (-0.4, 0.4, -0.1)
# print(len(dmulist))
# print(dmulist, xlist)
dmulist, xlist, namelist = modify(dmulist, xlist, db['HPES_O4'], mucut=100) 
colors = energy2color(dmulist)
# print(len(dmulist))
# print (namelist, xlist)
# print(dmulist)
# st = replic(LFPO3min.end)
st_ref = replic(st_ref, (2, 3, 3)) # include_boundary = (1,5)
# st.nn(12, 100, from_one=0) # to print 100 number of neighbors of 12 atom
# write_jmol_script(st_ref, colors, xlist, local_path = '/home/irina/LFP/xyz', num_of_H=3, shift = (0.3, 0.3, 0.3), shift1 = (0, 0, 0.0))
points_relax (xlist, st=st_ref, inum_of_H=3, submit = 0, namelist=namelist)
# print(xlist, dmulist, colors)









# separate_check

# PES_relax(xlist, submit = 1, cl=LFPsPvac.end)

# st = replic(LFPsPvac.end, (1, 2, 2))
# shift =(-0.4, 0.4, 0.40)

