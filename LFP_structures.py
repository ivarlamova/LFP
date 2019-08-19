from siman import header
from siman.calc_manage   import smart_structure_read
db = header.db
LFPsPvac= db['LiFePO4_Pvac_in1uc.if.1.1'] # LiFePO4 unit cell with P vacancy after full relaxation
LFP1 = db['LiFePO4122.1.1'] # LFP 122 after relaxation
LFPtet1 = db['LiFePO4122tet1.1.1'] # LFP 122 +H in free tetrahedral after relaxation
LFPO1min=db['101LFP-11-5-7.if.if.if.if.0.1'] # LFP structure with P vac and one of H (near 19 O) in relevant position after relax
LFPO2min=db['19LFP-8-7-5.if.if.if.if.0.1'] # LFP structure with P vac and two of H (near 19 and 12 O) in relevant position after relax
LFPO3min=db['71LFP-9-10-11.if.if.if.0.1'] # LFP structure with P vac and three of H (near 19 and 12 and 25 O) in relevant position after relax
# LFPsPvac_19 = smart_structure_read(input_geo_file ='xyz/LiFePO4_Pvac_in1uc.if.1.1.end_loc19.xyz')
# LFPsPvac_19_H = smart_structure_read(input_geo_file ='xyz/POSCAR_LiFePO4_Pvac_in1uc_if_1_1_end_loc19')
#LFPsPvac_pos = smart_structure_read (input_geo_file = 'xyz/POSCAR_LiFePO4_Pvac_in1uc_if_1_1_end')

#LFP_tet1 = smart_structure_read(input_geo_file = 'LFP_tet1/POSCAR') # POSCAR of LFP uc with H4 in free tetrahedral position (more stable)
