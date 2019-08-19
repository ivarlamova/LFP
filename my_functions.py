import os, copy, math
import  numpy as np
import matplotlib as mpl
from matplotlib import cm



from siman import header
from siman.geo import xred2xcart, image_distance
from siman.calc_manage import add, res
from siman.geo import replic
db = header.db

from LFP_structures import *


def PES_scan(submit = 0, cl = None, a = None, b = None, readfiles = 1):
    """
    Create spherical PES around oxygen for H and add
    
    INPUT:
    add (bool) - 1 - run add() and submit to cluster, 0 - run res()
    
    
    RETURN:
    namelist - list of calculation names
    """

    
    
    sc = cl.end
    sc_Hrep = copy.deepcopy(sc)
    sc_Hall = copy.deepcopy(sc)

    R = sc.rprimd
    R1n = np.linalg.norm(R[0])
    R2n = np.linalg.norm(R[1])
    R3n = np.linalg.norm(R[2])

    
    # xO1 = sc.xred[18] # function for choosing the atom with (No-1) identificator
    xO1 = cl.init.xred[18] # my error, i used O position from inital
    xO1c= sc.xcart[18] 
    xO2 = sc.xred[11] # second oxide atom
    xO2c= sc.xcart[11] 
    xO3 = sc.xred[24] # third oxygen atom position
    xO3c = sc.xcart[24] 
    xO4 = sc.xred[19] # fourth oxygen atom position
    xO4c = sc.xcart[19] 
    
    
    # print('Coordinates of O1 atom:', 'x =',  xO1[0], 'y =', xO1[1], 'z =', xO1[2], '\nParemeters of uc:', 'a =', R1n, 'b =', R2n, 'c =', R3n)
    # print('Coordinates of O2 atom:', 'x =',  xO2[0], 'y =', xO2[1], 'z =', xO2[2], '\nParemeters of uc:', 'a =', R1n, 'b =', R2n, 'c =', R3n)
    # print('Coordinates of O3 atom:', 'x =',  xO3[0], 'y =', xO3[1], 'z =', xO3[2], '\nParemeters of uc:', 'a =', R1n, 'b =', R2n, 'c =', R3n)
    print('Coordinates of O4 atom:', 'x =',  xO4[0], 'y =', xO4[1], 'z =', xO4[2], '\nParemeters of uc:', 'a =', R1n, 'b =', R2n, 'c =', R3n)
    
    #For second atom
    # x_list = [num / 100 for num in range(int(round((xO2[0] - 0.3)*100)), int(round((xO2[0] + 0.4)*100)), 3)]
    # y_list = [numy / 100 for numy in range(int(round((xO2[1] - 0.3)*100)), int(round((xO2[1] + 0.4)*100)), 3)]
    # z_list = [numz / 100 for numz in range(int(round((xO2[2] - 0.3)*100)), int(round((xO2[2] + 0.4)*100)), 4)]

    #For third atom:
    # x_list = [num / 100 for num in range(int(round((xO3[0] - 0.3)*100)), int(round((xO3[0] + 0.4)*100)), 3)]
    # y_list = [numy / 100 for numy in range(int(round((xO3[1] - 0.3)*100)), int(round((xO3[1] + 0.4)*100)), 4)]
    # z_list = [numz / 100 for numz in range(int(round((xO3[2] - 0.3)*100)), int(round((xO3[2] + 0.4)*100)), 4)]

    #For fourth atom:
    x_list = [num / 100 for num in range(int(round((xO4[0] - 0.3)*100)), int(round((xO4[0] + 0.4)*100)), 3)]
    y_list = [numy / 100 for numy in range(int(round((xO4[1] - 0.3)*100)), int(round((xO4[1] + 0.4)*100)), 4)]
    z_list = [numz / 100 for numz in range(int(round((xO4[2] - 0.3)*100)), int(round((xO4[2] + 0.4)*100)), 4)]

    
    
    num=0
    namelist = []
    for i, x in enumerate(x_list):
        for j, y in enumerate(y_list):
            for k, z in enumerate(z_list):
                d=math.sqrt(math.pow((x-xO4[0])*R1n, 2)+math.pow((y-xO4[1])*R2n, 2)+math.pow((z-xO4[2])*R3n,2))
                if (d >= a and d <= b):
                    num+=1
                    xH = xred2xcart([np.array([x, y, z])],R)[0]
                    d1, d2 = image_distance(xH, xO4c, R)
    #                 print( 'Distances are: {:3.2f} {:3.2f}'.format(d1, d))
                    sc_Hijk = sc_Hrep.add_atom([x / 1, y / 1, z / 1], 'H')
                    sc_Hall = sc_Hall.add_atom([x / 1, y / 1, z / 1], 'H')
    #                if (num>50):
                    id  = (str(num)+'LFP'+"-"+str(i)+"-"+str(j)+"-"+str(k), '0', 1)
                    if submit:
                        ''
                        # sc_Hijk.printme()
                        # print(id)
                        add(*id, up = 'up2', input_st = sc_Hijk, it_folder = 'LFP_scan_O4/', cluster='pardus', corenum=4)
                    else:
                        res(*id, up = 'up2', show = 'fo', readfiles = readfiles) # try show = 'en', 'conv', 'est'
                    namelist.append(id)

    sc_Hall.name+='Hall'
    sc_Hall.write_poscar()
    # if submit:
    #     add('sc_Hall', '0', 1, up = 'up2', input_st = sc_Hall, it_folder = 'LFP_scan_O2/', cluster='pardus', corenum=4)
    return namelist


def calc_solution_energies(namelist, inum_of_H=0, num_of_H=1, cl_ref = None, shift = None):
    """
    Calculate energy of Hydrogen solution in 
    inum_of_H - number of hydrogen in initial structure, 
    num_of_H - general number of hydrogens in structure under calculation
    """
    chem_pot_0 = LFPtet1.e0 - LFP1.e0
    
    dmulist = []
    xlist = []
    
    first_calculation = db[namelist[0]]
    z = 1 #Hydrogen
    
    list_of_positions_of_atoms = first_calculation.end.get_specific_elements([z])
    
    i = list_of_positions_of_atoms[inum_of_H] # hydrogen number in lists
    print (i)
    
    for name in namelist:
        dmulist.append( db[name].e0 - cl_ref.e0 - chem_pot_0) 
        
        st = db[name].end
        if shift is not None:
            st = st.shift_atoms(shift)


        xlist.append(st.xcart[i] )
    
    st_ref = cl_ref.end.copy()
    if shift is not None:
        st_ref = st_ref.shift_atoms(shift)
    # if replic:


    return dmulist, xlist, st_ref



def energy2color(elist, cmap = None):
    
    """
    INPUT:
    elist - list of energies
    cmap (str) - color map according www.map...
    """
    
    colors = []
    max_val = max(elist)
    min_val = min(elist)


#     if cmap = 'jet':
        


    norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val) # to normalize data from 0.0 (vmin) to 1.0 (vmax)
    m = cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet) # RGBA colors to normalization
    map_to_color = np.vectorize(m.to_rgba)


    for e in elist:
        color = mpl.colors.rgb2hex(  map_to_color( e )  )
        colors.append( color )

    return colors


def modify(dmulist, xlist, namelist, mucut=15):
    """
    to deletet the points with energies more than mucut
    """
    new_dmulist = []
    new_xlist = []
    new_namelist = []
    for dmu, x, n in zip(dmulist, xlist, namelist):

        if dmu < mucut:
            new_dmulist.append(dmu)
            new_xlist.append(x)
            new_namelist.append(n)
    new_xlist_n = []
    k=0
    for n in range(len(new_xlist)):
        k+=1
        new_xlist_n = new_xlist[n]
        # print(new_xlist_n, k) 
    return new_dmulist, new_xlist, new_namelist


# def around_o(st):
#     st.write_xyz(show_around = 19, nnumber = 50, analysis = 'imp_surrounding')

def write_jmol_script(st, colors, xlist, local_path, num_of_H=0, shift = None, shift1 = None):

    
    jmolfile = 'xyz/Henergy4.jmol'
    
    
    # st.write_xyz(show_around = 12, nnumber = 50, analysis = 'imp_surrounding')['st']
       
    stH = st.add_atoms(xlist, 'H')
    # stH.nn(12, 100, from_one=0)

    if shift:
        stH = stH.shift_atoms(shift)
        # stH.nn(12, 100, from_one=0)
    
    if 0:
        filename = stH.write_poscar()
    else:
        # filename, png = stH.write_xyz(show_around = 12, nnumber = 170, analysis = 'imp_surrounding')
        filename, png = stH.write_xyz()
        stH.write_poscar()



    list_of_positions_of_atoms = stH.get_specific_elements([1])
    
    i = list_of_positions_of_atoms[num_of_H]

    print(i)
    print(st.natom)
    print (list_of_positions_of_atoms)

    basename = os.path.basename(filename)
    if local_path:
        filename = local_path + '/' +basename
    
    colors_jmol = [c.replace('#', 'x') for c in colors]
    
    with open(jmolfile, 'w') as f:
        f.write('set frank off \nset autobond off \n')
        f.write('load "'+filename+'"'+'\n')
        f.write('select H* \ncpk 30 \n')
        f.write('set displayCellParameters False ; \n' )

        for i in range(len(xlist)):
            
            f.write('select atomno = '+str(i+st.natom+1)+'\n')
            # f.write('select atomno = '+str(i+2)+'\n')
            
            f.write('color ['+colors_jmol[i]+'];\n' )

        xlist = [stH.xcart[i] for i in list_of_positions_of_atoms]
        xc = sum(xlist)/len(xlist)
        print(xc)

        stH = stH.copy()
        # stH= replic(stH, (2, 2, 2))
        if shift1 is not None:
            stH = stH.shift_atoms(shift1)
        # stH = stH.return_atoms_to_cell()

        for i, x in enumerate(stH.xcart):
            d1, d2 = image_distance(xc, x, stH.rprimd)
            if d1>7:
                f.write('select atomno = '+str(i+1)+'; cpk 0\n')



    return jmolfile


def points_relax (xlist, st, inum_of_H=0, submit = 0, namelist=None):

    stH = st.add_atoms(xlist, 'H')

    list_of_positions_of_atoms = stH.get_specific_elements([1])
    
    i = list_of_positions_of_atoms[inum_of_H] 

    print(i)
    
    # to_calc1 = [119, 129, 131, 140, 146, 179, 196, 202, 207, 244, 253] # distribution around first-oxide
    
    # to_calc2 = [142, 131, 175, 148, 119, 163, 160, 129, 137, 183, 166, 155, 120] # distribution around second-oxide
     
    
    print (st.natom)
    n=[]
    energy=[]
    energy_diff=[]
    namelist_second=[]

    ''' # distribution around third-oxide

    to_calc3 = [556, 593, 608, 530, 641, 607, 578, 586, 570, 547, 680, 613] # in 3if SP calc
    for k in to_calc3:
        j1=k-st.natom-1
        print(j1)
        cl1=db[namelist[j1]]

        if submit:
            cl = db[cl1.id[0]+'.if'+'.if', '1l', 1]
            # cl.res()
            cl.run('0', iopt='full', i_child = 1, it_folder = cl.sfolder, add=1, up='up2')

        else:
        #     print (i, st.natom, j)

            cl1=db[cl1.id[0]+'.if'+'.if'+'.if', '0', 1]
            # cl1.res(up = 'up2')
            energy.append(cl1.e0)

    to_calc3a = [569, 685] # in 5if SP calc
    for l in to_calc3a:
        j2=l-st.natom-1
        print(j2)
        cl2=db[namelist[j2]]
        if submit:
            cl = db[cl1.id[0]+'.if'+'.if', '1l', 1]
            cl.run('0', iopt='full', i_child = 1, it_folder = cl.sfolder, add=1, up='up2')
        else:
            cl2=db[cl2.id[0]+'.if'+'.if'+'.if'+'.if'+'.if', '0', 1]
            # cl2.res(up = 'up2')
            energy.append(cl2.e0)
    '''
    '''
    to_calc4 = [557, 672, 574, 622, 567, 591, 684, 630, 565, 588, 680, 627] # distribution around fourth-oxide
    to_calc17 = [557] # 5if SP
    to_calc132 = [672] # 3if SP
    to_calc34 = [574] # 5if SP
    to_calc82 = [622] # 4if SP
    to_calc27 = [567, 591, 684] # 4if SP
    to_calc51 = [591] # 4if SP
    to_calc144 = [684] # 4if SP
    to_calc90 = [630] # 3if SP
    to_calc25 = [565] # 2if SP
    to_calc48 = [588] # 3if SP
    to_calc140 = [680] # 4if SP
    to_calc87 = [627] # 4if SP
    '''
    to_calc4_4 = [622, 680, 627, 567, 591, 684]
    to_calc4_3 = [672, 630, 588]
    to_calc4_5 = [557, 574]
    to_calc4_2 = [565]
    to_calc132 = [672] # 3if SP


    for k in to_calc4_5 :
        j1=k-st.natom-1
        # print(j1)
        cl1=db[namelist[j1]]
        # print(namelist)
        print(st.natom, k, j1, cl1)

        if submit:
            cl = db[cl1.id[0]+'.if'+'.if', '1l', 1]
            # cl.res(up = 'up2')

            cl.run('1l', iopt='full', i_child = None, it_folder = cl.sfolder, add=1, up='up2', cluster = 'cee')
            # cl2 = db[cl1.id[0]+'.if.if', '1l', 1]
            # cl2.res(up = 'up2')
        else:
        #     print (i, st.natom, j)

            cl=db[cl1.id[0]+'.if'+'.if'+'.if'+'.if'+'.if', '0', 1]
            cl.res(up='up2')
            # cl.res(show = 'en')
            energy.append(cl.e0)

    for k in to_calc4_4 :
        j1=k-st.natom-1
        # print(j1)
        cl1=db[namelist[j1]]
        # print(namelist)
        print(st.natom, k, j1, cl1)

        if submit:
            cl = db[cl1.id[0]+'.if'+'.if', '1l', 1]
            # cl.res(up = 'up2')

            cl.run('1l', iopt='full', i_child = None, it_folder = cl.sfolder, add=1, up='up2', cluster = 'cee')
            # cl2 = db[cl1.id[0]+'.if.if', '1l', 1]
            # cl2.res(up = 'up2')
        else:
        #     print (i, st.natom, j)
            cl=db[cl1.id[0]+'.if'+'.if'+'.if'+'.if', '0', 1]
            cl.res(up='up2')
            # cl.res(show = 'en')
            energy.append(cl.e0)

    for k in to_calc4_3 :
        j1=k-st.natom-1
        # print(j1)
        cl1=db[namelist[j1]]
        # print(namelist)
        print(st.natom, k, j1, cl1)

        if submit:
            cl = db[cl1.id[0]+'.if'+'.if', '1l', 1]
            # cl.res(up = 'up2')

            cl.run('1l', iopt='full', i_child = None, it_folder = cl.sfolder, add=1, up='up2', cluster = 'cee')
            # cl2 = db[cl1.id[0]+'.if.if', '1l', 1]
            # cl2.res(up = 'up2')
        else:
        #     print (i, st.natom, j)
            cl=db[cl1.id[0]+'.if'+'.if'+'.if', '0', 1]
            cl.res(up='up2')
            # cl.res(show = 'en')
            energy.append(cl.e0)

    for k in to_calc4_2 :
        j1=k-st.natom-1
        # print(j1)
        cl1=db[namelist[j1]]
        # print(namelist)
        print(st.natom, k, j1, cl1)

        if submit:
            cl = db[cl1.id[0]+'.if'+'.if', '1l', 1]
            # cl.res(up = 'up2')

            cl.run('1l', iopt='full', i_child = None, it_folder = cl.sfolder, add=1, up='up2', cluster = 'cee')
            # cl2 = db[cl1.id[0]+'.if.if', '1l', 1]
            # cl2.res(up = 'up2')
        else:
        #     print (i, st.natom, j)
            cl=db[cl1.id[0]+'.if'+'.if'+'.ifc', '0', 1]
            # cl1=db[cl1.id[0]+'.if'+'.if', '1l', 1]
            cl.res(up='up2')
            # cl.res(show = 'en')
            energy.append(cl.e0)

    print('Energy array:\n', energy)
    emin=min(energy)
    print('Minimal energy:\n', emin)
    for k in energy:
        dif=k-emin
        energy_diff.append(dif)
    print('Number of atom in system:\n', to_calc4_5, to_calc4_4, to_calc4_3, to_calc4_2)
    print('Energy difference:\n', energy_diff)
     
    # for i in to_calc2:
    #     j2=i-st.natom-1
    #     cl2=db[namelist[j2]]

    # for i in to_calc3:
    #     j3=i-st.natom-1
        # cl3=db[namelist[j3]]
     
    #     id_second  = (str(namelist[j2])+'.if', '1l', 1)
        # namelist_second.append(id_second)
        # cl_second=db[namelist_second[0]]
        # print(id_second)
        # if j==10:
        #     print(namelist[j])
        #     # ""
        #     cl.end.jmol()
        # id  = (str(j1)+'LFP'+'.if', '1l', 1)
        # if submit:
        #     cl = db[cl1.id[0]+'.if'+'.if', '1l', 1]
        #     # cl.res()
        #     cl.run('0', iopt='full', i_child = 1, it_folder = cl.sfolder, add=1, up='up2')
            
            
            # print(cl1.id[0]+'.if'+'.if', '1l', 1)
            # cl2.res()
            # cl=db[cl1.id[0]+'.if', '1l', 1]
            # cl.run('1l', iopt='full', it_folder = cl1.sfolder, add=0, up='up2')
            # cl2.run('1l', iopt='full', it_folder = cl2.sfolder, i_child = 1, add=0, up='up2')

            # add(*id, up = 'up2', input_st = stH, it_folder = 'LFP_scan_O2', cluster='pardus', corenum=4)
        # else:
        # #     print (i, st.natom, j)
        #     cl1=db[cl1.id[0]+'.if'+'.if'+'.if', '0', 1]
        #     cl1.res(up = 'up2')

            # cl2=db[cl2.id[0]+'.if'+'.if'+'.if'+'.if'+'.if', '0', 1]
            # cl.run('1l', iopt='full', i_child = 1, it_folder = cl1.sfolder, add=0, up='up2')

            # cl = db[cl1.id[0]+'.if', '1l', 1]

            # cl2 = db[cl2.id[0]+'.if'+'.if'+'.if'+'.if'+'.if', '0', 1]
            # cl2.res(up = 'up2')

            # print(cl.res)
            # cl.run('0', iopt='full', it_folder = cl1.sfolder, add=1, up='up2')
            # cl.res(comment=1, up = 'up2', show = 'path') # try show = 'en', 'conv', 'est'
            
            # cl2.res(up = 'up2')
            # print(cl.end.magmom)
            
            # print(cl2.end.magmom)
            # namelist.append(id)
            # energy.append(cl.e0)
    # print('Energy array:\n', energy)
    # emin=min(energy)
    # print('Minimal energy:\n', emin)
    # for k in energy:
    #     dif=k-emin
    #     energy_diff.append(dif)
    # print('Number of atom in system:\n', to_calc2)
    # print('Energy difference:\n', energy_diff)
    return n

def separate_check():
    # db['101LFP-11-5-7.if.if.1l.1'].end.get_mag_tran(to_ox=-2)
    # db['101LFP-11-5-7.if.if.if.0.1'].end.get_mag_tran(to_ox=-2)
    # print (db['101LFP-11-5-7.if.if.if.0.1'].end.magmom)
    print (db['25LFP-7-10-8.if.if.1l.1'].end.magmom)
    print (db['25LFP-7-10-8.if.if.if.0.1'].end.magmom)





# 