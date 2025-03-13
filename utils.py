import numpy as np

def read_xyz(file_xyz, cell_type):
    f  = open(file_xyz ,"r")
    
    natom = f.readline()
    natom = int(natom)
    print ("the number of atom:", natom)
    
    cell_3_3 = np.zeros(shape=(3,3))
    tmp = f.readline().split()
    
    if cell_type=="cell_3":
        for k in range(3):
            cell_3_3[k,k] = float(tmp[k])
    elif cell_type=="cell_9":
        cell_3_3[0,0] = float(tmp[0])
        cell_3_3[0,1] = float(tmp[1])
        cell_3_3[0,2] = float(tmp[2])
    
        cell_3_3[1,0] = float(tmp[3])
        cell_3_3[1,1] = float(tmp[4])
        cell_3_3[1,2] = float(tmp[5])
    
        cell_3_3[2,0] = float(tmp[6])
        cell_3_3[2,1] = float(tmp[7])
        cell_3_3[2,2] = float(tmp[8])
    elif cell_type=="NON_ORTHO":
        cell_3_3[0,0] = float(tmp[1])
        cell_3_3[0,1] = float(tmp[2])
        cell_3_3[0,2] = float(tmp[3])
    
        cell_3_3[1,0] = float(tmp[4])
        cell_3_3[1,1] = float(tmp[5])
        cell_3_3[1,2] = float(tmp[6])
    
        cell_3_3[2,0] = float(tmp[7])
        cell_3_3[2,1] = float(tmp[8])
        cell_3_3[2,2] = float(tmp[9])
    else:
        print ("unknown cell_type", cell_type)
        exit()
    
    atomList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(natom):
        tmp = f.readline()
        tmp = tmp.split()
        atomList.append(tmp[0])
        xyz[k,0] =  float(tmp[1])
        xyz[k,1] =  float(tmp[2])
        xyz[k,2] =  float(tmp[3])
    f.close
    return natom, cell_3_3, atomList, xyz


