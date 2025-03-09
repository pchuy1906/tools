import numpy as np

DFT_path = '/p/lscratchh/pham20/HN3/VASP/68-moles/HN3_room/FINAL_2_4/rdf_from_xyzf_general/'

files = ['gr_N_N.dat','gr_N_H.dat','gr_H_H.dat']

for file in files:
    DFT =  np.genfromtxt(DFT_path+file, dtype='float')
    ChIMES = np.genfromtxt(file, dtype='float')
    delta = ChIMES[:,1]-DFT[:,1]
    x = DFT[:,0]
    data = np.array([x,delta])
    data = data.T
    np.savetxt( 'delta_'+file, data)
