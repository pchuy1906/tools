import numpy as np
import os

element_dict_spdf = {\
'C' : "p",\
'H' : "s",\
'N' : "p",\
'O' : "p",\
'Al' : "p",\
'Fe' : "d" \
}

import argparse
parser = argparse.ArgumentParser(description='setup DFTB input for file_xyz')
# Arguments supported by the code.
parser.add_argument("--file_xyz",                          default='file.xyz',   help='file_xyz format xyz')
parser.add_argument("--cell_option",                       default='cell_3',     help='cell_option = "NON_ORTHO", "cell_3", or "cell_9" ')
# Driver
parser.add_argument("--calc_type",                         default="sp",         help='sp, opt, nvt')
parser.add_argument("--mdTemp",                type=float, default=300,          help='MD temperature')
# Hamiltonian
parser.add_argument("--eTemp",                 type=float, default=300,          help='electronic temperature')
parser.add_argument("--SKF",                               default="pbc-0-3",    help='pbc-0-3 / 3ob')
parser.add_argument("--ChIMES_params_file",                default="",           help='ChIMES.txt')
parser.add_argument("--PolyRep",                           default="Yes",        help='Yes (Polynomial repulsive from SKF)/No (Repulsive from Spline)')
parser.add_argument("--Kresol",                type=float, default=0.08,         help='Resolution for K-point grid ')
parser.add_argument("--Kpoint",     nargs='+',                                   help='K-point grid ')
parser.add_argument("--method",                            default="DFTB",       help='DFTB, DFTB_ChIMES')
parser.add_argument('--atom_types', nargs='+',                                   help='atom types: C H N O')
parser.add_argument('--spdf_types', nargs='+',                                   help='atom types: s p d f')
parser.add_argument("--dispersion",                        default="",           help='PBE-BJ')


args        = parser.parse_args()
file_xyz           = args.file_xyz
cell_option        = args.cell_option
Kresol             = args.Kresol
Kpoint             = args.Kpoint
method             = args.method
calc_type          = args.calc_type
eTemp              = args.eTemp
mdTemp             = args.mdTemp
SKF                = args.SKF
ChIMES_params_file = args.ChIMES_params_file
PolyRep            = args.PolyRep
atom_types         = args.atom_types
spdf_types         = args.spdf_types
dispersion         = args.dispersion

def cell_9_to_abc(cell_9):
    vec1= cell_9[0:3]
    vec2= cell_9[3:6]
    vec3= cell_9[6:9]
    a = np.sqrt(np.dot(vec1, vec1))
    b = np.sqrt(np.dot(vec2, vec2))
    c = np.sqrt(np.dot(vec3, vec3))
    alp = np.arccos(np.dot(vec2, vec3)/(b*c))
    bet = np.arccos(np.dot(vec1, vec3)/(a*c))
    gam = np.arccos(np.dot(vec1, vec2)/(a*b))
    return np.array([a,b,c,alp,bet,gam])
 
def KP(cell_9, Kresol):
    cell_xyz = np.zeros(shape=(3,3))
    cell_xyz[0,:]= cell_9[0:3]
    cell_xyz[1,:]= cell_9[3:6]
    cell_xyz[2,:]= cell_9[6:9]

    vol = np.linalg.det(cell_xyz)

    cell_abc = cell_9_to_abc(cell_9)
 
    dist = np.zeros(shape=(3))
    dist[2] = vol/abs(cell_abc[0]*cell_abc[1]*np.sin(cell_abc[5]))
    dist[1] = vol/abs(cell_abc[0]*cell_abc[2]*np.sin(cell_abc[4]))
    dist[0] = vol/abs(cell_abc[1]*cell_abc[2]*np.sin(cell_abc[3]))

    #Kpoints = np.zeros(shape=(3))
    #for k in range(3):
    #    Kpoints[k] = np.ceil(1.0/(dist[k]* float(Kresol)))
    Kpoints = np.ceil(1.0/(dist* float(Kresol))).astype(int)
    return Kpoints

def write_dftb_input(method, calc_type, eTemp, mdTemp, SKF, ChIMES_params_file, PolyRep, Kpoints, syms, spdf, dispersion):
    f2 = open("dftb_in.hsd", "w")
    f2.write("%-s\n" %( "Parallel { UseOmpThreads = Yes}"))
    f2.write("%-s\n" %( ""))
    f2.write("%-s\n" %( "Geometry = GenFormat {"))
    f2.write("%-s\n" %( "    <<<\"dftb.gen\""))
    f2.write("%-s\n" %( "}"))
    f2.write("%-s\n" %( ""))

    if (calc_type=="sp"):
        f2.write("%-s\n" %( "Driver = GeometryOptimization {"))
        f2.write("%-s\n" %( "    Optimizer = Rational {}"))
        f2.write("%-s\n" %( "    MaxSteps = 0"))
        f2.write("%-s\n" %( "}"))
    elif (calc_type=="opt"):
        f2.write("%-s\n" %( "Driver = GeometryOptimization {"))
        f2.write("%-s\n" %( "    Optimizer = Rational {}"))
        f2.write("%-s\n" %( "    MaxSteps = 10000"))
        f2.write("%-s\n" %( "}"))
    elif (calc_type=="vc-opt"):
        f2.write("%-s\n" %( "Driver = GeometryOptimization {"))
        f2.write("%-s\n" %( "    Optimizer = Rational {}"))
        f2.write("%-s\n" %( "    MaxSteps = 10000"))
        f2.write("%-s\n" %( "    LatticeOpt = Yes"))
        f2.write("%-s\n" %( "}"))
    elif (calc_type=="nvt"):
        f2.write("%-s\n"        %( "Driver = VelocityVerlet {"))
        f2.write("%-s\n"        %( "    Steps = 2000000"))
        f2.write("%-s\n"        %( "    TimeStep [Femtosecond] = 0.1"))
        f2.write("%-s\n"        %( "    MDRestartFrequency = 1"))
        f2.write("%-s\n"        %( "    Thermostat = NoseHoover {"))
        f2.write("%-s\n"        %( "        CouplingStrength [cm^-1] = 3500"))
        f2.write("%-s %12.3f\n" %( "        Temperature [Kelvin] = ", mdTemp))
        f2.write("%-s\n"        %( "    }"))
        f2.write("%-s\n"        %( "    MovedAtoms = 1:-1"))
        f2.write("%-s\n"        %( "}"))


    f2.write("%-s\n" %( ""))
    f2.write("%-s\n" %( "Options = {}"))
    f2.write("%-s\n" %( ""))
    f2.write("%-s\n" %( "ParserOptions = { ParserVersion = 3}"))
    f2.write("%-s\n" %( ""))

    if (method=="DFTB"):
        f2.write("%-s\n" %( "Hamiltonian = DFTB {"))
        f2.write("%-s\n" %( "    #ReadInitialCharges = Yes"))
        f2.write("%-s\n" %( "    #Solver = ELPA {"))
        f2.write("%-s\n" %( "    #    Mode = 2"))
        f2.write("%-s\n" %( "    #}"))

        f2.write("%-s\n" %( "    SlaterKosterFiles = Type2FileNames {"))
        f2.write("%-s\n" %( "        Prefix = \"" + SKF + "/\""))
        f2.write("%-s\n" %( "        Separator = \"-\""))
        f2.write("%-s\n" %( "        Suffix = \".skf\""))
        f2.write("%-s\n" %( "        LowerCaseTypeName = No"))
        f2.write("%-s\n" %( "    }"))
        f2.write("%-s\n" %( "    MaxAngularMomentum = {"))
        for i in range(len(spdf)):
            f2.write("%-s\n" %( "        "+syms[i]+" = \""+spdf[i]+"\""))
        f2.write("%-s\n" %( "    }"))
    
        if (PolyRep=="Yes"):
            f2.write("%-s\n" %( "    PolynomialRepulsive = {"))
            for i in range(len(syms)):
                for j in range(len(syms)):
                    f2.write("%-s\n" %( "        "+syms[i]+"-"+syms[j]+" = Yes"))
            f2.write("%-s\n" %( "    }"))
    elif (method=="GFN1-xTB"):
        f2.write("%-s\n" %( "Hamiltonian = xTB {"))
        f2.write("%-s\n" %( "    #ReadInitialCharges = Yes"))
        f2.write("%-s\n" %( "    #Solver = ELPA {"))
        f2.write("%-s\n" %( "    #    Mode = 2"))
        f2.write("%-s\n" %( "    #}"))
        f2.write("%-s\n" %( "    Method = \""+method+"\""))




    if ChIMES_params_file:
        print ("ChIMES")
        f2.write("%-s\n" %("    ChIMES {"))
        f2.write("%-s\n" %("        ParameterFile = \"" + ChIMES_params_file +"\""))
        f2.write("%-s\n" %("    }"))


    f2.write("%-s\n" %( "    SCC = YES"))
    f2.write("%-s\n" %( "    MaxSCCIterations = 200"))
    f2.write("%-s\n" %( "    Mixer = Broyden {"))
    f2.write("%-s\n" %( "        MixingParameter = 0.1"))
    f2.write("%-s\n" %( "    }"))

    if (eTemp > 10.0):
        f2.write("%-s\n"        %( "    Filling = Fermi {"))
        f2.write("%-s %12.3f\n" %( "        Temperature [Kelvin] = ", eTemp))
        f2.write("%-s\n"        %( "    }"))

    f2.write("%-s\n" %( "    KPointsAndWeights = SupercellFolding {"))
    f2.write("%-s\n" %( "        " + str(Kpoints[0]) + "  0  0"))
    f2.write("%-s\n" %( "        0  " + str(Kpoints[1]) + "  0"))
    f2.write("%-s\n" %( "        0  0  " + str(Kpoints[2])))
    f2.write("%-s\n" %( "        0  0  0"))
    f2.write("%-s\n" %( "    }"))

    if (dispersion=="PBE-BJ"):
        f2.write("%-s\n" %( "    Dispersion = DftD3 {"))
        f2.write("%-s\n" %( "        Damping = BeckeJohnson {"))
        f2.write("%-s\n" %( "            a1 = 0.4289"))
        f2.write("%-s\n" %( "            a2 = 4.4407"))
        f2.write("%-s\n" %( "        }"))
        f2.write("%-s\n" %( "        s6 = 1.0"))
        f2.write("%-s\n" %( "        s8 = 0.7875"))
        f2.write("%-s\n" %( "    }"))

    f2.write("%-s\n" %( "}"))
    f2.write("%-s\n" %( " "))
    f2.write("%-s\n" %( " "))
    f2.write("%-s\n" %( "Analysis {"))
    f2.write("%-s\n" %( "    CalculateForces = Yes"))
    f2.write("%-s\n" %( "}"))
    f2.close

f  = open(file_xyz ,"rt")
nframe = 0
while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break
    #v, e = [float(x) for x in line.split()[:2]]

    # line-1: number of atoms
    natom = int(tmp)
    
    # line-2: cell_parameter, stress tensor and energy
    tmp = f.readline().split()

    if (cell_option == "NON_ORTHO"):
        cell_9 = [float(x) for x in tmp[1:10]]
        stress_6 = [float(x) for x in tmp[10:16]]
        #print ("stress tensors ", stress_6)
        energy = float(tmp[16])
        #print ("energy ", energy)
    elif (cell_option == "cell_9"):
        cell_9 = [float(x) for x in tmp[0:9]]
    elif (cell_option == "cell_3"):
        cell_3 = [float(x) for x in tmp[0:3]]
        cell_9 = [0.0 for x in range(9)]
        cell_9[0] = cell_3[0]
        cell_9[4] = cell_3[1]
        cell_9[8] = cell_3[2]
    else:
        print ("Not implemented yet!")
        exit()
    #print ("cell parameters ", cell_9)

    # line-next: xyz coordinate
    atomList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(natom):
        tmp = f.readline().split()
        atomList.append(tmp[0])
        xyz[k,0],xyz[k,1],xyz[k,2] =  float(tmp[1]),float(tmp[2]),float(tmp[3])
    nframe += 1

    # print to dftb input file
    syms, counts_syms = np.unique(atomList, return_counts=True)
    asym_list = ' '.join(syms)
    spdf = []
    for asym in syms:
        spdf.append(element_dict_spdf.get(asym))
    #print (syms, spdf)

    f2 = open("dftb.gen", "w")
    f2.write("%-d %4s\n" %( natom, " S" ))
    f2.write("%-s\n" %(asym_list))
    for k in range(natom):
        ind = np.where(syms == atomList[k])[0][0]+1
        f2.write("%-5d %5d %15.9f %15.9f %15.9f\n" %(k+1, ind, xyz[k,0], xyz[k,1], xyz[k,2] ))
    f2.write("%15.9f %15.9f %15.9f\n" %( 0.0, 0.0, 0.0 ))
    f2.write("%15.9f %15.9f %15.9f\n" %( cell_9[0], cell_9[1], cell_9[2] ))
    f2.write("%15.9f %15.9f %15.9f\n" %( cell_9[3], cell_9[4], cell_9[5] ))
    f2.write("%15.9f %15.9f %15.9f\n" %( cell_9[6], cell_9[7], cell_9[8] ))
    f2.close

    if Kpoint is None:
        Kpoints = KP(cell_9, Kresol)
        print (type(Kpoints))
    else:
        Kpoints = Kpoint
    write_dftb_input(method, calc_type, eTemp, mdTemp, SKF, ChIMES_params_file, PolyRep, Kpoints, syms, spdf, dispersion)

    fold="runDFTB/run_"+str(nframe)
    cmd = "mkdir -p " + fold
    os.system(cmd)

    cmd = "mv dftb.gen dftb_in.hsd " + fold
    os.system(cmd)

f.close

