import numpy as np


import argparse

parser = argparse.ArgumentParser(description='make ChIMES parameter file')

# Arguments supported by the code.
parser.add_argument("--file_input",          default='dlars.log', help='file dlars.log')
parser.add_argument("--iteration", type=int, default=1000,        help='interation')
parser.add_argument("--nParaMax",  type=int, default=1000,        help='nParaMax')

args        = parser.parse_args()
file_input   = args.file_input
iteration    = args.iteration
nParaMax     = args.nParaMax

print ("The number of variables ", nParaMax)

DoneRead = False
with open( file_input) as f:
    while not DoneRead:
        content = f.readline()
        if "RMS" in content:
            RMS = content
            #print(content)
        if "[" in content:
            params = np.zeros(nParaMax)
            SubReadingParams = True
            while SubReadingParams:
                content = f.readline().split()
                if (len(content) == 2):
                    id = int(float(content[0]))
                    params[id] = float(content[1])
                    #print (id, content)
                else:
                    SubReadingParams = False
                    #np.savetxt('x0.dat', params, fmt='%15.9f')

        if "Finished iteration" in content:
            #print(content)
            it = int(float(content.split()[-1]))
            #print (it)
            if (it == iteration):
                print (RMS)
                print (content)
                np.savetxt('x0.dat', params, fmt='%15.9f')
                break

        if not content:
            break
