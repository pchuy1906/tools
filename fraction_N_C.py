import numpy as np
 
numCHNO = np.loadtxt('CHNO.dat')

rows = len(numCHNO)
print (rows)

half_rows = rows//2
print (half_rows)

numC = numCHNO[half_rows:rows,1]
numH = numCHNO[half_rows:rows,2]
numN = numCHNO[half_rows:rows,3]
numO = numCHNO[half_rows:rows,4]

numTot = numC + numH + numN + numO

fracC = numC/numTot
fracH = numH/numTot
fracN = numN/numTot
fracO = numO/numTot

aveC = np.mean(fracC)
aveH = np.mean(fracH)
aveN = np.mean(fracN)
aveO = np.mean(fracO)

print (aveC, aveH, aveN, aveO)
