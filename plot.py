import matplotlib.pyplot as plt
import numpy as np
import sys
from itertools import cycle

openfile = sys.argv[1]
filedirectory_to = sys.argv[2]
filedirectory = sys.argv[3]

FIELDfilenam = 'FIELD'

check1 = False
check2 = False
molecule_index = 0
MOLECULE_list = {}
order = {}
natmm_list = {}
total_number = 0
with open(filedirectory + FIELDfilenam,'r') as f:
    for i, line in enumerate(f):
        if not line.strip():
            continue
        else:
            if line.split()[0] == 'MOLECULES' or line.split()[0] == 'molecular':
                ntypemolecule = int(line.split()[-1])
                check1 = True
            if check1 and line.split()[0] != 'finish':
                if line.split()[0] == 'nummols':
                    MOLECULE_list[str(molecule_index+1)] = {'number':int(line.split()[-1])}
                    total_number += int(line.split()[-1])
                elif line.split()[0] == 'atoms':
                    natmm = int(line.split()[-1])
                    natmm_list[str(molecule_index+1)] = natmm
                    check2 = True
                    j = 0
                    continue
                if check2 and j < natmm:
                    MOLECULE_list[str(molecule_index+1)][line.split()[0]] = int(line.split()[3])
                    j += int(line.split()[3])
                    if 'O' in line.split()[0]:
                        order[str(molecule_index+1)] = j-1
            if check1 and line.split()[0] == 'finish':
                molecule_index += 1
                if molecule_index == ntypemolecule:
                    break
                else:
                    check2 = False


linestyles = ['-']

hbond = {}
check = False
with open(openfile,'r') as f:
	for i, line in enumerate(f):
		if line.split(':')[0] == 'pair':
			pair = line.split('\n')[0].split(':')[1]
			hbond[pair] = []
			check = True
		else:
			if check == False or not line.strip():
				continue
			hbond[pair].append(float(line.split()[0]))

for pair in hbond:
	hbond[pair] = np.array(hbond[pair])

eta0 = (hbond['0-0']/(hbond['0-0']+hbond['0-1']))*(float(total_number)/MOLECULE_list['1']['number'])
eta1 = (hbond['1-1']/(hbond['1-1']+hbond['0-1']))*(float(total_number)/MOLECULE_list['2']['number'])



eta0 = np.cumsum(eta0)/np.arange(1,len(eta0)+1,dtype = float)
eta1 = np.cumsum(eta1)/np.arange(1,len(eta1)+1,dtype = float)

eta = {r'$\eta_0$':eta0,r'$\eta_1$':eta1}

fig, ax = plt.subplots()

for item in iter(eta):
	linecycler = cycle(linestyles)
	ax.plot(np.arange(1,len(eta[item])+1),eta[item],next(linecycler),label = item)
	legend = ax.legend(loc = 'upper right')

plt.xlabel('time')
plt.ylabel(r'$\eta$')
plt.savefig(filedirectory_to)
plt.clf()