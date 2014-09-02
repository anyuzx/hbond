import sys
from math import *
import hbond
import numpy as np
from time import strftime
import os.path
from sys import stdout

if len(sys.argv) != 3:
    print 'Please check the arguments passed to the script.'
    sys.exit()
else:
    pass

filedirectory = sys.argv[1]
filedirectory_to = sys.argv[2]

trajfilenam = 'HISTORY'
FIELDfilenam = 'FIELD'
inffilenam = 'README'

ni = 0
count = 1
countMin = 1
countMax = 100000

if os.path.isfile(filedirectory+inffilenam):
    with open(filedirectory+inffilenam,'r') as f:
        contr_header = f.readline()
else:
    print "EAME file doesn't exist. Please check the directory or create a README file"
    sys.exit()
    
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


with open(filedirectory + trajfilenam,'r') as f:
    for i, line in enumerate(f):
        if i == 1:
            keytrj = int(line.split()[0])
            imcon = int(line.split()[1])
            if imcon == 0:
                print "ERROR: No periodic boundary condition is applied. Program is terminated"
                sys.exit()
        elif i > 1:
            break

boxvector = []
hb_avrg = {}
count_out = 0
with open(filedirectory + trajfilenam,'r') as f:
    for i, line in enumerate(f):
        if count <= countMax and count >= countMin:
            if (count-countMin)/500 == count_out + 1:
                count_out += 1
                stdout.write('\r****** Number of configurations used:%d' % int(count_out*500))
                stdout.flush()
            if line.split()[0] == "timestep":
                boxvector = []
                natomsc = int(line.split()[2])    # number of atoms in the CONFIGURATION
                ni = i
                nii = 0
                count += 1
                position = {}
                molecule_index = 1
                natmm = natmm_list[str(molecule_index)]
                for j in range(ntypemolecule):
                    position[str(j+1)] = []
                nmolecule = MOLECULE_list[str(molecule_index)]['number']
                continue
            if ni == 0:
                continue
            elif ni > 0 and (i == ni+1 or i == ni+2):
                boxvector.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
            elif ni > 0 and i == ni+3:
                boxvector.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
                LX = sqrt(boxvector[0][0]**2+boxvector[0][1]**2+boxvector[0][2]**2)
                LY = sqrt(boxvector[1][0]**2+boxvector[1][1]**2+boxvector[1][2]**2)
                LZ = sqrt(boxvector[2][0]**2+boxvector[2][1]**2+boxvector[2][2]**2)
                dim_array = np.array([LX,LY,LZ])
            elif ni > 0 and i >= ni + 4 and (i-ni-4) % (keytrj+2) == 1 and i < ni + nii + 3 + nmolecule*natmm*(keytrj+2):    # read the coordinate line
                position[str(molecule_index)].append(float(line.split()[0]))
                position[str(molecule_index)].append(float(line.split()[1]))
                position[str(molecule_index)].append(float(line.split()[2]))
            elif ni > 0 and i == ni + nii + 3 + nmolecule*natmm*(keytrj+2):                                                  # reach the end of one trajectory timestep
                if (i-ni-4) % (keytrj+2) == 1:
                    position[str(molecule_index)].append(float(line.split()[0]))
                    position[str(molecule_index)].append(float(line.split()[1]))
                    position[str(molecule_index)].append(float(line.split()[2]))
                molecule_index += 1
                if molecule_index <= ntypemolecule:
                    nmolecule = MOLECULE_list[str(molecule_index)]['number']
                    natmm = natmm_list[str(molecule_index)]
                    nii = i - (ni+3)
                elif molecule_index == ntypemolecule+1:
                    for item in position:
                        position[item] = np.array(position[item])
                    for j in range(ntypemolecule):
                        for k in range(j,ntypemolecule):
                            if j == k:
                                hb = hbond.hbondone(position[str(j+1)],dim_array,order[str(j+1)],natmm_list[str(j+1)])
                                if str(j)+'-'+str(k) not in hb_avrg:
                                    hb_avrg[str(j)+'-'+str(k)] = []
                                #fact = (MOLECULE_list[str(j+1)]['number']**2)/total_number
                                hb_avrg[str(j)+'-'+str(k)].append(hb)
                            else:
                                hb = hbond.hbondtwo(position[str(j+1)],position[str(k+1)],dim_array,order[str(j+1)],order[str(k+1)],natmm_list[str(j+1)],natmm_list[str(k+1)])
                                if str(j)+'-'+str(k) not in hb_avrg:
                                    hb_avrg[str(j)+'-'+str(k)] = []
                                #fact = (MOLECULE_list[str(j+1)]['number']*MOLECULE_list[str(k+1)]['number'])/total_number
                                hb_avrg[str(j)+'-'+str(k)].append(hb)
        elif count < countMin:
            if line.split()[0] == 'timestep':
                count += 1
        else:
            break

#for pair in iter(hb_avrg):
#    print hb_avrg[pair]
#    hb_avrg[pair] = np.array(hb_avrg[pair])
#    hb_avrg[pair] = np.cumsum(hb_avrg[pair])/np.arange(1,len(hb_avrg[pair])+1,dtype=float)

with open(filedirectory_to,'w') as f:
    f.write(contr_header)
    f.write(strftime("%Y-%m-%d %H:%H:%S")+'\n')
    for pair in iter(hb_avrg):
        f.write('pair:' + pair + '\n')
        for item in hb_avrg[pair]:
            f.write(str(item)+'\n')
        f.write('\n')

print '******  '+contr_header+'\n******  Hydrogen Bond analysis'
print '******  '+'Read the Trajectory File:'+filedirectory+trajfilenam
print '******  '+'Write the result to the file:'+filedirectory_to
print '\n******  Some analysis information:'
print '******  '+'# of configurations'.ljust(25)+'# of atoms'.ljust(15)
print '******  '+str(count-1).ljust(25)+str(natomsc).ljust(15)
print '******  '+'The code is successfully completed'