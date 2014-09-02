import hbond
import numpy as np

p = np.array([0,0,0,1,0,0,-0.333313,-0.942816,0,2.73,0,0,3.49604,0.642788,0,3.0807,-0.936489,0],dtype = float)
order = 0
dim_array = np.array([10,10,10],dtype = float)
natmm = 3

print hbond.hbondone(p,dim_array,order,natmm)