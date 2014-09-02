#include <math.h>

/*
c_hbond.cpp

c++ function that calculate the average hydrogen bonds between molecules

the standard criteria for hydrogen bond is that the distance between two oxgens is less then
3.5 Angstrom and the angle of HOO is less than 30 degree

Arguments:
r_array: flatten array of coordinates of all atoms need to be calculated

** I want this code to give the average number of hydrogen bonds per molecules between 
** any given kinds of molecules
** this can be between the same kind of molecules or two different kinds of molecules
** for instance: average hydrogen bonds per molecules between GT and GT molecules
** or average hydrogen bonds per molecules between GT and GTRC molecules

** the input arguments should be two lists of all molecules coordinates I want to calculate
** if it's between same molecules, then the input is the list of molecules coordinates
*/

/*
n: number of molecules
*/

void c_hbondone(double* r_array, double* dim_array, int order, int n, int natmm, int& hbond){
	int i,j,k;
	double dx, dy, dz, rOO, cos_angle;
	double rOH_array [3], rOO_array [3];

	for(i=0;i<n-1;i++){
		for(j=i+1;j<n;j++){
			dx = r_array[j*natmm*3 + order*3] - r_array[i*natmm*3 + order*3];
			dy = r_array[j*natmm*3 + order*3 + 1] - r_array[i*natmm*3 + order*3 + 1];
			dz = r_array[j*natmm*3 + order*3 + 2] - r_array[i*natmm*3 + order*3 + 2];
			if (fabs(dx) >= dim_array[0]/2) {dx = -dx*(dim_array[0]/fabs(dx)-1);}
			if (fabs(dy) >= dim_array[1]/2) {dy = -dy*(dim_array[1]/fabs(dy)-1);}
			if (fabs(dz) >= dim_array[2]/2) {dz = -dz*(dim_array[2]/fabs(dz)-1);}
			rOO_array[0] = dx;
			rOO_array[1] = dy;
			rOO_array[2] = dz;
			rOO = sqrt(pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0)); // calculate the distance between two oxygen atoms
            
			// calculate the angle between HO bond and OO bond 
			if (rOO < 3.5)
			{
				for(k=0;k<natmm;k++){
					if (k == order)
						continue;
					else
					{
						dx = r_array[i*natmm*3 + k*3] - r_array[i*natmm*3 + order*3];
				        dy = r_array[i*natmm*3 + k*3 + 1] - r_array[i*natmm*3 + order*3 + 1];
				        dz = r_array[i*natmm*3 + k*3 + 2] - r_array[i*natmm*3 + order*3 + 2];
				        if (fabs(dx) >= dim_array[0]/2) {dx = -dx*(dim_array[0]/fabs(dx)-1);}
			            if (fabs(dy) >= dim_array[1]/2) {dy = -dy*(dim_array[1]/fabs(dy)-1);}
			            if (fabs(dz) >= dim_array[2]/2) {dz = -dz*(dim_array[2]/fabs(dz)-1);}
				        rOH_array[0] = dx;
				        rOH_array[1] = dy;
				        rOH_array[2] = dz;

				        cos_angle = (rOO_array[0]*rOH_array[0] + rOO_array[1]*rOH_array[1] + rOO_array[2]*rOH_array[2])/rOO;
				        if (cos_angle > sqrt(3.0)/2.0) {hbond++;break;}

				        dx = r_array[j*natmm*3 + k*3] - r_array[j*natmm*3 + order*3];
				        dy = r_array[j*natmm*3 + k*3 + 1] - r_array[j*natmm*3 + order*3 + 1];
				        dz = r_array[j*natmm*3 + k*3 + 2] - r_array[j*natmm*3 + order*3 + 2];
				        if (fabs(dx) >= dim_array[0]/2) {dx = -dx*(dim_array[0]/fabs(dx)-1);}
			            if (fabs(dy) >= dim_array[1]/2) {dy = -dy*(dim_array[1]/fabs(dy)-1);}
			            if (fabs(dz) >= dim_array[2]/2) {dz = -dz*(dim_array[2]/fabs(dz)-1);}
				        rOH_array[0] = dx;
				        rOH_array[1] = dy;
				        rOH_array[2] = dz;

				        cos_angle = (-rOO_array[0]*rOH_array[0] - rOO_array[1]*rOH_array[1] - rOO_array[2]*rOH_array[2])/rOO;
				        if (cos_angle > sqrt(3.0)/2.0) {hbond++;break;}
					}
			    }
			}
		}
	}
	return;
}

void c_hbondtwo(double* r1_array, double* r2_array, double* dim_array, int order1, int order2, int n1, int n2, int natmm1, int natmm2, int& hbond){
	int i,j,k,hbond_temp;
	double dx, dy, dz, rOO, cos_angle;
	double rOH_array [3], rOO_array [3];

	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			dx = r2_array[j*natmm2*3 + order2*3] - r1_array[i*natmm1*3 + order1*3];
			dy = r2_array[j*natmm2*3 + order2*3 + 1] - r1_array[i*natmm1*3 + order1*3 + 1];
			dz = r2_array[j*natmm2*3 + order2*3 + 2] - r1_array[i*natmm1*3 + order1*3 + 2];
			if (fabs(dx) >= dim_array[0]/2) {dx = -dx*(dim_array[0]/fabs(dx)-1);}
			if (fabs(dy) >= dim_array[1]/2) {dy = -dy*(dim_array[1]/fabs(dy)-1);}
			if (fabs(dz) >= dim_array[2]/2) {dz = -dz*(dim_array[2]/fabs(dz)-1);}
			rOO_array[0] = dx;
			rOO_array[1] = dy;
			rOO_array[2] = dz;
			rOO = sqrt(pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0)); // calculate the distance between two oxygen atoms
            
			// calculate the angle between HO bond and OO bond 
			if (rOO < 3.5)
			{
				hbond_temp = hbond;
				for(k=0;k<natmm1;k++){
					if (k == order1)
						continue;
					else
					{
						dx = r1_array[i*natmm1*3 + k*3] - r1_array[i*natmm1*3 + order1*3];
				        dy = r1_array[i*natmm1*3 + k*3 + 1] - r1_array[i*natmm1*3 + order1*3 + 1];
				        dz = r1_array[i*natmm1*3 + k*3 + 2] - r1_array[i*natmm1*3 + order1*3 + 2];
				        if (fabs(dx) >= dim_array[0]/2) {dx = -dx*(dim_array[0]/fabs(dx)-1);}
			            if (fabs(dy) >= dim_array[1]/2) {dy = -dy*(dim_array[1]/fabs(dy)-1);}
			            if (fabs(dz) >= dim_array[2]/2) {dz = -dz*(dim_array[2]/fabs(dz)-1);}
				        rOH_array[0] = dx;
				        rOH_array[1] = dy;
				        rOH_array[2] = dz; 

				        cos_angle = (rOO_array[0]*rOH_array[0] + rOO_array[1]*rOH_array[1] + rOO_array[2]*rOH_array[2])/rOO;
				        if (cos_angle > sqrt(3.0)/2.0) {hbond++;break;}
					}
			    }
			    if (hbond_temp == hbond) {
			    	for(k=0;k<natmm2;k++){
					    if (k == order2)
						    continue;
					    else
					    {
						    dx = r2_array[j*natmm2*3 + k*3] - r2_array[j*natmm2*3 + order2*3];
				            dy = r2_array[j*natmm2*3 + k*3 + 1] - r2_array[j*natmm2*3 + order2*3 + 1];
				            dz = r2_array[j*natmm2*3 + k*3 + 2] - r2_array[j*natmm2*3 + order2*3 + 2];
				            if (fabs(dx) >= dim_array[0]/2) {dx = -dx*(dim_array[0]/fabs(dx)-1);}
			                if (fabs(dy) >= dim_array[1]/2) {dy = -dy*(dim_array[1]/fabs(dy)-1);}
			                if (fabs(dz) >= dim_array[2]/2) {dz = -dz*(dim_array[2]/fabs(dz)-1);}
				            rOH_array[0] = dx;
				            rOH_array[1] = dy;
				            rOH_array[2] = dz;

				            cos_angle = (-rOO_array[0]*rOH_array[0] - rOO_array[1]*rOH_array[1] - rOO_array[2]*rOH_array[2])/rOO;
				            if (cos_angle > sqrt(3.0)/2.0) {hbond++;break;}
					    }
			        }
			    }
			}
		}
	}
	return;
}