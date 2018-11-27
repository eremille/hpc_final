#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
// #include <Eigen/Sparse>
#include <time.h>
#include <sys/time.h>
//#include <unsupported/Eigen/SparseExtra>

//this program will assume a 98x98x98 grid with 2 cells of zero padding for the E fields
//the padded zeros act as PEC boundaries
//The H fields will be 99x99x99 (offset by half cell, inside the PEC boundary)

//This version is a very naive version without matrix versions of the calculations
using namespace std;

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int nx, ny, nz;
double dx, dy, dz;
double dt;
double ex[100][100][100] = {0};
double ey[100][100][100] = {0};
double ez[100][100][100] = {0};
double hx[99][99][99] = {0};
double hy[99][99][99] = {0};
double hz[99][99][99] = {0};

//calculate ex_{i,j,k} for the next time step
//depends on: ex_{i,j,k} for the current time step, hz of adjacent cells, hy of adj. cells,
//the time step, epsilon, and the cell steps
//ended up using this as the general calculation for all E and H components
double calc_exijk(double exn, double hzp,double hzn,double hyp,double hyn, double d1, double d2,double perm) {
	double update;
	double term1, term2;
	double t1, t2;
	t1 = hzp - hzn;
	term1 = t1/d1;
	term2 = (hyp - hyn)/d2;
	update = dt*(term1-term2)/perm+exn;
	return update;
};

// This is source term
// the argument is time value
// 
double source(double t) {
	double val;
	double expnum;
	expnum = pow(t-5e-7,2.0);
	// val = exp(-1*expnum/1e-15);
	val = exp(-1*expnum/1e-15);
	// cout << val<<endl;
	return val;
};

//use existing procedures for all calculations
//the various calc_enijk procedures are the exact same math, we will just use one
//depends on: desired quantity, Ex array, Ey array, Ez array, Hx array, Hy array,
//Hz array (all of which are pointers), dx, dy, dz, i, j, k
//type values: 0, 1, 2, 3, 4, 5 = ex, ey, ez, hx, hy, hz
// double calc_int(int type, double ex[][ny][nz], double ey[][ny][nz], double ez[][ny][nz],
// 	double hx[][ny-1][nz-1], double hy[][ny-1][nz-1], double hz[][ny-1][nz-1],
// 	double dx, double dy, double dz,
// 	double dt, int i, int j, int k)
double calc_int(int type, int i, int j, int k) {
	double next_val;
	double old, t1p, t1n, t2p, t2n;
	double eps = 8.85e-12;
	double mu = 1.257e-6;
	// double eps = 1.0;
	// double mu = 1.0;
	switch(type) {
		//case for ex
		case 0: old = ex[i][j][k];
			t1p = hz[i][j][k];
			t1n = hz[i][j-1][k];
			t2p = hy[i][j][k];
			t2n = hy[i][j][k-1];
			next_val = calc_exijk(old,t1p,t1n,t2p,t2n,dy,dx,eps);
			break;
		//case for ey
		case 1: old = ey[i][j][k];
			t1p = hx[i][j][k];
			t1n = hx[i][j][k-1];
			t2p = hz[i][j][k];
			t2n = hz[i-1][j][k];
			next_val = calc_exijk(old,t1p,t1n,t2p,t2n,dz,dx,eps);
			break;
		//case for ez
		case 2: old = ez[i][j][k];
			t1p = hy[i][j][k];
			t1n = hy[i-1][j][k];
			t2p = hx[i][j][k];
			t2n = hx[i][j-1][k];
			next_val = calc_exijk(old,t1p,t1n,t2p,t2n,dx,dy,eps);
			break;
		//case for hx
		case 3: old = hx[i][j][k];
			t1p = ey[i][j][k+1];
			t1n = ey[i][j][k];
			t2p = ez[i][j+1][k];
			t2n = ez[i][j][k];
			next_val = calc_exijk(old,t1p,t1n,t2p,t2n,dz,dy,mu);
			break;
		//case for hy (needs fixing)
		case 4: old = hx[i][j][k];
			t1p = ez[i+1][j][k];
			t1n = ez[i][j][k];
			t2p = ex[i][j][k+1];
			t2n = ex[i][j][k];
			next_val = calc_exijk(old,t1p,t1n,t2p,t2n,dz,dy,mu);
			break;
		//case for hz
		case 5: old = hx[i][j][k];
			t1p = ex[i][j+1][k];
			t1n = ex[i][j][k];
			t2p = ey[i+1][j][k];
			t2n = ey[i][j][k];
			next_val = calc_exijk(old,t1p,t1n,t2p,t2n,dy,dx,mu);
			break;
	}
	// cout << "type: " << type << endl;
	// cout << "indices: " << i <<" "<< j <<" "<< k << endl;
	// cout << next_val << endl;
	// cout << "t1p: " << t1p << endl;
	// cout << "t1n: " << t1n << endl;
	// cout << "t2p: " << t2p << endl;
	// cout << "t2n: " << t2n << endl << endl;
	return next_val;
};
//
//function to calculate the magnitude of a 3-vector
//used to write out results, not part of simulation
double magn(double x, double y, double z) {
	double mag;
	mag = sqrt(x*x+y*y+z*z);
	return mag;
}
// frunction to write out the magnitude of the E-field to a file
int write_to(ofstream& f, double t, int ind, int stride) {
	f << t;
	for (int i = 0; i <= nx; i+=stride) {
		f << "\t" << magn(ex[ind][i][49],ey[ind][i][49],ez[ind][i][49]);
	}
	f << endl;
	return 0;
}

// primary simulation chunk
int main() {
	nx = 100;
	ny = 100;
	nz = 100;
	dx = 1;
	dy = 1;
	dz = 1;
	dt = 1e-10;
	// cout << "middle element is: " << ex[49][49][49] << endl;
	//the courant condition for 1 meter is 1.9e-9
	//final time be 1e-6 (for 1000 time steps)
	double tf = 1e-6;
	double t =0.0;
	double tmp;
	int a = 0;

	ofstream outFiles[11];
	stringstream fname;
	for (int it = 0; it < 11; it++) {
		fname.str("");
		fname << "output" << it << ".txt";
		outFiles[it].open(fname.str());
	};

	double test;
	ofstream probef;
	probef.open("test.txt");
	ofstream probef2;
	probef2.open("test_h.txt");
	//set the bounds of the loops
	int lx, ly, lz;
	// cout << nx <<endl;
	lx = nx - 1;
	ly = ny - 1;
	lz = nz - 1;

	int outind;

	double newex, newey, newez, newhx, newhy, newhz;

	double difference, w_start, w_finish;

	// cout << lx;
	w_start = get_wall_time();
	while (t<tf) {
		cout << "t = " <<t <<endl;
		// set the source value for the incoming plane wave at x boundary
		for (int g = 0; g < ny; g++) {
			for (int h = 0; h < nz; h++) {
				ey[0][g][h] = source(t);
				// hz[]
			};
		};
		// old source term, not used anymore but might be useful later
		// ex[49][49][49] = source(t);
		// hy[49][49][49] = source(t)/377;
		// ex[49][49][49] += 1;
		// hy[49][49][49] += 1;

		// Every tenth time step, write out slices of e-field values to a set of files
		if (!(a%10)) {
			for (int fn = 0; fn < 11; fn++) {
				outind = fn*10;
				if (outind>lx) {
					outind = lx;
				};
				write_to(outFiles[fn], t, outind, 10);
			}
			probef << t;

			// write to a couple of debug probes placed in the center of the box
			for (int y = 45; y < 55; y+=1) {
				probef << "\t" << ex[49][49][y];
				probef2 << "\t" << hy[y][49][49];
			};
			probef << endl;
			probef2 << endl;
		};
		// Calculate the E-fields first (this allows )
		// loop bounds for E-fields are 1 to Nx-1 to avoid overwritting the boundary conditions
		for (int i = 1; i < lx; i++) {
			// ex[i][0][0] = source(t);
			for (int j = 1; j < ly; j++) {
				for (int k = 1; k < lz; k++) {
					newey = calc_int(1,i,j,k);
					ey[i][j][k] = newey;
					newez = calc_int(2,i,j,k);
					ez[i][j][k] = newez;

					ex[i][j][k] = calc_int(0,i,j,k);

					// if (!(i == 49 && j == 49 && k == 49)) {
					// 	tmp = calc_int(0,i,j,k);
					// 	ex[i][j][k] = calc_int(0,i,j,k);
						// if (tmp != 0) {
						// 	cout << i << "\t" << j << "\t" << k << endl;
						// 	cout << t << endl;
						// 	cout <<"tmp value: " <<tmp <<endl;
						// 	cout <<"matrix value: " <<ex[i][j][k] <<endl;
						// 	return 0;
						// };
						// cout << tmp << "\t" << ex[i][j][k] << endl;

						// cout << ex[i][j][k] << endl;
					// };
					// cout << i << "\t" << j << "\t" << k << endl;
					// cout << newex << "\t" << newey << "\t" << newez << endl;
					// cout << newhx << "\t" << newhy << "\t" << newhz << endl;
					// else {
					// 	ex[i][j][k] = calc_int(0,i,j,k);
					// };
					// cout <<"calculated E: " << i << "\t" << j <<"\t" << k << endl;


					// cout <<"calculated H: " << i << "\t" << j <<"\t" << k << endl;
				}
			}
		}

		// Calculate H-fields
		for (int l = 0; l < lx; l += 1) {
			for (int m = 0; m <ly; m += 1) {
				for (int n = 0; n < lz; n += 1) {
					hx[l][m][n] = calc_int(3,l,m,n);
					hy[l][m][n] = calc_int(4,l,m,n);
					hz[l][m][n] = calc_int(5,l,m,n);
				}
			}
		}


		t += dt;
		a += 1;
	}
	w_finish = get_wall_time();
	difference = w_finish - w_start;

    cout << "Naive: " << difference << " seconds\n";


	probef.flush();
	probef.close();
	probef2.flush();
	probef2.close();
	for (int it = 0; it < 11; it++) {
		outFiles[it].flush();
		outFiles[it].close();
	};
	// test = source(4e-7);
	// cout << test << endl;



	return 0;
}
