#include <vector>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define E_SIZE 100
#define H_SIZE 99

#define BLOCK 1024


//this program will assume a 98x98x98 grid with 2 cells of zero padding for the E fields
//the padded zeros act as PEC boundaries
//The H fields will be 99x99x99 (offset by half cell, inside the PEC boundary)

//This version is a very naive version without matrix versions of the calculations
using namespace std;

__global__ void InitWall(double* ey, double init) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    ey[tid] = init;
}

//calculate ex_{i,j,k} for the next time step
//depends on: ex_{i,j,k} for the current time step, hz of adjacent cells, hy of adj. cells,
//the time step, epsilon, and the cell steps
//ended up using this as the general calculation for all E and H components
__device__ double Calc(double exn, double hzp, double hzn, double hyp, double hyn, double d1, double d2, double perm, double dt) {
    double term1, term2;
    double t1;
    t1 = hzp - hzn;
    term1 = t1/d1;
    term2 = (hyp - hyn)/d2;
    return dt*(term1-term2)/perm+exn;
}

__device__ int E2H(int index) {
    int i = index / (E_SIZE*E_SIZE);
    index -= (E_SIZE*E_SIZE)*i;
    int j = index / E_SIZE;
    int k = index - E_SIZE*j;
    return i*H_SIZE*H_SIZE+j*H_SIZE+k;
}

__device__ int H2E(int index) {
    int i = index / (H_SIZE*H_SIZE);
    index -= (H_SIZE*H_SIZE)*i;
    int j = index / H_SIZE;
    int k = index - H_SIZE*j;
    return i*E_SIZE*E_SIZE+j*E_SIZE+k;
}

__global__ void Set_H_X(double* hx, double* ey, double* ez, double mu, int size, double dt) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  // don't do anything for any thread ids that are greater
  if (tid < size) {
    double old, t1p, t1n, t2p, t2n;
    int y_offset = 100;

    old = hx[tid];
    int edex = H2E(tid);
		t1p = ey[edex+1];
	  t1n = ey[edex];
	  t2p = ez[edex+y_offset];
		t2n = ez[edex];
    hx[tid] = Calc(old,t1p,t1n,t2p,t2n,1.0,1.0,mu,dt);
  }
}

__global__ void Set_H_Y(double* hy, double* ez, double* ex, double mu, int size, double dt) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  // don't do anything for any thread ids that are greater
  if (tid < size) {
    double old, t1p, t1n, t2p, t2n;
    int x_offset = 10000;
    old = hy[tid];
    int edex = H2E(tid);
		t1p = ez[edex+x_offset];
	  t1n = ez[edex];
	  t2p = ex[edex+1];
		t2n = ex[edex];
    hy[tid] = Calc(old,t1p,t1n,t2p,t2n,1.0,1.0,mu,dt);
  }
}

__global__ void Set_H_Z(double* hz, double* ex, double* ey, double mu, int size, double dt) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  // don't do anything for any thread ids that are greater
  if (tid < size) {
    double old, t1p, t1n, t2p, t2n;
    int x_offset = 10000;
    int y_offset = 100;

    old = hz[tid];
    int edex = H2E(tid);
		t1p = ex[edex+y_offset];
	  t1n = ex[edex];
	  t2p = ey[edex+x_offset];
		t2n = ey[edex];
    hz[tid] = Calc(old,t1p,t1n,t2p,t2n,1.0,1.0,mu,dt);
  }
}

// __global__ void Set_E_X(double* ex, double* hz, double* hy, double eps, int size, double dt) {
//   int tid = blockIdx.x * blockDim.x + threadIdx.x;
//   // respect the border
//   int n = 100;
//   int n2 = 10000;
//   int n3 = 1000000;
//   if (tid >= size ||((tid < n2) || (tid > (n3-n2))) || (((0 < tid%n2) && (tid%n2 < n)) || ((n2-n < tid%n2) && (tid%n2 < n2-1))) || ((tid%n == 0) || (tid%n == n-1))){}
//   // if ((tid < n2) || (tid > (n3-n2))){}
//   // else if (((0 < tid%n2) && (tid%n2 < n)) || ((n2-n < tid%n2) && (tid%n2 < n2-1))){}
//   // else if ((tid%n == 0) || (tid%n == n-1)){}
//   else {
//     // printf("%d\n", tid);
//     double old, t1p, t1n, t2p, t2n;
//
//     old = ex[tid];
// 		t1p = hz[tid];
// 	  t1n = hz[tid-H_SIZE];
// 	  t2p = hy[tid];
// 		t2n = hy[tid-1];
//     ex[tid] = Calc(old,t1p,t1n,t2p,t2n,1.0,1.0,eps,dt);
//   }
// }
//
// __global__ void Set_E_Y(double* ey, double* hx, double* hz, double eps, int size, double dt) {
//   int tid = blockIdx.x * blockDim.x + threadIdx.x;
//   // respect the border
//   // int n = E_SIZE;
//   // int n2 = E_SIZE*E_SIZE;
//   // int n3 = n2 * E_SIZE;
//
//   int n = 100;
//   int n2 = 10000;
//   int n3 = 1000000;
//   if (tid >= size ||((tid < n2) || (tid > (n3-n2))) || (((0 < tid%n2) && (tid%n2 < n)) || ((n2-n < tid%n2) && (tid%n2 < n2-1))) || ((tid%n == 0) || (tid%n == n-1))){}
//   //
//   // if ((tid < n2) || (tid > (n3-n2))){}
//   // else if (((0 < tid%n2) && (tid%n2 < n)) || ((n2-n < tid%n2) && (tid%n2 < n2-1))){}
//   // else if ((tid%n == 0) || (tid%n == n-1)){}
//   else {
//     double old, t1p, t1n, t2p, t2n;
//
//     old = ey[tid];
// 		t1p = hx[tid];
// 	  t1n = hx[tid-1];
// 	  t2p = hz[tid];
// 		t2n = hz[tid-H_SIZE*H_SIZE];
//     ey[tid] = Calc(old,t1p,t1n,t2p,t2n,1.0,1.0,eps,dt);
//   }
// }
//
// __global__ void Set_E_Z(double* ez, double* hy, double* hx, double eps, int size, double dt) {
//   int tid = blockIdx.x * blockDim.x + threadIdx.x;
//   // respect the border
//   // int n = E_SIZE;
//   // int n2 = E_SIZE*E_SIZE;
//   // int n3 = n2 * E_SIZE;
//
//   int n = 100;
//   int n2 = 10000;
//   int n3 = 1000000;
//   if (tid >= size ||((tid < n2) || (tid > (n3-n2))) || (((0 < tid%n2) && (tid%n2 < n)) || ((n2-n < tid%n2) && (tid%n2 < n2-1))) || ((tid%n == 0) || (tid%n == n-1))){}
//   //
//   // if ((tid < n2) || (tid > (n3-n2))){}
//   // else if (((0 < tid%n2) && (tid%n2 < n)) || ((n2-n < tid%n2) && (tid%n2 < n2-1))){}
//   // else if ((tid%n == 0) || (tid%n == n-1)){}
//   else {
//     double old, t1p, t1n, t2p, t2n;
//
//     old = ez[tid];
// 		t1p = hy[tid];
// 	  t1n = hy[tid-H_SIZE*H_SIZE];
// 	  t2p = hx[tid];
// 		t2n = hx[tid-H_SIZE];
//     ez[tid] = Calc(old,t1p,t1n,t2p,t2n,1.0,1.0,eps,dt);
//   }
// }

__global__ void Set_E_X(double* ex, double* hz, double* hy, double eps, int size, double dt, int* inner_indices) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  // respect the border

  if (tid < size) {
    // printf("%d\n", tid);
    double old, t1p, t1n, t2p, t2n;
    int index = inner_indices[tid];
    old = ex[index];
    int hdex = E2H(index);
		t1p = hz[hdex];
	  t1n = hz[hdex-99];
	  t2p = hy[hdex];
		t2n = hy[hdex-1];
    ex[index] = Calc(old,t1p,t1n,t2p,t2n,1.0,1.0,eps,dt);
  }
}

__global__ void Set_E_Y(double* ey, double* hx, double* hz, double eps, int size, double dt, int* inner_indices) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  // respect the border
  if (tid < size) {
    double old, t1p, t1n, t2p, t2n;
    int index = inner_indices[tid];
    int hdex = E2H(index);
    old = ey[index];
		t1p = hx[hdex];
	  t1n = hx[hdex-1];
	  t2p = hz[hdex];
		t2n = hz[hdex-99*99];
    ey[index] = Calc(old,t1p,t1n,t2p,t2n,1.0,1.0,eps,dt);
  }
}

__global__ void Set_E_Z(double* ez, double* hy, double* hx, double eps, int size, double dt, int* inner_indices) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  // respect the border
  if (tid < size) {
    double old, t1p, t1n, t2p, t2n;
    int index = inner_indices[tid];
    int hdex = E2H(index);
    old = ez[index];
		t1p = hy[hdex];
	  t1n = hy[hdex-99*99];
	  t2p = hx[hdex];
		t2n = hx[hdex-99];
    ez[index] = Calc(old,t1p,t1n,t2p,t2n,1.0,1.0,eps,dt);
  }
}


// Used for time keeping independent of the clock
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
// double ex[E_SIZE][E_SIZE][E_SIZE] = {0};
// double ey[E_SIZE][E_SIZE][E_SIZE] = {0};
// double ez[E_SIZE][E_SIZE][E_SIZE] = {0};
// double hx[H_SIZE][H_SIZE][H_SIZE] = {0};
// double hy[H_SIZE][H_SIZE][H_SIZE] = {0};
// double hz[H_SIZE][H_SIZE][H_SIZE] = {0};

// enum Field {e_x, e_y, e_z, h_x, h_y, h_z};


// This is source term
// the argument is time value
//
double source(double t) {
	double expnum;
	expnum = pow(t-5e-7,2.0);
	return exp(-1*expnum/1e-15);
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


//function to calculate the magnitude of a 3-vector
//used to write out results, not part of simulation
double magn(double x, double y, double z) {
  // cout << x << ' ' << y << ' ' << z << endl;
	double mag;
	mag = sqrt(x*x+y*y+z*z);
	return mag;
}
// frunction to write out the magnitude of the E-field to a file
int write_to(ofstream& f, double t, int ind, int stride, double* ex, double* ey, double* ez) {
	f << t;
	int i;
	for (i = 0; i < nx; i+=stride) {
    int index = ind*100*100+i*100+49; // middle index for 100
		f << "\t" << magn(ex[index],ey[index],ez[index]);
	}
	f << "\t" << ind << "\t" << i << "\t" << 49;
	f << endl;
	return 0;
}

// primary simulation chunk
int main() {
  double eps = 8.85e-12;
	double mu = 1.257e-6;

  double *ex, *ey, *ez, *hx, *hy, *hz;
  int *inner_indices;
  int e_size = E_SIZE*E_SIZE*E_SIZE;
  int h_size = H_SIZE*H_SIZE*H_SIZE;
  int i_size = 98*98*98;

  ex = (double *)malloc((e_size)*sizeof(double));
  ey = (double *)malloc((e_size)*sizeof(double));
  ez = (double *)malloc((e_size)*sizeof(double));
  hx = (double *)malloc((h_size)*sizeof(double));
  hy = (double *)malloc((h_size)*sizeof(double));
  hz = (double *)malloc((h_size)*sizeof(double));

  inner_indices = (int *)malloc((i_size)*sizeof(int));

  // initialize to zero
  for (int i = 0; i < e_size; i++) {
    ex[i] = 0.0;
    ey[i] = 0.0;
    ez[i] = 0.0;
  }

  for (int i = 0; i < h_size; i++) {
    hx[i] = 0.0;
    hy[i] = 0.0;
    hz[i] = 0.0;
  }

  // cuda variables
  double *d_ex, *d_ey, *d_ez, *d_hx, *d_hy, *d_hz;
  int *d_inner;

  cudaMalloc((void **)&d_ex, sizeof(double) * (e_size));
  cudaMalloc((void **)&d_ey, sizeof(double) * (e_size));
  cudaMalloc((void **)&d_ez, sizeof(double) * (e_size));
  cudaMalloc((void **)&d_hx, sizeof(double) * (h_size));
  cudaMalloc((void **)&d_hy, sizeof(double) * (h_size));
  cudaMalloc((void **)&d_hz, sizeof(double) * (h_size));
  cudaMalloc((void **)&d_inner, sizeof(int) * (i_size));

  cudaMemcpy(d_ex, ex, sizeof(double) * (e_size), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ey, ey, sizeof(double) * (e_size), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ez, ez, sizeof(double) * (e_size), cudaMemcpyHostToDevice);
  cudaMemcpy(d_hx, hx, sizeof(double) * (h_size), cudaMemcpyHostToDevice);
  cudaMemcpy(d_hy, hy, sizeof(double) * (h_size), cudaMemcpyHostToDevice);
  cudaMemcpy(d_hz, hz, sizeof(double) * (h_size), cudaMemcpyHostToDevice);

	nx = E_SIZE;
	ny = E_SIZE;
	nz = E_SIZE;
	dx = 1;
	dy = 1;
	dz = 1;
	dt = 1e-9;
	// cout << "middle element is: " << ex[49][49][49] << endl;
	//the courant condition for 1 meter is 1.9e-9
	//final time be 1e-6 (for 1000 time steps)
	double tf = 1e-6;
	double t = 0.0;
	int a = 0;

  ofstream outFiles[11];
  stringstream fname;
  for (int it = 0; it < 11; it++) {
    fname.str("");
    fname << "paraOut/output" << it << ".txt";
    outFiles[it].open(fname.str());
  };

  int outind;

  ofstream probef;
  probef.open("paraOut/test.txt");
  ofstream probef2;
  probef2.open("paraOut/test_h.txt");

	double difference, w_start, w_finish;

	w_start = get_wall_time();

  int numBlocksH = h_size/BLOCK+1;
  int numBlocksI = i_size/BLOCK+1;
  dim3 threadsPerBlock(BLOCK, 1); // Max one dimensional block


  int count = 0;


  for (int i = 1; i < H_SIZE; i++) {
    for (int j = 1; j < H_SIZE; j++) {
      for (int k = 1; k < H_SIZE; k++) {
        inner_indices[count] = i*E_SIZE*E_SIZE+j*E_SIZE+k;
        count++;
      }
    }
  }

  cudaMemcpy(d_inner, inner_indices, sizeof(int) * (i_size), cudaMemcpyHostToDevice);

	while (t<tf) {
		cout << "t = " <<t <<endl;
		// set the source value for the incoming plane wave at x boundary
		double ey_init = source(t);

    // InitWall<<<10,1000>>>(d_ey, ey_init);

    for (int g = 0; g < ny; g++) {
      for (int h = 0; h < nz; h++) {
        int index = g*100 + h;
        ey[index] = ey_init;
      }
    }

    // Every tenth time step, write out slices of e-field values to a set of files
    if (!(a%10)) {
      cudaMemcpy(ex, d_ex, sizeof(double) * e_size, cudaMemcpyDeviceToHost);
      cudaMemcpy(ez, d_ez, sizeof(double) * e_size, cudaMemcpyDeviceToHost);
      cudaMemcpy(hy, d_hy, sizeof(double) * h_size, cudaMemcpyDeviceToHost);
      for (int fn = 0; fn < 11; fn++) {
        outind = fn*10;
        if (outind > 99) {
          outind = 99;
        }
        write_to(outFiles[fn], t, outind, 10, ex, ey, ez);
      }
      probef << t;

      // write to a couple of debug probes placed in the center of the box
      for (int y = 45; y < 55; y+=1) {
        int ex_index = 49*100*100+49*100+y;
        int hy_index = y*99*99+49*99+49;
        probef << "\t" << ex[ex_index];
        probef2 << "\t" << hy[hy_index];
      };
      probef << endl;
      probef2 << endl;
    };

    // after wall, might be better to do wall in CUDA
    cudaMemcpy(d_ey, ey, sizeof(double) * (e_size), cudaMemcpyHostToDevice);

    Set_H_X<<<numBlocksH, threadsPerBlock>>>(d_hx, d_ey, d_ez, mu, h_size, dt);
    cudaDeviceSynchronize();
    Set_H_Y<<<numBlocksH, threadsPerBlock>>>(d_hy, d_ez, d_ex, mu, h_size, dt);
    cudaDeviceSynchronize();
    Set_H_Z<<<numBlocksH, threadsPerBlock>>>(d_hz, d_ex, d_ey, mu, h_size, dt);
    cudaDeviceSynchronize();


    Set_E_X<<<numBlocksI, threadsPerBlock>>>(d_ex, d_hz, d_hy, eps, i_size, dt, d_inner);
    cudaDeviceSynchronize();
    Set_E_Y<<<numBlocksI, threadsPerBlock>>>(d_ey, d_hx, d_hz, eps, i_size, dt, d_inner);
    cudaDeviceSynchronize();
    Set_E_Z<<<numBlocksI, threadsPerBlock>>>(d_ez, d_hy, d_hx, eps, i_size, dt, d_inner);
    cudaDeviceSynchronize();

    // bring back for wall
    cudaMemcpy(ey, d_ey, sizeof(double) * e_size, cudaMemcpyDeviceToHost);
		t += dt; // time step counter
		a += 1; // printing counter
	}
	w_finish = get_wall_time();
	difference = w_finish - w_start;

    cout << "Parallel: " << difference << " seconds\n";

	probef.flush();
	probef.close();
	probef2.flush();
	probef2.close();
	for (int it = 0; it < 11; it++) {
		outFiles[it].flush();
		outFiles[it].close();
	};

	return 0;
}
