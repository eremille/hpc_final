//this program will assume a 98x98x98 grid with 2 cells of zero padding for the E fields
//the padded zeros act as PEC boundaries
//The H fields will be 99x99x99 (offset by half cell, inside the PEC boundary)

#include <vector>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define EX_SIZE 100
#define EY_SIZE 100
#define EZ_SIZE 100
#define HX_SIZE 99
#define HY_SIZE 99
#define HZ_SIZE 99

#define DX 1.0
#define DY 1.0
#define DZ 1.0

#define BLOCK 1024


using namespace std;

// __global__ functions are called by the host and invoke a kernel (must be void)
// __device__ functions are called by the device and are local to the gpu (can have a return value)

// int tid = blockIdx.x * blockDim.x + threadIdx.x; is the most common way to keep track of thread id
// each block can have at most 1024 threads where multiple of 32 threads are allocated per block are ideal
// blockIdx.x depends on how many numBlocks are passed
// blockDim.x refers to the size of each block that was assigned
// threadIdx.x can range from 0 to threadsPerBlock that was assigned, can have multiple dimensions

__global__ void InitWall(double* ey, double init, int size) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (tid < size) {
      ey[tid] = init;
    }
}

//calculate e/h{i,j,k} for the next time step
//depends on: e/h{i,j,k} for the current time step, hz of adjacent cells, hy of adj. cells,
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

// conversion for the 1D array
__device__ int E2H(int index) {
    int i = index / (EY_SIZE*EZ_SIZE);
    index -= (EY_SIZE*EZ_SIZE)*i;
    int j = index / EZ_SIZE;
    int k = index - EZ_SIZE*j;
    return i*HY_SIZE*HZ_SIZE+j*HZ_SIZE+k;
}

// conversion for the 1D array 
__device__ int H2E(int index) {
    int i = index / (HY_SIZE*HZ_SIZE);
    index -= (HY_SIZE*HZ_SIZE)*i;
    int j = index / HZ_SIZE;
    int k = index - HZ_SIZE*j;
    return i*EY_SIZE*EZ_SIZE+j*EZ_SIZE+k;
}

__global__ void Set_H_X(double* hx, double* ey, double* ez, double mu, int size, double dt) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  // don't do anything for any thread ids that are greater
  if (tid < size) {
    double old, t1p, t1n, t2p, t2n;
    int edex = H2E(tid);
    int y_offset = EZ_SIZE;

    old = hx[tid];
		t1p = ey[edex+1];
	  t1n = ey[edex];
	  t2p = ez[edex+y_offset];
		t2n = ez[edex];
    hx[tid] = Calc(old,t1p,t1n,t2p,t2n,DZ,DY,mu,dt);
  }
}

__global__ void Set_H_Y(double* hy, double* ez, double* ex, double mu, int size, double dt) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  // don't do anything for any thread ids that are greater
  if (tid < size) {
    double old, t1p, t1n, t2p, t2n;
    int x_offset = EY_SIZE*EZ_SIZE;

    old = hy[tid];
    int edex = H2E(tid);
		t1p = ez[edex+x_offset];
	  t1n = ez[edex];
	  t2p = ex[edex+1];
		t2n = ex[edex];
    hy[tid] = Calc(old,t1p,t1n,t2p,t2n,DX,DZ,mu,dt);
  }
}

__global__ void Set_H_Z(double* hz, double* ex, double* ey, double mu, int size, double dt) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  // don't do anything for any thread ids that are greater
  if (tid < size) {
    double old, t1p, t1n, t2p, t2n;
    int x_offset = EY_SIZE*EZ_SIZE;
    int y_offset = EZ_SIZE;

    old = hz[tid];
    int edex = H2E(tid);
		t1p = ex[edex+y_offset];
	  t1n = ex[edex];
	  t2p = ey[edex+x_offset];
		t2n = ey[edex];
    hz[tid] = Calc(old,t1p,t1n,t2p,t2n,DY,DX,mu,dt);
  }
}

__global__ void Set_E_X(double* ex, double* hz, double* hy, double eps, int size, double dt, int* inner_indices) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  // respect the border
  if (tid < size) {
    double old, t1p, t1n, t2p, t2n;
    int index = inner_indices[tid];
    int hdex = E2H(index);
    int y_offset = HZ_SIZE;

    old = ex[index];
		t1p = hz[hdex];
	  t1n = hz[hdex-y_offset];
	  t2p = hy[hdex];
		t2n = hy[hdex-1];
    ex[index] = Calc(old,t1p,t1n,t2p,t2n,DY,DX,eps,dt);
  }
}

__global__ void Set_E_Y(double* ey, double* hx, double* hz, double eps, int size, double dt, int* inner_indices) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  // respect the border
  if (tid < size) {
    double old, t1p, t1n, t2p, t2n;
    int index = inner_indices[tid];
    int hdex = E2H(index);
    int x_offset = HY_SIZE * HZ_SIZE;

    old = ey[index];
		t1p = hx[hdex];
	  t1n = hx[hdex-1];
	  t2p = hz[hdex];
		t2n = hz[hdex-x_offset];
    ey[index] = Calc(old,t1p,t1n,t2p,t2n,DZ,DX,eps,dt);
  }
}

__global__ void Set_E_Z(double* ez, double* hy, double* hx, double eps, int size, double dt, int* inner_indices) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  // respect the border
  if (tid < size) {
    double old, t1p, t1n, t2p, t2n;
    int index = inner_indices[tid];
    int hdex = E2H(index);
    int y_offset = HZ_SIZE;
    int x_offset = HY_SIZE*HZ_SIZE;

    old = ez[index];
		t1p = hy[hdex];
	  t1n = hy[hdex-x_offset];
	  t2p = hx[hdex];
		t2n = hx[hdex-y_offset];
    ez[index] = Calc(old,t1p,t1n,t2p,t2n,DX,DY,eps,dt);
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

// This is source term
// the argument is time value
//
double source(double t) {
	double expnum = pow(t-5e-7,2.0);
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
	for (i = 0; i < EX_SIZE; i+=stride) {
    int index = ind*EY_SIZE*EZ_SIZE+i*EZ_SIZE+EZ_SIZE/2-1; // middle index for 100
		f << "\t" << magn(ex[index],ey[index],ez[index]);
	}
	f << "\t" << ind << "\t" << i << "\t" << EZ_SIZE/2-1;
	f << endl;
	return 0;
}

// primary simulation chunk
int main() {
  double eps = 8.85e-12;
  double mu = 1.257e-6;
  double dt = 1e-9;
  double tf = 1e-6;
	double t = 0.0;
	int out_index = 0;

  double *ex, *ey, *ez, *hx, *hy, *hz;
  int *inner_indices;
  int e_size = EX_SIZE*EY_SIZE*EZ_SIZE;
  int h_size = HX_SIZE*HY_SIZE*HZ_SIZE;
  int i_size = (EX_SIZE-2)*(EY_SIZE-2)*(EZ_SIZE-2);
  int s_size = EY_SIZE*EZ_SIZE;

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

  // allocate memory for the cuda variables
  cudaMalloc((void **)&d_ex, sizeof(double) * (e_size));
  cudaMalloc((void **)&d_ey, sizeof(double) * (e_size));
  cudaMalloc((void **)&d_ez, sizeof(double) * (e_size));
  cudaMalloc((void **)&d_hx, sizeof(double) * (h_size));
  cudaMalloc((void **)&d_hy, sizeof(double) * (h_size));
  cudaMalloc((void **)&d_hz, sizeof(double) * (h_size));
  cudaMalloc((void **)&d_inner, sizeof(int) * (i_size));

  // copy memory from host to device
  cudaMemcpy(d_ex, ex, sizeof(double) * (e_size), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ey, ey, sizeof(double) * (e_size), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ez, ez, sizeof(double) * (e_size), cudaMemcpyHostToDevice);
  cudaMemcpy(d_hx, hx, sizeof(double) * (h_size), cudaMemcpyHostToDevice);
  cudaMemcpy(d_hy, hy, sizeof(double) * (h_size), cudaMemcpyHostToDevice);
  cudaMemcpy(d_hz, hz, sizeof(double) * (h_size), cudaMemcpyHostToDevice);

	// cout << "middle element is: " << ex[49][49][49] << endl;
	//the courant condition for 1 meter is 1.9e-9
	//final time be 1e-6 (for 1000 time steps)

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

  int numBlocksH = h_size/BLOCK+1; // set the numblock size to at least one
  int numBlocksI = i_size/BLOCK+1;
  int numBlocksS = s_size/BLOCK+1;
  dim3 threadsPerBlock(BLOCK, 1); // Max one dimensional block

  int count = 0;

  // keep track of the inner indices to respect the boundaries of the E-field
  for (int i = 1; i < HX_SIZE; i++) {
    for (int j = 1; j < HY_SIZE; j++) {
      for (int k = 1; k < HZ_SIZE; k++) {
        inner_indices[count] = i*EY_SIZE*EZ_SIZE+j*EZ_SIZE+k;
        count++;
      }
    }
  }

  cudaMemcpy(d_inner, inner_indices, sizeof(int) * (i_size), cudaMemcpyHostToDevice);

	while (t<tf) {
		cout << "t = " <<t <<endl;
		// set the source value for the incoming plane wave at x boundary
		double ey_init = source(t);

    InitWall<<<numBlocksS,threadsPerBlock>>>(d_ey, ey_init, s_size);

    // Every tenth time step, write out slices of e-field values to a set of files
    if (!(out_index%10)) {
      cudaMemcpy(ex, d_ex, sizeof(double) * e_size, cudaMemcpyDeviceToHost);
      cudaMemcpy(ey, d_ey, sizeof(double) * e_size, cudaMemcpyDeviceToHost);
      cudaMemcpy(ez, d_ez, sizeof(double) * e_size, cudaMemcpyDeviceToHost);
      cudaMemcpy(hy, d_hy, sizeof(double) * h_size, cudaMemcpyDeviceToHost);
      for (int fn = 0; fn < 11; fn++) {
        outind = fn*10;
        if (outind > HY_SIZE) {
          outind = HY_SIZE;
        }
        write_to(outFiles[fn], t, outind, 10, ex, ey, ez);
      }
      probef << t;

      // write to a couple of debug probes placed in the center of the box
      for (int y = 45; y < 55; y+=1) {
        int ex_index = 49*EY_SIZE*EZ_SIZE+49*EZ_SIZE+y;
        int hy_index = y*HY_SIZE*HZ_SIZE+49*HZ_SIZE+49;
        probef << "\t" << ex[ex_index];
        probef2 << "\t" << hy[hy_index];
      };
      probef << endl;
      probef2 << endl;
    };

    Set_H_X<<<numBlocksH, threadsPerBlock>>>(d_hx, d_ey, d_ez, mu, h_size, dt);
    Set_H_Y<<<numBlocksH, threadsPerBlock>>>(d_hy, d_ez, d_ex, mu, h_size, dt);
    Set_H_Z<<<numBlocksH, threadsPerBlock>>>(d_hz, d_ex, d_ey, mu, h_size, dt);
    cudaDeviceSynchronize(); // waits for kernels to return before continuing on the CPU

    Set_E_X<<<numBlocksI, threadsPerBlock>>>(d_ex, d_hz, d_hy, eps, i_size, dt, d_inner);
    Set_E_Y<<<numBlocksI, threadsPerBlock>>>(d_ey, d_hx, d_hz, eps, i_size, dt, d_inner);
    Set_E_Z<<<numBlocksI, threadsPerBlock>>>(d_ez, d_hy, d_hx, eps, i_size, dt, d_inner);
    cudaDeviceSynchronize();


		t += dt; // time step counter
		out_index += 1; // printing counter
	}
	w_finish = get_wall_time();
	difference = w_finish - w_start;

  cout << "Parallel: " << difference << " seconds\n";

  // probe clean up 
	probef.flush();
	probef.close();
	probef2.flush();
	probef2.close();
	for (int it = 0; it < 11; it++) {
		outFiles[it].flush();
		outFiles[it].close();
  };
  
  // memory clean up
  free(ex);
  free(ey);
  free(ez);
  free(hx);
  free(hy);
  free(hz);
  free(inner_indices);

  cudaFree(d_ex);
  cudaFree(d_ey);
  cudaFree(d_ez);
  cudaFree(d_hx);
  cudaFree(d_hy);
  cudaFree(d_hz);
  cudaFree(d_inner);

	return 0;
}
