#include <vector>
#include <iostream>
#include <fstream>
#include <Eigen/Sparse>
//#include <unsupported/Eigen/SparseExtra>


using namespace std;

//create sparse matrix with a +1 on diagonals starting at the ith column
//and create -1 on the diagonal starting at the jth column
//negative indices are interpreted as a row indication instead
//all indices start at 0

typedef Eigen::Triplet<int> T;
Eigen::SparseMatrix<int> createSparse(int i, int j, int x, int y, int z)
{
	std::vector<T> triples;
	int N;
	N = x*y*z;
	for (int a = 0; a<N; a += 1) {
		//add the +1 element
		if (a+i>=0 && a+i<N) {
			triples.push_back(T(a,a+i,1));
		}
		//add the -1 element
		if (a+j>=0 && a+j < N) {
			triples.push_back(T(a,a+j,-1));
		}
	}
	Eigen::SparseMatrix<int> mat(N,N);
	mat.setFromTriplets(triples.begin(),triples.end());

	return mat;


}

int main() {
	int x,y,z;
	int i = 0;
	int j = 9;
	x = y = z = 3;

	Eigen::SparseMatrix<int> dez,dhz;

	dez = createSparse(j,i,x,y,5);
	dhz = createSparse(0,-9,x,y,5);
	for (int k = 0; k< dhz.outerSize(); ++k) {
		for (Eigen::SparseMatrix<int>::InnerIterator it(dhz,k); it; ++it) {
			cout<<it.value()<<"\t";
			cout << it.row()<<"\t";
			cout << it.col() <<endl;
		}
	}
	

	return 0;
}

