#ifndef MMATRIX_H // the 'include guard'
#define MMATRIX_H
#include <vector>
#include <math.h>

// Class that represents a mathematical matrix
// Remember class uses C++ convention of zero-based indexing (i.e. (0,0) for top left entry rather than (1,1) in mathematical convention)
class MMatrix
{
public:
	// constructors
	MMatrix() : nRows(0), nCols(0) {}
	MMatrix(int n, int m, double x = 0) : nRows(n), nCols(m), A(n * m, x){}
	
	// access element, indexed by (row, column) [rvalue]
	double operator()(int i, int j) const
	{
		return A[j + i * nCols];
	}
	
	// access element, indexed by (row, column) [lvalue]
	double &operator()(int i, int j)
	{
		return A[j + i * nCols];
	}
	
	// size of matrix
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }
	
private:
	unsigned int nRows, nCols;
	std::vector<double> A;
};


// Overload the << operator to output MMatrix to screen or file
std::ostream& operator<<(std::ostream& os, const MMatrix& m){
	int rowsize = m.Rows();
	int columnsize = m.Cols();
	
	for(int i = 0; i < rowsize; i++){
		for(int j = 0; j < columnsize; j++){
			std::cout.width(5);
			os << m(i,j);
		}
		os << std::endl;
	}
	return os;
}


//Overloading Matrix/Vector Product
inline MVector operator*(const MMatrix& A, const MVector& v){
	MVector returnVect(A.Rows());
	double temp = 0.0;
	
	for(int i = 0; i < A.Rows(); i++){
		temp = 0.0;
		for(int j = 0; j < v.size(); j++){
			temp  = temp + A(i,j) * v[j];
		}
		returnVect[i] = temp;
	}
	return returnVect;
}

//Overloading Vector/Matrix Product
inline MVector operator*(const MVector& v, const MMatrix& A){
	MVector returnVect(A.Cols());
	double temp = 0.0;
	
	for(int i = 0; i < A.Cols(); i++){
		temp = 0.0;
		for(int j = 0; j < v.size(); j++){
			temp  = temp + v[j] * A(j,i);
		}
		returnVect[i] = temp;
	}
	return returnVect;
}

//Overloading Matrix/Matrix Product
inline MMatrix operator*(const MMatrix& A, const MMatrix& B){
	MMatrix returnMatrix(A.Rows(),B.Cols());
	double temp = 0.0;
	
	for(int i = 0; i < A.Rows(); i++){
		for(int j = 0; j < B.Cols(); j++){
			temp = 0.0;
			for(int k = 0; k < A.Cols(); k++){
				temp = temp + A(i,k) * B(k,j);
			}
			returnMatrix(i,j) = temp;
		}
	}
	return returnMatrix;
}


#endif
