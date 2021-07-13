#ifndef MBANDEDMATRIX_H // the 'include guard'
#define MBANDEDMATRIX_H

#include <vector>

class MBandedMatrix{
public:
	// constructors
	MBandedMatrix() : nRows(0), nCols(0) {}
	MBandedMatrix(int n, int m, int lband, int rband, double x = 0) : nRows(n), nCols(m), A(n * (lband + rband + 1), x), l(lband), r(rband) {}
	
	
	// access element [rvalue]
	double operator()(int i, int j) const
	{
		return A[j + i * nCols];
	}
	
	// access element [lvalue]
	double &operator()(int i, int j)
	{
		return A[j + i * nCols];
	}
	
	// size of matrix
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }
	int Bands() const { return r + l + 1; } // total number of bands
	int LBands() const { return l; } // number of left bands
	int RBands() const { return r; } // number of right bands
	
private:
	unsigned int nRows, nCols;
	std::vector<double> A;
	int l, r; // number of left/right diagonals
};


std::ostream& operator<<(std::ostream& output, const MBandedMatrix& banded){
	int r = banded.Rows(), c = banded.Cols();
	
	for (int i = 0; i < banded.Rows(); i++){
		
		// calculate position of lower and upper band
		int jmin = std::max(std::min(i-banded.LBands(), banded.Cols()),0);
		int jmax = std::min(i+banded.RBands()+1, banded.Cols());
		
		output << "( ";
		for (int j=0; j<jmin; j++)
			output << 0 << "\t ";
		for (int j=jmin; j<jmax; j++)
			output << banded(i,j) << "\t ";
		for (int j=jmax; j<c; j++)
			output << 0 << "\t ";
		output << ")\n";
	}
	return output;
}


#endif
