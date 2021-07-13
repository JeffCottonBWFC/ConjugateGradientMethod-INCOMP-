#include <iostream>
#include <fstream>
#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"

//Initialise the tridiagonal matrix for MMatrix
MMatrix TridiagonalMatrixInitialisation(MMatrix& ReturnMatrix){
	
	for(int i = 0; i < ReturnMatrix.Rows(); i++){
		for(int j = 0; j < ReturnMatrix.Cols(); j++){
			if(i==j){
				ReturnMatrix(i,j) = 2;
			}
			else if (std::abs(i-j) == 1){
				ReturnMatrix(i,j) = -1;
			}
			else{
				ReturnMatrix(i,j) = 0;
			}
		}
	}
	return ReturnMatrix;
}

//Initialise the b vector
MVector VectorBInitialisation(MVector& ReturnVector){
	
	for(int i = 0; i < ReturnVector.size(); i++){
		ReturnVector[i] = 1/pow((ReturnVector.size()+1),2);
	}
	
	return ReturnVector;
}


//Initialise the Altered matrix
MMatrix AlteredMatrixInitialisation(MMatrix& ReturnMatrix, double m){
	
	for(int i = 0; i < ReturnMatrix.Rows(); i++){
		for(int j = 0; j < ReturnMatrix.Cols(); j++){
			if(i==j){
				ReturnMatrix(i,j) = 2*pow(i+1,2) + m;
			}
			else if (std::abs(i-j) == 1){
				ReturnMatrix(i,j) = -pow(i+1,2);
			}
			else{
				ReturnMatrix(i,j) = 0;
			}
		}
	}
	return ReturnMatrix;
}

//Initialise the b vector
MVector AlteredVectorBInitialisation(MVector& ReturnVector){
	
	for(int i = 0; i < ReturnVector.size(); i++){
		ReturnVector[i] = 2.5;
	}
	
	return ReturnVector;
}


//Implementation of the Conjugate Gradient Method
int CGM(MMatrix& A, MVector& b, MVector& x){
	int maxIterations = 1000;
	double tolerance = 1e-6;
	
	//Vector intialisation
	MVector rold = b - A*x;
	MVector p = rold;
	MVector rnew;
	
	//Variable declaration
	int TotalIterations = 0;
	double alpha = 0;
	double beta = 0;
	
	//CGM method
	for (int CurrentIter=0; CurrentIter<maxIterations; CurrentIter++){
		alpha = dot(rold,rold)/(p*A*p);
		x = x + alpha * p;
		rnew = rold - alpha * (A * p);
		
		// check if solution is accurate enough
		if (rnew.LTwoNorm() < tolerance) break;
		
		beta = dot(rnew,rnew)/dot(rold,rold);
		p = rnew + beta * p;
		rold = rnew;
		TotalIterations = CurrentIter+1;
	}
	
	//Test whether system converges onto a solution
	if(TotalIterations == maxIterations){
		std::cout << "Does not converge" << "\n";
		return 0;
	}
	
	return TotalIterations;
}


int main(){
	
	MBandedMatrix B(5,5,1,2);
	B(0,0) = 1;
	B(0,1) = 6;
	B(0,2) = 10;
	B(1,0) = 13;
	B(1,1) = 2;
	B(1,2) = 0;
	B(1,3) = 11;
	B(2,0) = 14;
	B(2,1) = 3;
	B(2,2) = 8;
	B(2,3) = 12;
	B(3,2) = 0;
	B(3,3) = 4;
	B(3,4) = 9;
	B(4,3) = 16;
	B(4,4) = 5;
	std::cout << B;
	
//
//	MBandedMatrix BStar(5,4,1,2);
//
//	for(int i = 0; i < BStar.Rows(); i++){
//		for(int j = 0; j < BStar.Cols(); j++){
//			BStar(i,j+BStar.LBands()-i) = B(i,j);
//		}
//	}
//
//	std::cout << BStar << "\n";
	

}










//OLD CODE


//	//Code which calculates the number of iterations to converge with the altered matrix A and vector b. The variable m can be altered.
//	//Test values
//	//
//	//m = 10.0
//	//n = 10 -> 38, n = 25 -> 113, n = 100 -> 632
//	//
//	//m = 20.0
//	//n = 10 -> 28, n = 25 -> 82, n = 100 -> 385
//	//
//	//m = 50.0
//	//n = 10 -> 18, n = 25 -> 52, n = 100 -> 241
//	//
//	//m = 5.5
//	//n = 10 -> 52, n = 25 -> 145 , n = 100 -> No Convergence
//
//
//	int n = 100;
//	double m = 5.5;
//
//	int TotalIterations = 0;
//
//	MMatrix A(n,n);
//	AlteredMatrixInitialisation(A,m);
//	MVector b(n);
//	AlteredVectorBInitialisation(b);
//
//	MVector x0(n);
//
//	TotalIterations = CGM(A,b,x0);
//
//	std::cout << TotalIterations << "\n";



//Initialise discretisation of Poisson's Equation
//	int n = 25;
//	int TotalIterations = 0;
//
//	MMatrix A(n,n);
//	TridiagonalMatrixInitialisation(A);
//	MVector b(n);
//	VectorBInitialisation(b);
//
//	MVector x0(n);
//
//	TotalIterations = CGM(A,b,x0);
//
//	std::cout << TotalIterations << "\n";



//SIZE OF MATRIX VS NUMBER OF ITERATIONS TO CONVERGE
//	std::ofstream FileStream;
//	FileStream.open("5to100IterationConvergenceVals.txt");
//
//	for(int n=5; n<101; n+=5){
//		int MatrixSize = n;
//		MMatrix A(MatrixSize,MatrixSize);
//		MVector b(MatrixSize), x0(MatrixSize), r0(MatrixSize);
//
//		int convergenceIterations = 0;
//
//		TridiagonalMatrixInitialisation(A);
//		VectorBInitialisation(b);
//
//		convergenceIterations = CGM(A,b,x0);
//
//		FileStream << n << " " << convergenceIterations << "\n";
//
//		}
//	FileStream.close();



//	//TEST FILESTREAM
//	std::ofstream FileStream;
//	FileStream.open("100by100XVals.txt");
//	for(int i = 0; i < n; i++){
//		FileStream << x0[i] << "\n";
//	}
//	FileStream.close();








