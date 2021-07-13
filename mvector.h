#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2
#include <vector>
#include <math.h>

// Class that represents a mathematical vector
class MVector{
public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}
	MVector(std::initializer_list<double> l) : v(l) {}
	
	// access element (lvalue)
	double &operator[](int index){
		return v[index];
	}
	
	// access element (rvalue)
	double operator[](int index) const{
		return v[index];
	}
	
	// access size of vector
	int size() const { return v.size(); } // number of elements
	
	// Evaluate the infinity-norm of the vector (Maximal value after absolute value taken)
	double LInfNorm() const{
		double maxAbs = 0;
		int VectorSize = size();
		for(int i=0; i<VectorSize; i++){
			maxAbs = std::max(std::abs(v[i]), maxAbs);
		}
		return maxAbs;
	}
	
	// Evaluate the L2-norm of the vector (Squareroot of sum of squares)
	double LTwoNorm() const{
		double SumOfSquares = 0;
		int VectorSize = size();
		for(int i=0; i<VectorSize; i++){
			SumOfSquares += pow(v[i], 2);
		}
		return sqrt(SumOfSquares);
	}
	
	
private:
	std::vector<double> v;
};

double dot(const MVector& Vect1, const MVector& Vect2){
	double dotProduct = 0.0;
	if(Vect1.size() == Vect2.size()){
		for(int i = 0; i < Vect1.size(); i++){
			dotProduct += (Vect1[i]*Vect2[i]);
		}
		return dotProduct;
	}
	else{
		std::cout << "ERROR: Vector lengths do not match \n";
		exit(1);
	}
}

// Overload the << operator to output MVectors to screen or file
std::ostream& operator<<(std::ostream& os, const MVector& v)
{
	int VectorSize = v.size();
	os << "(";
	for (int i=0; i < VectorSize-1; i++){
		os << v[i] << ", ";
	}
	os << v[VectorSize-1] << ")" << std::endl;
	return os;
}

//Overloading addition of 2 vectors
inline MVector operator+(MVector Vect1, MVector Vect2){
	MVector returnVect(Vect1.size());
	if(Vect1.size() == Vect2.size()){
		for(int i = 0; i < Vect1.size(); i++){
			returnVect[i] = Vect1[i] + Vect2[i];
		}
	}
	else{
		std::cout << "ERROR: Vector lengths do not match \n";
		exit(1);
	}
	return returnVect;
}

//Overloading subtraction of 2 vectors
inline MVector operator-(MVector Vect1, MVector Vect2){
	MVector returnVect(Vect1.size());
	if(Vect1.size() == Vect2.size()){
		for(int i = 0; i < Vect1.size(); i++){
			returnVect[i] = Vect1[i] - Vect2[i];
		}
	}
	else{
		std::cout << "ERROR: Vector lengths do not match \n";
		exit(1);
	}
	return returnVect;
}

// Overloading operator for 'scalar * vector'
inline MVector operator*(const double& lhs, const MVector& rhs){
	MVector temp(rhs);
	for (int i=0; i < temp.size(); i++) temp[i]*=lhs;
	return temp;
}

// Overloading operator for 'vector / scalar'
inline MVector operator/(const MVector& lhs, const double& rhs){
	MVector temp(lhs);
	for (int i=0; i < temp.size(); i++) temp[i]/=rhs;
	return temp;
}

//Overloading vector inner/dot product (vector*vector)
inline double operator*(MVector Vect1, MVector Vect2){
	double returnScalar = 0;
	if(Vect1.size() == Vect2.size()){
		for(int i = 0; i < Vect1.size(); i++){
			returnScalar += Vect1[i] * Vect2[i];
		}
	}
	else{
		std::cout << "ERROR: Vector lengths do not match \n";
		exit(1);
	}
	return returnScalar;
}







#endif
