#ifndef GAUSS_SEIDEL_HPP
#define GAUSS_SEIDEL_HPP 1

class LinearSystem {
	private:
		int variableCount;
		double** A;
		double* x;
		double* b;
	public:
		LinearSystem();
		LinearSystem(int varCount);
		~LinearSystem();
		
}

#endif