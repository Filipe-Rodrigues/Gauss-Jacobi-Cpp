#ifndef GAUSS_SEIDEL_HPP
#define GAUSS_SEIDEL_HPP 1

#define UNFORMATTED 0
#define COMPACT_MATRIX 1
#define EXTENDED_MATRIX 2
#define EQUATION 3

#include <iostream>
#include <string>

using namespace std;

wchar_t BEGIN_BRACKET_1 = L"\u23a1"[0];
wchar_t BEGIN_BRACKET_2 = L"\u23a2"[0];
wchar_t BEGIN_BRACKET_3 = L"\u23a3"[0];
wchar_t END_BRACKET_1 = L"\u23a4"[0];
wchar_t END_BRACKET_2 = L"\u23a5"[0];
wchar_t END_BRACKET_3 = L"\u23a6"[0];
wchar_t FUNCTION_1 = L"\u23a7"[0];
wchar_t FUNCTION_2 = L"\u23a8"[0];
wchar_t FUNCTION_3 = L"\u23a9"[0];
wchar_t FUNCTION_A = L"\u23aa"[0];

wchar_t A0 = L"\u2080"[0];
wchar_t MULT = L"\u00d7"[0];
wchar_t EQUALS = L"\u003d"[0];
wchar_t BLANK = L"\u0020"[0];

class LinearSystem {
	private:
		int variableCount;
		double** A;
		double* x;
		double* b;
		bool solved;
		void loadInstance(string fileName);
		void initialize();
		void printVerticalArray(wostream& output), double* array, int size;
		void printUnformatted(wostream& output);
		void printCompactMatrix(wostream& output);
		void printExtendedMatrix(wostream& output);
		void printEquationMatrix(wostream& output);
	public:
		LinearSystem(string fileName);
		~LinearSystem();
		void setXValues(double* x) { this -> x = x; solved = true; }
		double* getEquation(int line);
		double* getXValues() { return x; }
		void print(wostream& output, int mode);
};

class Gauss_Seidel {
	private:
		LinearSystem* system;
		double tolerance;
		void computeAllRoots();
		double computeError(double* xPrev);
	public:
		Gauss_Seidel(string fileName, double tolerance);
		void computeRootsSequential();
		void print(ostream& output);
};

#endif