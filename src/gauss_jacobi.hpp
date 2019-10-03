#ifndef GAUSS_JACOBI_HPP
#define GAUSS_JACOBI_HPP 1

#define RESULT_ONLY 0
#define UNFORMATTED 1
#define COMPACT_MATRIX 2
#define EXTENDED_MATRIX 3
#define EQUATION 4

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

class InputConfiguration {
	public:
		int format;
		int precision;
		int tolerance;
		string fileName;
		bool error;
		InputConfiguration() { format = 0; precision = 2; fileName = ""; tolerance = 5; error = true; }
		InputConfiguration(const InputConfiguration& conf2);
		bool checkConfiguration();
};

class LinearSystem {
	private:
		int variableCount;
		int precision;
		double** A;
		double* x;
		double* b;
		bool solved;
		void loadInstance(string fileName);
		void initialize();
		double* subs();
		void printResults(wostream& output);
		void printUnformatted(wostream& output);
		void printCompactMatrix(wostream& output);
		void printExtendedMatrix(wostream& output);
		void printEquationMatrix(wostream& output);
	public:
		LinearSystem(string fileName, int precision);
		~LinearSystem();
		void setPrecision(int newPrecision) { precision = newPrecision; }
		void setSolved(bool solved) { this -> solved = solved; }
		double* getEquation(int line);
		double* getXValues() { return x; }
		double* getB() { return b; }
		double* getResidualVector();
		int getVariableCount() { return variableCount; }
		bool ensureNonZeroDiagonal();
		void printVerticalArray(wostream& output, double* array, int size, wstring label);
		void print(wostream& output, int mode);
};

class Gauss_Jacobi {
	private:
		LinearSystem* system;
		double tolerance;
		InputConfiguration configuration;
		void computeAllRoots();
		void computeRoots(int tid);
		double computeRoot(int x_id);
		double computeError(double* xPrev);
		void printIntro(wostream& output);
		void printFullHeader(wostream& output);
		void printBasicHeader(wostream& output);
	public:
		Gauss_Jacobi(LinearSystem* system, double tolerance);
		Gauss_Jacobi(InputConfiguration config);
		~Gauss_Jacobi() { delete system; }
		void computeRootsSequential() { if (configuration.checkConfiguration()) computeAllRoots(); }
};

#endif