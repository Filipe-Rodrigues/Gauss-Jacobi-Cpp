#include "gauss_jacobi.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

//////////////////////////         AUXILIAR FUNCTIONS  ////////////////////////
///////////////////////////////////////////////////////////////////////////////

double* copy(double* original, int vars) {
	double* newArr = new double[vars];
	for (int i = 0; i < vars; i++) {
		newArr[i] = original[i];
	}
	return newArr;
}

void printError() {

}

int getColumnWidth(double* values, int size, int precision) {
	int width = 0;
	double numLength;
	for (int i = 0; i < size; i++) {
		numLength = floor(log10(abs(values[i]))) + 1;
		if (values[i] < 1) {
			numLength = 1;
		}
		if (values[i] < 0) {
			numLength++;
		}
		numLength += precision + 1;
		if (numLength > width) {
			width = numLength;
		}
	}
	return width + 1;
}

int getColumnWidth(double** values, int size, int precision) {
	int width = 0;
	double numLength;
	for (int i = 0; i < size; i++) {
		numLength = getColumnWidth(values[i], size, precision) - 1;
		if (numLength > width) {
			width = numLength;
		}
	}
	return width + 1;
}

wstring getSubsFromIndex(int index) {
	int size = (int) log10(index);
	if (index == 0) {
		size = 0;
	}
	wstring subs = L"";
	for (int i = (int) pow(10, size); i >= 1; i /= 10) {
		int c = index / i;
		wchar_t newChar = A0 + c;
		subs += newChar;

		index -= c * i;
	}
	return subs;
}

void printIntro(wostream& output) {
	wchar_t i = L"\u1d62"[0];
	wchar_t j = L"\u2c7c"[0];
	wchar_t sum = L"\u2211"[0];
	wchar_t uneq = L"\u2260"[0];
	output << "######## Hello! This is the GAUSS-JACOBI method for Linear Systems!\n";
	output << "###################################################################\n\n";
	output << "    This work was made by the team FHRR, hope you enjoy it as we do ;)\n\n";
	output << "    The Gauss-Jacobi method is very simple, it consists in isolating each variable from the system,\n";
	output << "then using each of the obtained equations to calculate new values for their respective variables.\n";
	output << "For each iteration it does, the value of the variables gets more and more close to the solution\n";
	output << "of the system, and the number of iterations is directly proportional to the adopted tolerance, so\n";
	output << "the greater the precision desired, more iterations will be needed.\n\n";
	output << "    There's two conditions to assure the convergence of the method: the main diagonal must be free\n";
	output << "from zero elements, and the matrix must be diagonally dominant, that is:\n\n";
	output << L"\t\t\tFor each line i, |a" << i << i << "| > " << sum << "|a";
	output << i << j << "|; i" << uneq << "j\n\n";

}

//////////////////////////            IMPLEMENTATIONS  ////////////////////////
///////////////////////////////////////////////////////////////////////////////

//////////////////////////#### InputConfiguration:

bool InputConfiguration::checkConfiguration() {
	return (format >= 0 and precision >= 0 and fileName != "");
}

//////////////////////////#### LinearSystem:

LinearSystem::LinearSystem(string fileName, int precision) {
	loadInstance(fileName);
	this -> precision = precision;
	solved = false;
}

LinearSystem::~LinearSystem() {
	for (int i = 0; i < variableCount; i++) {
		delete[] A[i];
	}
	delete[] A;
	delete[] x;
	delete[] b;
}

void LinearSystem::loadInstance(string fileName) {
	ifstream input(fileName.c_str());
	if (input.good()) {
		input >> variableCount;
		initialize();
		for (int i = 0; i < variableCount; i++) {
			for (int j = 0; j < variableCount; j++) {
				input >> A[i][j];
			}
			input >> b[i];
		}
	}
}

void LinearSystem::initialize() {
	A = new double*[variableCount];
	x = new double[variableCount];
	b = new double[variableCount];
	for (int i = 0; i < variableCount; i++) {
		A[i] = new double[variableCount];
		x[i] = 0;
	}
}

double* LinearSystem::getEquation(int line) {
	double* equation = new double[variableCount + 1];
	for (int i = 0; i < variableCount; i++) {
		equation[i] = A[line][i];
	}
	equation[variableCount] = b[line];
	return equation;
}

bool LinearSystem::ensureNonZeroDiagonal() {
	bool ok = true;
	for (int i = 0; i < variableCount; i++) {
		if (A[i][i] == 0) {
			return false;
		}
	}
	return true;
}

void LinearSystem::printVerticalArray(wostream& output, double* array, int size, wstring label) {
	int varLine = size / 2;
	wchar_t beg, content, end;
	wstring pre;
	int cw = getColumnWidth(array, size, precision);
	output << endl;
	for (int i = 0; i < size; i++) {
		if (i == 0) {
			beg = BEGIN_BRACKET_1;
			end = END_BRACKET_1;
		} else if (i < size - 1) {
			beg = BEGIN_BRACKET_2;
			end = END_BRACKET_2;
		} else {
			beg = BEGIN_BRACKET_3;
			end = END_BRACKET_3;
		}

		if (i == varLine) {
			pre = label + L" = ";
		} else {
			pre = L"    ";
		}

		output << pre << beg << right << setw(cw) << setprecision(precision) << fixed << array[i];

		output << " " << end << endl;
	}
}

void LinearSystem::printUnformatted(wostream& output) {
	output << "A:" << endl << endl;
	for (int i = 0; i < variableCount; i++) {
		for (int j = 0; j < variableCount; j++) {
			output << A[i][j] << '\t';
		}
		output << endl;
	}
	output << endl;
	output << "--------------------------" << endl;
	output << endl;
	if (solved) {
		output << "x:" << endl << endl;
		for (int i = 0; i < variableCount; i++) {
			output << x[i] << endl;
		}
	}
	output << endl;
	output << "--------------------------" << endl;
	output << endl;
	output << "b:" << endl << endl;
	for (int i = 0; i < variableCount; i++) {
		output << b[i] << endl;
	}
}

void LinearSystem::printCompactMatrix(wostream& output) {
	int varLine = variableCount / 2;
	wchar_t beg, content, end;
	wstring pre;
	pre = BLANK;
	wchar_t inter1, inter2;
	int coeffCW = getColumnWidth(A, variableCount, precision);
	int varCW;
	if (solved) {
		varCW = getColumnWidth(x, variableCount, precision);
	} else {
		varCW = (int) floor(log10(variableCount)) + 1;
	}
	int resCW = getColumnWidth(b, variableCount, precision);

	for (int i = 0; i < variableCount; i++) {

		output << pre;

		if (i == 0) {
			beg = BEGIN_BRACKET_1;
			end = END_BRACKET_1;
		} else if (i < variableCount - 1) {
			beg = BEGIN_BRACKET_2;
			end = END_BRACKET_2;
		} else {
			beg = BEGIN_BRACKET_3;
			end = END_BRACKET_3;
		}
		if (i == varLine) {
			inter1 = MULT;
			inter2 = EQUALS;
		} else {
			inter1 = inter2 = BLANK;
		}

		output << beg;

		for (int j = 0; j < variableCount; j++) {
			output << right << setw(coeffCW) << setprecision(precision) << fixed << A[i][j];
		}

		output << " " << end << " " << inter1 << " " << beg;
		if (solved) {
			output << right << setw(varCW) << x[i];
		} else {
			output << left << setw(varCW) << "a" << getSubsFromIndex(i);
		}
		output << " " << end << " " << inter2 << " " << beg;
		output << right << setw(resCW) << b[i] << " " << end << endl;
	}

	if (solved) {
		printVerticalArray(output, x, variableCount, L"x");
	}
}

void LinearSystem::printExtendedMatrix(wostream& output) {
	int varLine = variableCount / 2;
	wchar_t beg, content, end;
	wstring pre;
	int coeffCW = getColumnWidth(A, variableCount, precision);

	for (int i = 0; i < variableCount; i++) {

		if (i == varLine) {
			pre = L"A = ";
		} else {
			pre = L"    ";
		}

		output << pre;

		if (i == 0) {
			beg = BEGIN_BRACKET_1;
			end = END_BRACKET_1;
		} else if (i < variableCount - 1) {
			beg = BEGIN_BRACKET_2;
			end = END_BRACKET_2;
		} else {
			beg = BEGIN_BRACKET_3;
			end = END_BRACKET_3;
		}

		output << beg;

		for (int j = 0; j < variableCount; j++) {
			output << right << setw(coeffCW) << setprecision(precision) << fixed << A[i][j];
		}

		output << " " << end << endl;
	}
	if (solved) {
		printVerticalArray(output, x, variableCount, L"x");
	}
	printVerticalArray(output, b, variableCount, L"b");

}

void LinearSystem::printEquationMatrix(wostream& output) {
	int varLine = variableCount / 2;
	wchar_t beg, content;
	wstring pre = L" ";
	int coeffCW = getColumnWidth(A, variableCount, precision);
	int resCW = getColumnWidth(b, variableCount, precision);
	int varCW = (int) floor(log10(variableCount));
	int columnWid = coeffCW + varCW;

	for (int i = 0; i < variableCount; i++) {

		output << pre;

		if (i == 0) {
			beg = FUNCTION_1;
		} else if (i == varLine) {
			beg = FUNCTION_2;
		} else if (i < variableCount - 1) {
			beg = FUNCTION_A;
		} else {
			beg = FUNCTION_3;
		}

		output << beg;

		for (int j = 0; j < variableCount; j++) {
			output << right << setw(columnWid) << setprecision(precision) << fixed << A[i][j] << "a" << getSubsFromIndex(j);
		}

		output << " = " << setw(resCW) << setprecision(precision) << fixed << b[i] << endl;
	}

	if (solved) {
		int xCW = getColumnWidth(x, variableCount, precision);
		output << endl << " [";
		for (int i = 0; i < variableCount - 1; i++) {
			output << left << setw(varCW) << "a" << getSubsFromIndex(i) << ", ";
		}

		output << setw(varCW) << "a" << getSubsFromIndex(variableCount - 1) << "] = [";

		for (int i = 0; i < variableCount - 1; i++) {
			output << setw(xCW) << setprecision(precision) << fixed << x[i] << ", ";
		}
		output << setw(xCW) << setprecision(precision) << fixed << x[variableCount - 1] << "]" << endl;
	}
	output << endl;
}

void LinearSystem::print(wostream& output, int mode) {
	output << endl;
	
	switch(mode) {
		case UNFORMATTED:
			printUnformatted(output);
			break;
		case COMPACT_MATRIX:
			printCompactMatrix(output);
			break;
		case EXTENDED_MATRIX:
			printExtendedMatrix(output);
			break;
		case EQUATION:
			printEquationMatrix(output);
			break;
	}

	output << endl;
}

//////////////////////////#### Gauss_Jacobi:

Gauss_Jacobi::Gauss_Jacobi(LinearSystem* system, double tolerance) {
	this -> system = system;
	this -> tolerance = tolerance;
}

Gauss_Jacobi::Gauss_Jacobi(InputConfiguration config) {
	if (config.checkConfiguration()) {
		configuration = config;
		system = new LinearSystem(config.fileName, config.precision);
		tolerance = config.tolerance;
	} else {
		printError();
	}
}

void Gauss_Jacobi::computeAllRoots() {
	setlocale(LC_ALL, "");
	printIntro(wcout);
	if (!system -> ensureNonZeroDiagonal()) {
		wcout << "The linear system you gave me isn't adequate, it must have non-zero diagonal!" << endl;
		return;
	}
	int vars = system -> getVariableCount();
	double* xValues = system -> getXValues();
	double* xPrev;
	double* line;
	double error;
	double tol = pow(10, -tolerance);
	int itCount = 0;
	do {
		itCount++;
		xPrev = copy(xValues, vars);
		for (int i = 0; i < vars; i++) {
			line = system -> getEquation(i);
			double sum = 0;
			for (int j = 0; j < vars; j++) {
				if (i != j) {
					sum += line[j] * xPrev[j];
				}
			}
			xValues[i] = (line[vars] - sum) / (line[i]);
		}
		error = computeError(xPrev);
		delete[] xPrev;
	} while (error >= tol);
	system -> setSolved(true);
	system -> print(wcout, configuration.format);
}

double computeRoot(int x_id) {

}

double Gauss_Jacobi::computeError(double* xPrev) {
	double* xValues = system -> getXValues();
	double num = 0;
	double den = 0;
	for (int i = 0; i < system -> getVariableCount(); i++) {
		num += pow((xValues[i] - xPrev[i]), 2);
		den += pow(xValues[i], 2);
	}
	return (sqrt(num) / sqrt(den));
}