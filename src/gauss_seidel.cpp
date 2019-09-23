#include "gauss_seidel.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

//////////////////////////         AUXILIAR FUNCTIONS  ////////////////////////
///////////////////////////////////////////////////////////////////////////////

int getColumnWidth(double* values, int size) {
	int width = 0;
	double numLength;
	for (int i = 0; i < size; i++) {
		numLength = floor(log10(abs(values[i]))) + 1;
		if (values[i] < 0) {
			numLength++;
		}
		if (abs(values[i]) - floor(abs(values[i])) >= 0) {
			numLength += 3;
		}
		if (numLength > width) {
			width = numLength;
		}
	}
	return width + 1;
}

int getColumnWidth(double** values, int size) {
	int width = 0;
	double numLength;
	for (int i = 0; i < size; i++) {
		numLength = getColumnWidth(values[i], size) - 1;
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

//////////////////////////            IMPLEMENTATIONS  ////////////////////////
///////////////////////////////////////////////////////////////////////////////

LinearSystem::LinearSystem(string fileName) {
	loadInstance(fileName);
	solved = false;
	x = NULL;
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

void LinearSystem::printVerticalArray(wostream& output, double* array, int size, wstring label) {
	int varLine = size / 2;
	wchar_t pre, beg, content, end;
	int cw = getColumnWidth(array, size);
	for (int i = 0; i < size; i++) {

		output << pre;

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
			inter1 = MULT;
			inter2 = EQUALS;
		} else {
			inter1 = inter2 = BLANK;
		}

		output << beg;

		for (int j = 0; j < variableCount; j++) {
			output << right << setw(coeffCW) << setprecision(2) << fixed << A[i][j];
		}

		output << " " << end << " " << inter1 << " " << beg << " ";
		if (solved) {
			output << right << setw(varCW) << x[i];
		} else {
			output << left << setw(varCW) << "a" << getSubsFromIndex(i);
		}
		output << " " << end << " " << inter2 << " " << beg;
		output << right << setw(resCW) << b[i] << " " << end << endl;
	}
}

void LinearSystem::printUnformatted(wostream& output) {
	output << "A:" << endl;
	for (int i = 0; i < variableCount; i++) {
		for (int j = 0; j < variableCount; j++) {
			cout << A[i][j] << '\t';
		}
		cout << endl;
	}
	if (solved) {
		output << "x:" << endl;
		for (int i = 0; i < variableCount; i++) {
			cout << x[i] << endl;
		}
	}
	output << "b:" << endl;
	for (int i = 0; i < variableCount; i++) {
		cout << b[i] << endl;
	}
}

void LinearSystem::printCompactMatrix(wostream& output) {
	int varLine = variableCount / 2;
	wchar_t pre, beg, content, end;
	pre = BLANK;
	wchar_t inter1, inter2;
	int coeffCW = getColumnWidth(A, variableCount);
	int varCW;
	if (solved) {
		varCW = getColumnWidth(x, variableCount);
	} else {
		varCW = (int) floor(log10(variableCount)) + 1;
	}
	int resCW = getColumnWidth(b, variableCount);

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
			output << right << setw(coeffCW) << setprecision(2) << fixed << A[i][j];
		}

		output << " " << end << " " << inter1 << " " << beg << " ";
		if (solved) {
			output << right << setw(varCW) << x[i];
		} else {
			output << left << setw(varCW) << "a" << getSubsFromIndex(i);
		}
		output << " " << end << " " << inter2 << " " << beg;
		output << right << setw(resCW) << b[i] << " " << end << endl;
	}
}

void LinearSystem::printExtendedMatrix(wostream& output) {
	int varLine = variableCount / 2;
	wchar_t pre, beg, content, end;
	int coeffCW = getColumnWidth(A, variableCount);

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
			output << right << setw(coeffCW) << setprecision(2) << fixed << A[i][j];
		}

		output << " " << end << endl;
	}
	if (solved) {
		output << endl;
		printVerticalArray(output, x, variableCount);
	}
	output << endl;
	printVerticalArray(output, b, variableCount);
	output << endl;

}

void LinearSystem::printEquationMatrix(wostream& output) {

}

void LinearSystem::print(wostream& output, int mode) {
	output << endl;
	
	switch(mode) {
		case COMPACT_MATRIX:
			printCompactMatrix(output);
			break;
		case EXTENDED_MATRIX:

	}

	output << endl;
}