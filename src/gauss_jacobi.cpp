#include "gauss_jacobi.hpp"
#include <cmath>

#include <fstream>
#include <iomanip>
#include <omp.h>
#include <pthread.h>
#include <unistd.h>

using namespace std;

//////////////////////////         AUXILIAR FUNCTIONS  ////////////////////////
///////////////////////////////////////////////////////////////////////////////

void copy(double* copied, double* original, int vars) {
	for (int i = 0; i < vars; i++) {
		copied[i] = original[i];
	}
}

void printError() {
	wcout << "\nYou gave me incorrect parameters!\n";
	wcout << "Here, learn your commands for me:\n\n";
	wcout << left << "\nBASIC OPTIONS:\n\n";
	wcout << left << setw(35) << "--instance [INSTANCE_NAME.txt]";
	wcout << left << "With this you give me the name of the file which have the linear system you want to play with.\n";
	wcout << left << setw(35) << "--precision [DECIMAL_PLACES]";
	wcout << left << "Here you specify the precision you want me to show on the output. It doesn't really impact on the calculations...\n";
	wcout << left << setw(35) << "--tolerance [MAX_ERROR]";
	wcout << left << "And here you specify what's the maximum tolerated error. A bigger value shows better results, but it's slower!\n";
	wcout << left << "\nOUTPUT OPTIONS:\n\n";
	wcout << left << setw(35) << "--compact";
	wcout << left << "It means you want me to show the results in a compact matrix multiplication view. Don't use it if the input is too large!\n";
	wcout << left << setw(35) << "--extended";
	wcout << left << "It means you want me to show each matrix separately. Don't use it if the input is too large!\n";
	wcout << left << setw(35) << "--equation";
	wcout << left << "It means you want me to show the linear system in the analytic form. Don't use it if the input is too large!\n";
	wcout << left << setw(35) << "--results";
	wcout << left << "Use it to show only the results of the operations, that is, the x values and the residual vector. Perfect for large inputs!\n";
	wcout << left << "\nPARALLELIZATION OPTIONS:\n\n";
	wcout << left << setw(35) << "--sequential";
	wcout << left << "I'll calculate everything the way you're familiar with, everything in one big thread, no parallelization at all.\n";
	wcout << left << setw(35) << "--openmp";
	wcout << left << "I'll split the work of calculations between different and concurrent lines of processing, using the OpenMP library.\n";
	wcout << left << setw(35) << "--pthread [NUMBER_OF_THREADS]";
	wcout << left << "I'll split the work of calculations between different and concurrent lines of processing, using the POSIX Thread library.\n";
	wcout << left << "\nTIME OPTIONS:\n\n";
	wcout << left << setw(35) << "--secs";
	wcout << left << "I'll count the time in seconds.\n";
	wcout << left << setw(35) << "--millis";
	wcout << left << "I'll count the time in milliseconds.\n";
	wcout << left << setw(35) << "--micros";
	wcout << left << "I'll count the time in microseconds.\n";
	wcout << left << setw(35) << "--nanos";
	wcout << left << "I'll count the time in nanoseconds.\n\n\n";
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

wstring getTimeUnit(double resolution) {
	if (resolution == SECOND) return L"s";
	else if (resolution == MILLISECOND) return L"ms";
	else if (resolution == MICROSECOND) return L"\u03bcs";
	else if (resolution == NANOSECOND) return L"ns";
}

//////////////////////////            IMPLEMENTATIONS  ////////////////////////
///////////////////////////////////////////////////////////////////////////////

//////////////////////////#### Stopwatch:

void Stopwatch::mark() {
	if (t0 == -1) {
		t0 = clock();
	} else if (t1 == -1) {
		t1 = clock();
	} else {
		reset();
		t0 = clock();
	}
}

double Stopwatch::getElapsedTime(double resolution) {
	if (good()) {
		return double(resolution * (t1 - t0)) / CLOCKS_PER_SEC;
	} else {
		return -1;
	}
}

//////////////////////////#### InputConfiguration:

InputConfiguration::InputConfiguration() {
	format = RESULT_ONLY;
	strategy = SEQUENTIAL;
	threads = 1;
	resolution = MICROSECOND;
	precision = 2;
	fileName = "";
	tolerance = 5;
	error = true;
}

InputConfiguration::InputConfiguration(const InputConfiguration& conf2) {
	format = conf2.format;
	strategy = conf2.strategy;
	threads = conf2.threads;
	resolution = conf2.resolution;
	precision = conf2.precision;
	tolerance = conf2.tolerance;
	fileName = conf2.fileName;
	error = conf2.error;
}

bool InputConfiguration::checkConfiguration() {
	return !error;
}

//////////////////////////#### LinearSystem:

LinearSystem::LinearSystem(string fileName, int precision) {
	loadInstance(fileName);
	this -> precision = precision;
	solved = false;
}

LinearSystem::~LinearSystem() {
	delete[] x;
	delete[] b;
	delete[] allocation;
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
	allocation = new double[variableCount * variableCount];
	A = new double*[variableCount];
	x = new double[variableCount];
	b = new double[variableCount];
	for (int i = 0; i < variableCount; i++) {
		A[i] = &allocation[i * variableCount];
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

double* LinearSystem::subs() {
	double* result = new double[variableCount];
	for (int i = 0; i < variableCount; i++) {
		result[i] = 0;
		for (int j = 0; j < variableCount; j++) {
			result[i] += A[i][j] * x[j];
		}
	}
	return result;
}

double* LinearSystem::getResidualVector() {
	double* xSubs = subs();
	double* result = new double[variableCount];
	for (int i = 0; i < variableCount; i++) {
		result[i] = abs(xSubs[i] - b[i]);
	}
	return result;
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

void LinearSystem::printResults(wostream& output) {
	output << endl;
	double* residue = getResidualVector();
	printVerticalArray(output, x, variableCount, L"x");
	printVerticalArray(output, residue, variableCount, L"r");
	delete[] residue;
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
	double* xSubs = subs();
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
		output << right << setw(resCW) << xSubs[i] << " " << end << endl;
	}
	delete[] xSubs;
	if (solved) {
		printVerticalArray(output, x, variableCount, L"x");
	}
}

void LinearSystem::printExtendedMatrix(wostream& output) {
	int varLine = variableCount / 2;
	wchar_t beg, content, end;
	wstring pre;
	int coeffCW = getColumnWidth(A, variableCount, precision) + 1;

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
		double* xSubs = subs();
		printVerticalArray(output, x, variableCount, L"x");
		printVerticalArray(output, xSubs, variableCount, L"s");
		delete[] xSubs;
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
	int sig = 3;

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

		bool firstPassed = false;
		output << beg;

		for (int j = 0; j < variableCount; j++) {
			output << internal << setw(sig);
			if (A[i][j] != 0) {
				if (!firstPassed) {
					if (A[i][j] < 0) {
						output << "-";
					} else {
						output << "";
					}
					firstPassed = true;
				} else {
					if (A[i][j] < 0) {
						output << "-";
					} else {
						output << "+";
					}
				}
			} else {
				//output << "";
			}
			if (A[i][j] != 0) {
				output << right << setw(columnWid) << setprecision(precision) << fixed << abs(A[i][j]) << "a" << getSubsFromIndex(j);
			} else {
				//output << "";
			}
		}

		output << "  = " << setw(resCW) << setprecision(precision) << fixed << b[i] << endl;
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
	
	switch (mode) {
		case RESULT_ONLY:
			printResults(output);
			break;
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
		xPrev = new double[system -> getVariableCount()];
	} else {
		printError();
	}
}

void Gauss_Jacobi::printIntro(wostream& output) {
	wchar_t i = L"\u1d62"[0];
	wchar_t j = L"\u2c7c"[0];
	wchar_t sum = L"\u2211"[0];
	wchar_t uneq = L"\u2260"[0];
	output << "######## Hello! This is the GAUSS-JACOBI method for Linear Systems!\n";
	output << "###################################################################\n\n";
	output << "    This work was made by the team FHRR, hope you enjoy it as we do ;)\n\n";
	output << "    The Gauss-Jacobi method is very simple, it consists in isolating each variable from the system,\n";
	output << "then using each of the obtained equations to calculate new values for their respective variables.\n";
	output << "For each iteration it does, the value of the variables gets more and more close to the solution of\n";
	output << "the system, and the number of iterations is directly proportional to the adopted tolerance, so the\n";
	output << "greater the precision desired, the more iterations will be needed.\n\n";
	output << "    There's two conditions to assure the convergence of the method: the main diagonal must be free\n";
	output << "from zero elements, and the matrix must be diagonally dominant, that is:\n\n";
	output << L"\t\t\tFor each line i, |a" << i << i << "| > " << sum << "|a";
	output << i << j << "|; i" << uneq << "j\n\n";
}

void Gauss_Jacobi::printFullHeader(wostream& output) {
	output << "    The input file you gave me had this linear system:\n";
	system -> setPrecision(0);
	system -> print(output, EQUATION);
	system -> setPrecision(configuration.precision);
	output << "After applying the Gauss-Jacobi method, this is what I've got for you:\n";
}

void Gauss_Jacobi::printBasicHeader(wostream& output) {
	wchar_t dot = L"\u22c5"[0];
	wchar_t note = L"\u2042"[0];
	output << "    As per your request, I'll give you only the results from the operations.\n";
	output << "    First you have the x vector, then you have the residue vector ";
	output << "r = |A" << dot << "x - b|\n\n    ";
	output << note << " Note that the pipes here means absolute value, not norm.\n";
}

void Gauss_Jacobi::findSolution() {
	setlocale(LC_ALL, "");
	if (configuration.checkConfiguration()) {
		printIntro(wcout);
		if (configuration.format == RESULT_ONLY) {
			printBasicHeader(wcout);
		} else {
			printFullHeader(wcout);
		}
		if (!system -> ensureNonZeroDiagonal()) {
			wcout << "The linear system you gave me isn't adequate, it must have non-zero diagonal!" << endl;
			return;
		}
		switch (configuration.strategy) {
			case SEQUENTIAL:
				computeRootsSequential();
				break;
			case PARALLEL_OPENMP:
				computeRootsParallelOpenMP();
				break;
			case PARALLEL_PTHREAD:
				computeRootsParallelPThread();
				break;
			default:
				return;
		}
		system -> setSolved(true);
		system -> print(wcout, configuration.format);
		wstring unit = getTimeUnit(configuration.resolution);
		wcout << "\nThe time taken to do everything was " << stopwatch.getElapsedTime(configuration.resolution) << unit << ".\n\n";
	}
}

void Gauss_Jacobi::computeRootsSequential() {
	int vars = system -> getVariableCount();
	double* xValues = system -> getXValues();
	double** matrix = system -> getA();
	double* secHand = system -> getB();
	double error;
	double tol = pow(10, -tolerance);
	double totalTime = 0;
	double sum;
	int itCount = 1;
	int i, j;
	stopwatch.mark();
	do {
		copy(xPrev, xValues, vars);
		for (i = 0; i < vars; i++) {
			sum = 0;

			for (j = 0; j < vars; j++) {
				if (i != j) {
					sum += matrix[i][j] * xPrev[j];
				}
			}
			xValues[i] = (secHand[i] - sum) / (matrix[i][i]);
		}
		error = computeError();
		if (itCount <= 500) {
			stopwatch.mark();
			totalTime += stopwatch.getElapsedTime(configuration.resolution);
			wcerr << "SOLUTION #" << itCount << " ||| ~TIME = " << (totalTime / itCount) << "ms" << ".\n";
			itCount++;
			for (int i = 0; i < vars; i++) {
				xValues[i] = 0;
			}
			copy(xPrev, xValues, vars);
			error = 1;
			stopwatch.mark();
		}
	} while (error >= tol);
	delete[] xPrev;
	stopwatch.mark();
}

void Gauss_Jacobi::computeRootsParallelOpenMP() {
	int vars = system -> getVariableCount();
	double* xValues = system -> getXValues();
	double** matrix = system -> getA();
	double* secHand = system -> getB();
	double error;
	double tol = pow(10, -tolerance);
	double sum;
	double totalTime = 0;
	int itCount = 1;
	int i, j, tNumber, tid;
	int ll, ul, mod;
	stopwatch.mark();
	copy(xPrev, xValues, vars);
	#pragma omp parallel shared(vars, xValues, matrix, secHand, mod, tNumber, error) private (sum, tid, i, j, ll, ul)
	{
		#pragma omp master
		{
			tNumber = omp_get_num_threads();
			mod = vars % tNumber;
		}
		tid = omp_get_thread_num();
		#pragma omp barrier
		if (tid < mod) {
			ll = tid * ((vars + 1) / tNumber);
			ul = (tid + 1) * ((vars + 1) / tNumber);
		} else {
			ll = tid * (vars / tNumber) + mod;
			ul = (tid + 1) * (vars / tNumber) + mod;
		}
		do {
			for (i = ll; i < ul; i++) {
				sum = 0;
				for (j = 0; j < vars; j++) {
					if (i != j) {
						sum += matrix[i][j] * xPrev[j];
					}
				}
				xValues[i] = (secHand[i] - sum) / (matrix[i][i]);
			}
			//wcerr << "TA INDO OU NN? :::" << tid << ":::\t";
			#pragma omp barrier
			#pragma omp master
			{
				error = computeError();
				if (error >= tol) {
					copy(xPrev, xValues, vars);
				} else if (itCount <= 500) {
					stopwatch.mark();
					totalTime += stopwatch.getElapsedTime(configuration.resolution);
					wcerr << "SOLUTION #" << itCount << " ||| ~TIME = " << (totalTime / itCount) << "ms" << ".\n";
					itCount++;
					for (int i = 0; i < vars; i++) {
						xValues[i] = 0;
					}
					copy(xPrev, xValues, vars);
					error = 1;
					stopwatch.mark();
				}
			}
			#pragma omp barrier
		} while (error >= tol);
	}
	delete[] xPrev;
	stopwatch.mark();
}

void Gauss_Jacobi::computeRootsParallelPThread() {
	pthread_t* threadIDs = new pthread_t[configuration.threads];
	pthread_barrier_init(&barrier, NULL, configuration.threads);
	int rc;
	for (int i = 0; i < configuration.threads; i++) {
		//wcerr << "CREATING THREAD #" << i << "\n";
		rc = pthread_create(&threadIDs[i], NULL, Gauss_Jacobi::worker_wrapper, (void*) new ThreadParameters(this, i));
        if (rc){
            wcout << "ERROR: pthread_create() returned the code " << rc << "\n\n";
            exit(-1);
        }
	}
	for (int i = 0; i < configuration.threads; i++) {
		pthread_join(threadIDs[i], NULL);
	}
	delete[] xPrev;
	delete[] threadIDs;
	running = false;
	pthread_barrier_destroy(&barrier);
}

void* Gauss_Jacobi::worker_wrapper(void* param) {
	((ThreadParameters*) param) -> instance -> gauss_jacobi_worker(((ThreadParameters*) param) -> threadID, ((ThreadParameters*) param) -> itCounter);
}

void Gauss_Jacobi::gauss_jacobi_worker(int threadID, int& itCount) {
	int vars = system -> getVariableCount();
	int tid = (int) threadID;
	int tNumber = configuration.threads;
	double* xValues = system -> getXValues();
	double** matrix = system -> getA();
	double* secHand = system -> getB();
	double tol = pow(10, -tolerance);
	double sum;
	double totalTime = 0;
	int i, j;
	int lowerLimit, upperLimit;
	int mod = vars % tNumber;
	if (tid < mod) {
		lowerLimit = tid * ((vars + 1) / tNumber);
		upperLimit = (tid + 1) * ((vars + 1) / tNumber);
	} else {
		lowerLimit = tid * (vars / tNumber) + mod;
		upperLimit = (tid + 1) * (vars / tNumber) + mod;
	}

	if (tid == 0) {
		copy(xPrev, system -> getXValues(), system -> getVariableCount());
	}

	pthread_barrier_wait(&barrier);

	if (tid == 0) {
		stopwatch.mark();
	}

	//wcerr << "THREAD #" << tid << "\n";
	do {
		//wcerr << "a";
		for (i = lowerLimit; i < upperLimit; i++) {
			sum = 0;
			for (j = 0; j < vars; j++) {
				if (i != j) {
					sum += matrix[i][j] * xPrev[j];
				}
			}
			xValues[i] = (secHand[i] - sum) / (matrix[i][i]);
		}
		pthread_barrier_wait(&barrier);
		if (tid == 0) {
			running = computeError() > tol;
			if (running) {
				copy(xPrev, system -> getXValues(), system -> getVariableCount());
			} else if (itCount <= 500) {
				stopwatch.mark();
				totalTime += stopwatch.getElapsedTime(configuration.resolution);
				wcerr << "SOLUTION #" << itCount << " ||| ~TIME = " << (totalTime / itCount) << "ms" << ".\n";
				itCount++;
				for (int i = 0; i < vars; i++) {
					xValues[i] = 0;
				}
				copy(xPrev, xValues, vars);
				running = true;
				stopwatch.mark();
			}
		}
		pthread_barrier_wait(&barrier);
	} while (running);
	if (tid == 0) {
		stopwatch.mark();
	}
	pthread_exit(NULL);
}

double Gauss_Jacobi::computeError() {
	double* xValues = system -> getXValues();
	double num = 0;
	double den = 0;
	for (int i = 0; i < system -> getVariableCount(); i++) {
		num += pow((xValues[i] - xPrev[i]), 2);
		den += pow(xValues[i], 2);
	}
	return (sqrt(num) / sqrt(den));
}