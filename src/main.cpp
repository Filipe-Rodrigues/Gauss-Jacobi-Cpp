#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif
#include "gauss_jacobi.cpp"
#include <cstdlib>

using namespace std;

InputConfiguration getConfiguration(int argc, char** argv) {
	InputConfiguration config;
	if (argc > 1) {
		int count = 1;
		while (count < argc) {
			string arg(argv[count]);
			if (arg == "--unformatted") {
				config.format = UNFORMATTED;
				count++;
			} else if (arg == "--compact") {
				config.format = COMPACT_MATRIX;
				count++;

			} else if (arg == "--extended") {
				config.format = EXTENDED_MATRIX;
				count++;
			} else if (arg == "--equation") {
				config.format = EQUATION;
				count++;
			} else if (arg == "--results") {
				config.format = RESULT_ONLY;
				count++;
			} else if (arg == "--sequential") {
				config.strategy = SEQUENTIAL;
				count++;
			} else if (arg == "--openmp") {
				config.strategy = PARALLEL_OPENMP;
				count++;
			} else if (arg == "--pthread") {
				if (++count < argc) {
					config.strategy = PARALLEL_PTHREAD;
					config.threads = stoi(argv[count++]);
				}
			} else if (arg == "--secs") {
				config.resolution = SECOND;
				count++;
			} else if (arg == "--millis") {
				config.resolution = MILLISECOND;
				count++;
			} else if (arg == "--micros") {
				config.resolution = MICROSECOND;
				count++;
			} else if (arg == "--nanos") {
				config.resolution = NANOSECOND;
				count++;
			} else if (arg == "--instance") {
				if (++count < argc) {
					config.fileName = string(argv[count]);
					count++;
					config.error = false;
				}
			} else if (arg == "--precision") {
				if (++count < argc) {
					config.precision = stoi(argv[count]);
					count++;
				}
			} else if (arg == "--tolerance") {
				if (++count < argc) {
					config.tolerance = stoi(argv[count]);
					count++;
				}
			} else {
				config.error = true;
				break;
			}
		}
	}
	return config;
}

void testLS(InputConfiguration config) {
	if (config.checkConfiguration()) {
		setlocale(LC_ALL, "");
		LinearSystem ls(config.fileName, config.precision);
		ls.print(wcout, config.format);
	} else {
		printError();
	}
}

void testGJ(InputConfiguration config) {
	if (config.checkConfiguration()) {
		#ifdef _WIN32
		_setmode(_fileno(stdout), 0x00040000);
		system("cls");
		#elif __linux
		system("clear");
		#endif
		Gauss_Jacobi gj(config);
		gj.findSolution();
	} else {
		printError();
	}
}

int main (int argc, char** argv) {
	//testLS(getConfiguration(argc, argv));
	testGJ(getConfiguration(argc, argv));
	return 0;
}