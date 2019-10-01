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
			} else if (arg == "--instance") {
				if (++count < argc) {
					config.fileName = string(argv[count]);
					count++;
				}
			} else if (arg == "--precision") {
				if (++count < argc) {
					config.precision = stoi(argv[count]);
					count++;
				}
			} else {
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
	Gauss_Jacobi gj(config);
	gj.computeRootsSequential();
}

int main (int argc, char** argv) {
	#ifdef _WIN32
	_setmode(_fileno(stdout), 0x00040000);
	system("cls");
	#elif __linux
	system("reset");
	#endif
	//testLS(getConfiguration(argc, argv));
	testGJ(getConfiguration(argc, argv));
	return 0;
}