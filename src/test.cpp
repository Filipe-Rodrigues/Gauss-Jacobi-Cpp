#include "gauss_seidel.cpp"

using namespace std;

void printError() {

}

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

int main (int argc, char** argv) {
	testLS(getConfiguration(argc, argv));
	return 0;
}