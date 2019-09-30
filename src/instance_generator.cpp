#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

#define CHANCE_OF_ZERO 50

using namespace std;

void generateInstance(int size) {
	if (size <= 1) {
		cout << "\nInvalid size. Size must be greater than 2.\n\n";
	} else {
		time_t rawtime = time(nullptr);
		struct tm* timeinfo;
		timeinfo = localtime (&rawtime);
		stringstream ss;
		ss << "instances/matrix_" << put_time(timeinfo, "%Y%m%d%H%M%S") << ".txt";

		string fileName =ss.str();
		ofstream output(fileName.c_str());
		output << size << endl;
		srand(time(NULL));
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (i == j) {
					output << rand() % 50 + 31 << " ";
				} else {
					int roulette = rand() % 100;
					if (roulette < CHANCE_OF_ZERO) {
						output << 0 << " ";
					} else {
						output << rand() % 10 + 1 << " ";
					}
				}
			}
			output << rand() % 100 + 21 << endl;
		}
		output.close();
	}
}

void printError(string exeName) {
	cout << "\nUsage:\t" << exeName << " --size [NUMBER OF VARIABLES]" << endl << endl;
}

int main(int argc, char** argv) {
	if (argc == 3) {
		if (string(argv[1]) == "--size") {
			int size = stoi(argv[2]);
			generateInstance(size);
			return 0;
		}
	}
	printError(string(argv[0]));
	return -1;
}