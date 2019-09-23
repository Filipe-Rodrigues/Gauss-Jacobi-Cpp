#include "gauss_seidel.cpp"

using namespace std;

int main () {
	setlocale(LC_ALL, "");
	LinearSystem ls("instances/test1.txt");
	ls.print(wcout, COMPACT_MATRIX);
	return 0;
}