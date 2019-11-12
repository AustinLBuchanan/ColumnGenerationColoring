#include "GRB.h"
#include "kgraph.h"
#include <sstream>
#include <string>

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cerr << "ERROR: Not enough arguments.";
	}
	else if (strcmp(argv[1], "degeneracy") == 0) // compute degeneracy(G)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		vector<long> ordering = g.FindDegeneracyOrdering();
	}
	else if (strcmp(argv[1], "degeneracy_coloring") == 0) // compute degeneracy_coloring(G). It will use at most degeneracy+1 colors.
	{
		KGraph g(argv[3], argv[3], argv[2]);
		vector<long> ordering = g.FindDegeneracyOrdering();
		vector<long> coloring = GreedyColoring(g, ordering);
		cout << "Is valid coloring? " << IsValidColoring(g, coloring);
	}
	else if (strcmp(argv[1], "chromatic_fractional") == 0) // compute chi_f(G)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		double chi_f = FractionalChromaticNumber(g);
		cout << "Fractional chromatic number \chi_f(G) = " << chi_f << endl;
	}
	else
	{
		cout << "ERROR: Your command is not valid." << endl;
	}
	return EXIT_SUCCESS;
}
