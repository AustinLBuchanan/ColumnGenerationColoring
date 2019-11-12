#include "gurobi_c++.h"
#include "KGraph.h"
#include <algorithm>

double FractionalChromaticNumber(KGraph &g);
vector< vector<long> > initMaxIndSets(KGraph &g);

vector<long> GreedyColoring(KGraph &g);
vector<long> GreedyColoring(KGraph &g, vector<long> &ordering);

bool IsValidColoring(KGraph &g, vector<long> &coloring);
void makeMaximal(KGraph &g, vector<long> &S); // make S a maximal ind set by adding to it
