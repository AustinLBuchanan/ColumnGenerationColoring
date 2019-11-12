#include "GRB.h"

vector< vector<long> > initMaxIndSets(KGraph &g)
{
	vector<long> ordering = g.FindDegeneracyOrdering();
	vector<long> coloring = GreedyColoring(g, ordering);

	long numColors = 0;
	for (long i = 0; i < g.n; ++i)
		numColors = max(numColors, coloring[i] + 1);

	vector< vector<long> > S_prime(numColors);
	for (long i = 0; i < g.n; ++i)
	{
		long color = coloring[i];
		S_prime[color].push_back(i);
	}

	for (long color = 0; color < numColors; ++color)
	{
		makeMaximal(g, S_prime[color]);
	}

	return S_prime;
}

void makeMaximal(KGraph &g, vector<long> &S)
{
	vector<bool> available(g.n, true);

	// mark S and all neighbors as unavailable
	for (long i = 0; i < S.size(); ++i)
	{
		long s = S[i];
		for (long j = 0; j < g.degree[s]; ++j)
		{
			long v = g.adj[s][j];
			available[v] = false;
		}
		available[s] = false;
	}

	// add to S
	for (long i = 0; i < g.n; ++i)
	{
		if (!available[i]) continue;

		available[i] = false;
		S.push_back(i);
		for (long j = 0; j < g.degree[i]; ++j)
		{
			long v = g.adj[i][j];
			available[v] = false;
		}
	}
	if (count(available.begin(), available.end(), true) > 0)
	{
		cerr << "ERROR: there are still available vertices to add to ind set S." << endl;
	}

	// sort S, just in case
	sort(S.begin(), S.end());
}

double FractionalChromaticNumber(KGraph &g)
{
	vector< vector<long> > S_prime = initMaxIndSets(g);
	double chi_f = S_prime.size();

	try {
		cerr << "Creating initial master problem" << endl;

		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		GRBModel master = GRBModel(env);
		vector<GRBVar> x;// = vector<GRBVar>(S_prime.size());// = master.addVars(S_prime.size(), GRB_CONTINUOUS);
		master.update();

		for (long S = 0; S < S_prime.size(); S++)
			x.push_back(master.addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS));

		master.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		master.update();

		vector<GRBConstr> coverConstr(g.n); // = vector<GRBConstr>(g.n);
		for (long i = 0; i < g.n; ++i)
		{
			GRBLinExpr expr = 0;
			for (long s = 0; s < S_prime.size(); ++s)
			{
				for (long p = 0; p < S_prime[s].size(); ++p)
				{
					long v = S_prime[s][p];
					if (i == v) expr += x[s];
				}
			}
			coverConstr[i] = master.addConstr(expr >= 1);
		}
		master.optimize();
		chi_f = min(chi_f, master.get(GRB_DoubleAttr_ObjVal));

		vector<double> w(g.n);
		for (long i = 0; i < g.n; ++i)
			w[i] = coverConstr[i].get(GRB_DoubleAttr_Pi);

		cout << "Creating initial subproblem" << endl;
		vector<long> indset;
		GRBModel subproblem = GRBModel(env);
		GRBVar *y = subproblem.addVars(g.n, GRB_BINARY);
		subproblem.update();

		for (long i = 0; i < g.n; i++)
			y[i].set(GRB_DoubleAttr_Obj, w[i]);
		subproblem.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		subproblem.set(GRB_DoubleParam_MIPGap, 0.);
		subproblem.update();

		for (long i = 0; i < g.n; i++)
		{
			for (long j = 0; j < g.degree[i]; j++)
			{
				long v = g.adj[i][j];
				if (i < v) continue; //add edge constraints once
				subproblem.addConstr(y[i] + y[v] <= 1);
			}
		}
		subproblem.optimize();

		for (long i = 0; i < g.n; i++)
			if (y[i].get(GRB_DoubleAttr_X) > 0.5)
				indset.push_back(i);

		makeMaximal(g, indset);
		double weight = 0;
		for (long p = 0; p < indset.size(); ++p)
		{
			long v = indset[p];
			weight += w[v];
		}
		cout << "weight of MWIS = " << weight << endl;

		while (weight > 1)
		{
			// add new column to master, and re-solve
			S_prime.push_back(indset);
			x.push_back(master.addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS));
			for (long p = 0; p < indset.size(); ++p)
			{
				long v = indset[p];
				master.chgCoeff(coverConstr[v], x[x.size() - 1], 1);
			}
			master.optimize();
			chi_f = min(chi_f, master.get(GRB_DoubleAttr_ObjVal));
			cout << "master obj = " << chi_f << ", \t";

			// update obj coefficients in subproblem, and re-solve
			for (long i = 0; i < g.n; ++i)
			{
				w[i] = coverConstr[i].get(GRB_DoubleAttr_Pi);
				y[i].set(GRB_DoubleAttr_Obj, w[i]);
			}
			subproblem.optimize();

			indset.clear();
			for (long i = 0; i < g.n; i++)
				if (y[i].get(GRB_DoubleAttr_X) > 0.5)
					indset.push_back(i);

			makeMaximal(g, indset);
			weight = 0;
			for (long p = 0; p < indset.size(); ++p)
			{
				long v = indset[p];
				weight += w[v];
			}
			cout << " subproblem obj = " << weight << endl;
		}

		// solve IP over the current columns to get an upper bound on chi(G)
		for (long s = 0; s < x.size(); ++s)
			x[s].set(GRB_CharAttr_VType, GRB_BINARY);
		master.set(GRB_DoubleParam_MIPGap, 0.);
		master.optimize();
		double UB = master.get(GRB_DoubleAttr_ObjVal);
		cout << chi_f << "= chi_f(G) <= chi(G) <= " << UB << endl;

	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return chi_f;
}

vector<long> GreedyColoring(KGraph &g)
{
	vector<long> ordering(g.n);
	for (long i = 0; i < g.n; ++i)
		ordering[i] = i;
	return GreedyColoring(g, ordering);
}

vector<long> GreedyColoring(KGraph &g, vector<long> &ordering)
{
	// starting from back of ordering, color each vertex with smallest available color from {0, 1, 2, ..., n-1 }
	vector<bool> available(g.n, true);
	vector<long> coloring(g.n, -1);			// -1 denotes not colored yet

	for (long i = g.n - 1; i >= 0; --i)
	{
		long v = ordering[i];

		// what colors have been taken by v's neighbors?
		for (long j = 0; j < g.degree[v]; ++j)
		{
			long w = g.adj[v][j];
			long color = coloring[w];
			if (color >= 0)
				available[color] = false;
		}

		// color v with smallest available color
		for (long color = 0; coloring[v] == -1; ++color)
			if (available[color])
				coloring[v] = color;

		// reset available colors for the next iteration
		for (long j = 0; j < g.degree[v]; ++j)
		{
			long w = g.adj[v][j];
			long color = coloring[w];
			if (color >= 0)
				available[color] = true;
		}
	}
	long numColors = 0;
	for (long i = 0; i < g.n; ++i)
		numColors = max(numColors, coloring[i] + 1);
	cout << "Number of colors used = " << numColors << endl;

	return coloring;
}

bool IsValidColoring(KGraph &g, vector<long> &coloring)
{
	for (long i = 0; i < g.n; ++i)
	{
		for (long j = 0; j < g.degree[i]; ++j)
		{
			long v = g.adj[i][j];
			if (coloring[i] == coloring[v])
				return false;
		}
	}
	return true;
}
