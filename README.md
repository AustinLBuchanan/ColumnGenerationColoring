# Column Generation for Graph Coloring
A quick implementation of column generation for graph coloring. The IP formulation is:

min \sum_s x_s
s.t. \sum_{s:i \in s} x_s >= 1
     x_s \in {0,1} for all maximal independent sets s.
     
The implementation given here only solves the LP relaxation. It begins by computing a degeneracy coloring of the graph. Each color class is extended to a maximal independent set s. This gives an initial set of variables: x_s for s \in S'. The LP (restricted to these columns S') is solved. New columns are added as needed until the restricted LP gives a solution that is optimal for the full LP. The pricing problem is a maximum-weight independent set problem and is solved via the usual IP formulation:

max \sum_i w_i x_i 
s.t. x_i + x_j <= 1 for all edges {i,j}
     x_i \in {0,1} for all vertices i.

This implementation given here is not great; I did it just as a coding exercise in an afternoon. The links below offer the theory and a much better implementation.

Mehrotra and Trick paper:
https://pubsonline.informs.org/doi/abs/10.1287/ijoc.8.4.344

Held, Cook, and Sewell paper:
https://link.springer.com/article/10.1007/s12532-012-0042-3

Held, Cook, and Sewell implementation:
https://code.google.com/archive/p/exactcolors/
