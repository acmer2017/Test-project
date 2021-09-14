// Copyright 2017 Alan Kuhnle.

// This file is part of mim.

// mim is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// mim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with mim.  If not, see <http://www.gnu.org/licenses/>.
#ifndef INFLUENCE_CPP
#define INFLUENCE_CPP
using namespace std;
using namespace mygraph;

/*
 * Conducts forward BFS on g from S.
 * returns reachable set R
 */
void forwardBFS ( tinyGraph& g, set< node_id >& S, set< node_id >& R ) {
   queue< node_id > Q; //BFS Queue
   R.clear();
   for (auto it = S.begin(); it != S.end(); ++it) {
      Q.push( *it );
   }

   while (!Q.empty()) {
      node_id curr = Q.front();
      Q.pop();

      set< node_id > neis;
      g.getOutNeis( curr, neis );
      for (auto it = neis.begin(); it != neis.end(); ++it) {
	 if (R.find( *it ) == R.end()) {
	    R.insert( *it );
	    Q.push( *it );
	 } 
      }
   }
}

void sampleGraphIC( tinyGraph& g, tinyGraph& h ) {
   h.adjList.clear();
   h.n = g.n;
   h.init_empty_graph();

   uniform_int_distribution< unsigned > dist(0, 254);
   for (node_id i = 0; i < g.n; ++i) { 
      tinyNode& I = g.adjList[ i ];
      for (size_t j = 0; j < I.neis.size(); ++j) {
	 tinyEdge& e = I.neis[j];
	 unsigned rn = dist( gen );
	 if (rn < static_cast< unsigned >( e.weight )) {
	    //Add this edge to h
	    h.add_directed_edge_immediate( i, e.getId() );
	 }
      }
   }
}

/*
 * Estimates reachability of seed set S
 * According to IC model
 * g is directed, weighted graph
 */
double monteCarloIC( tinyGraph& g, set< node_id >& S, size_t nIter = 10000 ) {
   set< node_id > R;   //Reachable set
   tinyGraph h;
   double est = 0.0;
   for (size_t i = 0; i < nIter; ++i) {
      sampleGraphIC( g, h ); //h is instance of IC from g
      forwardBFS( h, S, R );
      est += static_cast< double >( R.size() );
   }
   est /= nIter;
   return est;
}

#endif
