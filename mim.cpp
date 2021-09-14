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
#include "mygraph.cpp"
#include "influence.cpp"
#include "multiplex.cpp"
#include "mckp.cpp"

#include <iostream>
#include <set>
#include <unistd.h>
#include <sstream>
#include <ctime>

//imm functions
#include "imm/imm_util.cpp"

using namespace std;
using namespace mygraph;

void print_help() {
   cout << "Options:" << endl;
   cout << "\t-M [multiplex dir]\n"
	<< "\t-I (run ISF)\n"
      	<< "\t-K (run KSN)\n"
      	<< "\t-k [number of seeds]\n"
	<< "\t-x [max number of threads]\n";
        
}

struct npair {
   pair< node_id, node_id > P;
   double rank;
};

class compare_pair {
public:
   bool operator()( npair& a, npair& b ) {
      return (a.rank < b.rank);
   }
};

typedef pair<node_id, node_id > MultiNode;


void algKSN( Multiplex& M, size_t k, set< MultiNode >& S, unsigned nThreads = 5 ) {
   vector< vector< vector< node_id > > > ss (M.nLayers, vector< vector<node_id> >());
   vector< vector< double > > profit( M.nLayers, vector< double >() );

   unsigned layersProcessed = 0;

   while( layersProcessed < M.nLayers ) {
     unsigned layersIter = M.nLayers - layersProcessed;
     if (layersIter > nThreads)
       layersIter = nThreads;
   
     thread* immThreads = new thread[ layersIter ];
   
     for (size_t i = layersProcessed; i < layersProcessed + layersIter; ++i) {
       string diri = "tmp_imm" + to_string( i ) + "/";
       tinyGraph& G = M.Layers[i];
       string icORlt;
       switch (M.layerModels[i]) {
       case IC:
	 icORlt = "ic";
	 break;
       case LT:
	 icORlt = "lt";
	 break;
       default:
	 icORlt = "ic";
	 break;
       }
       immThreads[i - layersProcessed] = thread(
						run_imm_up_to_k,
						diri,
						ref(G),
						ref(ss[i]),
						ref(profit[i]),
						k,
						icORlt );
     }

     for (size_t i = 0; i < layersIter; ++i) {
       immThreads[i].join();
     }

     delete [] immThreads;

     layersProcessed += layersIter;
   }

   //Now, for each layer, ss has a set of k seed sets with corresponding profits in profit
   stringstream sInput;
   sInput << k << endl;
   for (size_t i = 0; i < M.nLayers; ++i) {
      for (size_t j = 0; j < k; ++j) {
	 sInput << profit[i][j] << ' ';
      }
      sInput << endl;
   }

   vector< vector< mckp::item > > sys;

   double b;
   mckp::read_input( sInput, sys, b );
   
   for (unsigned i = 0; i < sys.size(); ++i) {
      for (unsigned j = 0; j < sys[i].size(); ++j) {
	 cerr << sys[i][j].c << ' ';
      }
      cerr << endl;
   }

   vector< unsigned > choice;
   mckp::BM( choice, sys, b );

   cout << "size of choice: " << choice.size() << endl;
   for (unsigned i = 0; i < choice.size(); ++i) {
      cout << i + 1 << ' ' << choice[i] << endl;
      //Picking seed set choice[i] of layer i = ss[i][choice[i]]
      if (choice[i] > 0) {
	 for (unsigned j = 0; j < ss[i][choice[i] - 1].size(); ++j) {
	    MultiNode P;
	    P.first = i;
	    P.second = ss[i][choice[i]-1][j];
	    S.insert( P );
	 }
      }
   }
   
}


void algISF( Multiplex& M, size_t k, set< MultiNode >& S ) {
   S.clear();

   priority_queue< npair, vector< npair >, compare_pair > Q;

   map< MultiNode, bool > b_valid;
   
   for (unsigned j = 0; j < M.nLayers; ++j) {
      size_t nj = M.Layers[j].n;
      for (unsigned l = 0; l < nj; ++l) {
	 npair tmp;
	 pair< node_id, node_id > N( j, l );
	 tmp.P = N;
	 tmp.rank = M.nNodesAllLayers; //upper bound on marginal gain

	 Q.push( tmp );

	 b_valid[ N ] = false;
      }
   }

   for (size_t i = 0; i < k; ++i) {
      npair pmax = Q.top();
      Q.pop();

      double curr_act = M.monteCarloInfluence( S );
      
      while (true) {
	 if ( b_valid[ pmax.P ] ) {
	    S.insert( pmax.P );
	    break;
	 } else {
	    S.insert( pmax.P );
	    double marg = M.monteCarloInfluence( S ) - curr_act;
	    pmax.rank = marg;
	    Q.push( pmax );
	    b_valid[ pmax.P ] = true;
	    S.erase( pmax.P );
	 }
      }
   }

   //   double curr_act = M.monteCarloInfluence( S );
   //   M.logg << "ISF terminating, final activation: " << curr_act << endL;
}

enum Algo {ISF, KSN};

int main(int argc, char** argv) {
   extern char *optarg;
   //  extern int optint, optopt;

   if (argc == 1) {
      print_help();
      return 1;
   }

   int c;
   string s_arg;

   Multiplex M;
   Algo A;
   size_t k = 10;

   unsigned nThreads = 1;
   
   Logger logg;
   
   while ((c = getopt( argc, argv, ":M:IKk:x:") ) != -1) {
      switch(c) {
      case 'x':
	 s_arg.assign( optarg );
	 nThreads = stoi( s_arg );
	 break;
      case 'k':
	 s_arg.assign( optarg );
	 k = stoi( s_arg );
	 break;
      case 'I':
	 A = ISF;
	 break;
      case 'K':
	 A = KSN;
	 break;
      case 'M':
	 s_arg.assign( optarg );
	 M.init( s_arg );
	 break;
      case '?':
	 print_help();
	 return 1;
	 break;
      }
   }
   
   logg << "Running Alg. ";
   switch( A ) {
   case KSN:
      logg << "KSN";
      break;
   case ISF:
      logg << "ISF";
      break;
   }

   logg << " with k = " << k << endL;
   set< MultiNode > S;



   std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
   clock_t tStart = clock();
   
   switch( A ) {
   case KSN:
     algKSN( M, k, S, nThreads );
      break;
   case ISF:
      algISF( M, k, S );
      break;
   }

   

   double cpuTime = ( (double) clock() - tStart ) / CLOCKS_PER_SEC;
   
   std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

   
   size_t timeMillis = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
   logg << "Time taken (Wall, ms): " << timeMillis << endL;
   logg << "Time taken (CPU, s): " << cpuTime << endL;
   logg << "Seed set (layer id, node id): ";
   for (auto i = S.begin(); i != S.end(); ++i) {
      logg << "(" << i->first << "," << i->second << "),";
   }
   logg << endL;
   logg << "Calculating influence of seed set by 10000 Monte Carlo iterations... " << endL;
   

   //   logg << "Computing activation of seed set:" << endL;
   
   logg << "Total activation: " << M.monteCarloInfluence( S ) << endL;
   
   return 0;
}
