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
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;
using namespace mygraph;

void write_model( unsigned i, string name ) {
   ofstream ofile( name.c_str() );
   ofile << i % 2;
   ofile.close();
}

void write_ov( unsigned i, unsigned k, string name, double beta, unsigned n ) {
   ofstream ofile( name.c_str() );

   unsigned N = beta * n;

   for (unsigned node = 0; node < N; ++node) {
      for (unsigned j = 0; j < k; ++j) {
	 if (j != i) {
	    ofile << node << ' ' << j << ' ' << node << endl;
	 }
      }
   }
   
   
   ofile.close();
}

void process_layer( unsigned n,  //ER n
		    double p,    //ER p
		    unsigned i,  //layer index
		    double beta, //overlap beta
		    unsigned k,
		    string fname ) {   //number of layers ) {
  tinyGraph g;
  g.erdos_renyi_directed( n, p, 25 );

  string name = fname + "/layer" + to_string( i );
      
  g.write_directed_edge_list( name );

  name = fname + "/layer" + to_string( i ) + "model";
  write_model( i, name );
  
  name = fname + "/layer" + to_string( i ) + "ov";

  write_ov( i, k, name, beta, n );
}


int main( int argc, char** argv ) {
   if ( argc < 4 ) {
      cout << "Usage: " << argv[0] << " <number of nodes in each layer> <er edge prob p> <number of layers> <fraction of overlapping users> <directory to write files>" << endl;
      exit(1);
   }
   
   unsigned n = stoi( argv[1] );
   double p = stod( argv[2] );
   unsigned k = stoi( argv[3] );
   double beta = stod( argv[ 4 ] );
   string fname ( argv[5] );

   thread* layerThreads = new thread[ k ];
   
   for (unsigned i = 0; i < k; ++i) {
     layerThreads[i] = thread(
			      process_layer,
			      n,
			      p,
			      i,
			      beta,
			      k,
			      fname );

   }

   for (unsigned i = 0; i < k; ++i) {
     layerThreads[i].join();
   }

   delete [] layerThreads;
   
   return 0;
}
