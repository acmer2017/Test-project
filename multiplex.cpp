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

#ifndef MULTIPLEX_CPP
#define MULTIPLEX_CPP

#include "mygraph.cpp"
#include <string>
#include <fstream>
#include <cstdlib>

using namespace std;
using namespace mygraph;

enum Model { IC = 0, LT = 1, DLT = 2 };

bool file_exists (const std::string& name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }   
}

class Multiplex {
public:
   Logger logg;
   
   unsigned nLayers; //number of layers
   unsigned nNodesAllLayers; //no overlap
   
   vector< tinyGraph > Layers;
   vector< Model > layerModels;
   
   vector< map< node_id, vector< pair< node_id, node_id > > > > Overlap;

   /*
    * Overlap file contains lines
    * i j k
    * where i is node_id mapped to layer j, node k
    */
   void addOverlap( string fname ) {
      map< node_id, vector< pair< node_id, node_id > > > tmp;

      ifstream ifile( fname.c_str() );
      node_id i, j, k;


      
      while (ifile >> i) {
	 vector< pair< node_id, node_id > >& v = tmp[ i ];
	 ifile >> j;
	 ifile >> k;
	 pair< node_id, node_id > pTmp;
	 pTmp.first = j;
	 pTmp.second = k;
	 v.push_back( pTmp );
      }

      Overlap.push_back( tmp );
   }

   //public:

   Multiplex() {
      
   }
   
   Multiplex( string fname ) {
      init( fname );
   }
   
   void init( string input ) {
      logg << "Constructing multiplex from " << input << "..." << endL;

      input += "/layer";
      
      ifstream ifile (input.c_str());
      
      nLayers = 0;

      string layerFname;
      string layerOverlapName;
      Model layerModel;
      unsigned tmp;
      string lfname = input + to_string( nLayers );
      nNodesAllLayers = 0;
      while (file_exists( lfname )) {
	 ++nLayers;
	 tinyGraph g;
	 g.read_directed_edge_list( lfname );
	 Layers.push_back( g );
	 nNodesAllLayers += g.n;
	 lfname += "ov";
	 addOverlap( lfname );

	 lfname = input + to_string( nLayers - 1 ) + "model";
	 ifstream ifile( lfname.c_str() );
	 ifile >> tmp;
	 layerModel = (Model) tmp;
	 ifile.close();
	 layerModels.push_back( layerModel );

	 logg << "Model of layer " << nLayers - 1 << " is ";
	 switch (layerModel) {
	 case IC:
	    logg << "IC" << endL;
	    break;
	 case LT:
	    logg << "LT" << endL;
	    break;
	 case DLT:
	    logg << "DLT" << endL;
	    break;
	 }
	 
	 lfname = input + to_string( nLayers );
      }

      logg << "Multiplex has " << nLayers << " layers." << endL;
   }

   bool activate( vector< tinyGraph >& detLayers,
		  unsigned char edgeWeight,
		  pair< node_id, node_id > next,
		  set < pair< node_id, node_id > >& R ) {

      if (R.find( next ) != R.end() ) {
	 return false; //already activated
      }
      
      unsigned layer = next.first;
      switch (layerModels[ layer ]) {
      case IC:
	 return true;
	 break;
      case LT:
	 tinyNode& N = detLayers[ layer ].adjList[ next.second ];
	 N.weight += edgeWeight;

	 if (N.weight >= N.thresh)
	    return true;

	 return false;
	 
	 break;
      }

      return false;
   }
   
   unsigned forwardProp( vector< tinyGraph >& detLayers,
			 set< pair< node_id, node_id > >& S ) {
      unsigned totAct = 0;
      
      queue< pair< node_id, node_id > > Q;
      set< pair< node_id, node_id > > R;
      for (auto it = S.begin(); it != S.end(); ++it) {
	 //	 if (R.find( *it ) == R.end()) {
	 //	    ++totAct;
	    Q.push( *it );
	    // R.insert( *it );
	    //Get all overlapping nodes for this one
	    //We know none of them have been activated yet
	    //	    vector< pair< node_id, node_id > >& vOv = (Overlap[ it->first ])[ it->second ];
	    //	    for (size_t j = 0; j < vOv.size(); ++j) {
	    //    R.insert( vOv[ j ] );
	    //   Q.push( vOv[j] );
	    //}
	    //	 }
      }

      while (!Q.empty()) {
	 pair< node_id, node_id > curr = Q.front();
	 Q.pop();

	 
	 //detLayers[ curr.first ].getOutNeis( curr.second, neis );


	 if (activate( detLayers, static_cast<unsigned char>(255), curr, R )) {
	    ++totAct;

	    vector< pair< node_id, node_id > >& vOv = (Overlap[ curr.first ])[ curr.second ];
	    for (size_t j = 0; j < vOv.size(); ++j) {
	       R.insert( vOv[ j ] );
	       Q.push( vOv[j] );
	    }
	 }

	 //cerr << curr.first << ' ' << curr.second << endl;
	 
	 vector< tinyEdge >& neis = detLayers[curr.first].adjList[curr.second].neis;
	 
	 for (auto it = neis.begin(); it != neis.end(); ++it) {
	    pair < node_id, node_id > pTmp;
	    pTmp.first = curr.first;
	    pTmp.second = it->getId();

	    if (activate( detLayers, it->weight, pTmp, R )) {
	       ++totAct;
	       
	       R.insert( pTmp );
	       Q.push( pTmp );

	       //Get all overlapping nodes for this one
	       //We know none of them have been activated yet
	       vector< pair< node_id, node_id > >& vOv = (Overlap[ pTmp.first ])[ pTmp.second ];
	       for (size_t j = 0; j < vOv.size(); ++j) {
		  R.insert( vOv[ j ] );
		  Q.push( vOv[j] );
	       }
	    }
	 }
      }

      return totAct;
   }

   void sampleMultiplex( vector< tinyGraph >& detLayers ) {
      detLayers.clear();
      for (unsigned i = 0; i < nLayers; ++i) {
	 tinyGraph h;
	 switch (layerModels[i]) {
	 case IC:
	    {
	       tinyGraph& g = Layers[i];
	       h.n = g.n;
	       h.init_empty_graph();
	       uniform_int_distribution< unsigned > dist(0, 255);
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
	    break;
	 case LT:
	    {
	       tinyGraph& g = Layers[i];
	       h.assign( g );
	       uniform_int_distribution< unsigned > dist(0, 255);
	       for (node_id i = 0; i < h.n; ++i) { 
		  unsigned rn = dist( gen );
		  h.adjList[i].thresh = rn;
		  h.adjList[i].weight = 0;
	       }
	    }
	    break;

	 case DLT:
	    {
	       tinyGraph& g = Layers[i];
	       h.assign( g );
	    }
	    break;
	 }
	 detLayers.push_back( h );
      }
   }
   
   /*
    * Estimates reachability of seed set S
    * According to multiplex models
    */
  double monteCarloInfluence( set< pair< node_id, node_id > >& S, size_t nIter = 10000, bool bVerbose = false ) {
    if (bVerbose) {
      logg << "Estimating influence with " << nIter << " MC samples..." << endL;
    }
      double est = 0.0;
      for (size_t i = 0; i < nIter; ++i) {
	 //	cout << "\r                       \r" << i;
	 //	cout.flush();
	  vector< tinyGraph > detLayers;
	sampleMultiplex( detLayers );
	est += static_cast< double >( forwardProp( detLayers, S ) );
      }

      est /= nIter;

      return est;
   }
   
};

#endif
