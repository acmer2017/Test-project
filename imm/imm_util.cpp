#ifndef IMM_UTIL_CPP
#define IMM_UTIL_CPP

#include "../mygraph.cpp"
#include <string>
#include <fstream>
#include <cstdlib>

////IMM stuff

#define HEAD_INFO

#include "sfmt/SFMT.h"
#include "head.h"

class Argument{
public:
    int k;
    string dataset;
    double epsilon;
    string model;
};

#include "graph.h"
#include "infgraph.h"
#include "imm.h"


void write_imm_input( string imm_dir,
		      tinyGraph& G, //directed, weighted
		      string icORlt = "ic" ) {

   string s_attr = imm_dir + "attribute.txt";

   node_id n = G.n;
   node_id m = G.m;

   ofstream of_atr( s_attr.c_str() );
   of_atr << "n=" << n << endl;
   of_atr << "m=" << m << endl;
   of_atr.close();

   string s_graph = imm_dir + "graph_" + icORlt + ".inf";
   ofstream of_gr( s_graph.c_str() );

   for (node_id i = 0; i < n; ++i) {
      tinyNode& I = G.adjList[i];
      for (size_t j = 0; j < I.neis.size(); ++j) {
	 of_gr << i << ' ' << I.neis[j].getId() << ' ' << static_cast< double >(I.neis[j].weight) / 255.0 << '\n';
      }
   }
}

void run_imm( string in_dir,
	      tinyGraph& G,
	      vector< node_id >& ss,
	      double& influence,
	      unsigned k,
	      string icORlt = "ic" ) {
   string scomm = "mkdir " + in_dir;
   system( scomm.c_str() );

   write_imm_input( in_dir, G, icORlt );

   Argument arg;
   arg.dataset = in_dir;
   arg.epsilon = 0.1;
   arg.k = k;
   if (icORlt == "ic")
      arg.model = "IC";
   else
      arg.model = "LT";

   string graph_file = arg.dataset + "graph_" + icORlt + ".inf";
   InfGraph g( arg.dataset, graph_file );
   if (icORlt == "ic")
      g.setInfuModel( InfGraph::IC );
   else
      g.setInfuModel( InfGraph::LT );

   ss.clear();

   g.init_hyper_graph();
   Imm::InfluenceMaximize( g, arg );

   cout << "IMM: " <<  g.InfluenceHyperGraph() << endl;

   influence = g.InfluenceHyperGraph();

   ss.assign( g.seedSet.begin(), g.seedSet.end() );

   scomm = "rm -rf " + in_dir;
   system( scomm.c_str() );

}

void run_imm_up_to_k( string in_dir,
		      tinyGraph& G,
		      vector< vector< node_id > >& ss, //seed sets, size 1 to k
		      vector< double >& profit, //influence of each seed set
		      unsigned k,
		      string icORlt = "ic" ) {
   ss.assign( k, vector< node_id >() );
   profit.assign( k, 0.0 );

   for (unsigned i = 0; i < k; ++i) {
      run_imm( in_dir,
	       G,
	       ss[ i ],
	       profit[ i ],
	       i + 1,
	       icORlt );
   }
}

#endif
