#define HEAD_INFO

#include "sfmt/SFMT.h"
#include "head.h"

class Argument{
public:
    int k;
    string dataset;
    double epsilon;
    string model;
    double T;
};

#include "graph.h"
#include "infgraph.h"
#include "imm.h"
#include <vector>
#include <iostream>

bool bin_search( InfGraph& g, Argument& arg, int first, int last, int& res ) {
  std::cout << first << ' ' << last << endl;
  g.init_hyper_graph(); 
  if (first == last) {
    arg.k = first;
    Imm::InfluenceMaximize(g, arg);
    if (g.InfluenceHyperGraph() >= arg.T) {
      res = first;
      return true;
    }
    else {
      return false;
    }
  }

  int mid = (first + last) / 2;
  arg.k = mid;
  Imm::InfluenceMaximize(g, arg);

  if ( g.InfluenceHyperGraph() >= arg.T ) {
    if ( bin_search( g, arg, first, mid - 1, res ) ) {
      return true;
    } else {
      //the answer must be mid, couldn't find a k value excluding it
      res = mid;
      return true;
    }
  } else {
    return bin_search( g, arg, mid + 1, last, res );
    
  }
}

void run_with_parameter(InfGraph &g, Argument & arg)
{
        cout << "--------------------------------------------------------------------------------" << endl;
        cout << arg.dataset << " k=" << arg.k << " epsilon=" << arg.epsilon <<   " " << arg.model << " " << arg.T << endl;

	// arg.k = 0;
	// do{
	//   arg.k = arg.k + 1;
	//   INFO( arg.k );
	//   Imm::InfluenceMaximize(g, arg);
	// } while ( g.InfluenceHyperGraph() < arg.T );

	//        INFO(g.seedSet);
	//        INFO(g.InfluenceHyperGraph());

	int res;
	bin_search( g, arg, 1, g.n, res );

	cout << "res = " << res << endl;
	arg.k = res;
	Imm::InfluenceMaximize( g, arg );
        INFO(g.seedSet);
        INFO(g.InfluenceHyperGraph());
    Timer::show();
}
void Run(int argn, char **argv)
{
    Argument arg;


    for (int i = 0; i < argn; i++)
    {
        if (argv[i] == string("-help") || argv[i] == string("--help") || argn == 1)
        {
            cout << "./tim -dataset *** -epsilon *** -k ***  -model IC|LT|TR|CONT " << endl;
            return ;
        }
        if (argv[i] == string("-dataset")) 
            arg.dataset = argv[i + 1];
        if (argv[i] == string("-epsilon")) 
            arg.epsilon = atof(argv[i + 1]);
        if (argv[i] == string("-T")) 
            arg.T = atof(argv[i + 1]);
        if (argv[i] == string("-k")) 
            arg.k = atoi(argv[i + 1]);
        if (argv[i] == string("-model"))
            arg.model = argv[i + 1];
    }
    ASSERT(arg.dataset != "");
    ASSERT(arg.model == "IC" || arg.model == "LT" || arg.model == "TR" || arg.model=="CONT");

    string graph_file;
    if (arg.model == "IC")
        graph_file = arg.dataset + "graph_ic.inf";
    else if (arg.model == "LT")
        graph_file = arg.dataset + "graph_lt.inf";
    else if (arg.model == "TR")
        graph_file = arg.dataset + "graph_tr.inf";
    else if (arg.model == "CONT")
        graph_file = arg.dataset + "graph_cont.inf";
    else
        ASSERT(false);

    InfGraph g(arg.dataset, graph_file);


    if (arg.model == "IC")
        g.setInfuModel(InfGraph::IC);
    else if (arg.model == "LT")
        g.setInfuModel(InfGraph::LT);
    else if (arg.model == "TR")
        g.setInfuModel(InfGraph::IC);
    else if (arg.model == "CONT")
        g.setInfuModel(InfGraph::CONT);
    else
        ASSERT(false);

    INFO(arg.T);

    run_with_parameter(g, arg);
}


int main(int argn, char **argv)
{
    __head_version = "v1";
    OutputInfo info(argn, argv);


    Run( argn, argv );
}


