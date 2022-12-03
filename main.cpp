#include "meshreader.h"
#include "meshwriter.h"
#include "dsa.h"
using namespace std;

int main(int argc, char* argv[]) {
	string input_fname = argv[1];
    string output_fname = argv[2];
	
	cout << endl;
	for(int i = 0; i < 75; ++i) cout << "*"; cout << endl;
	MeshReader reader;
	MeshWriter writer;
	Graph graph;
	DualGraph dualgraph;
	GraphSolver graphsolver;

	reader.read(input_fname, graph);
	dualgraph.Assign_Weight(graph);
	dualgraph.Shortest_Path(graph);
	graphsolver.Init(&graph, &dualgraph, 0);
	graphsolver.Solve();
	writer.write(output_fname, graph);
	
	for(int i = 0; i < 75; ++i) cout << "*"; cout << endl;
	cout << endl;
	cout << "done!!!" << endl;
	return 0;
}