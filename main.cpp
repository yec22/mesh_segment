#include "meshreader.h"
#include "meshwriter.h"
#include "dsa.h"
#include <ctime>
using namespace std;

int main(int argc, char* argv[]) {
	string input_fname = argv[1];
    string output_fname = argv[2];
	
	cout << endl;
	for(int i = 0; i < 75; ++i) cout << "*"; cout << endl;
	clock_t start = clock();

	MeshReader reader;
	MeshWriter writer;
	Graph graph;
	DualGraph dualgraph;
	GraphSolver graphsolver;
	vector<int> null;

	reader.read(input_fname, graph);
	dualgraph.Assign_Weight(graph);
	dualgraph.Shortest_Path(graph);
	graphsolver.Init(&graph, &dualgraph, 0, null);
	graphsolver.Solve();
	writer.write(output_fname, graph);
	
	clock_t end = clock();
	for(int i = 0; i < 75; ++i) cout << "*"; cout << endl;
	cout << endl;
	cout << "mesh decomposition finished!!!" << endl;
	cout << "time cost: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}