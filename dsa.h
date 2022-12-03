// Reference: https://github.com/fornorp/Mesh-Segmentation
#pragma once
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

typedef std::pair<int, double> PII;

const double EPS = 1e-10;
const double INF = 1e10;
const double ETA = 0.2;
const double DELTA = 0.8;
const double PROB_THR = 0.05;
const double REP_DIST_RATIO = 0.1;
const double AVG_DIST_RATIO = 0.2;
const int MAX_ITER = 10;
const int FUZZY = -1234;
const int MAX_DEPTH = 2;
const int K = 2;
const int A_FLAG = 0;
const int B_FLAG = 1;
const int OTHER_FLAG = 777;

static int LABEL_BASE = 0;
static double GLOBAL_MAX_DIST;
static double GLOBAL_AVG_DIST;

struct Vec3d {
	double x[3];
    
    Vec3d() {}
	Vec3d(const Vec3d& v){
        x[0] = v.x[0]; x[1] = v.x[1]; x[2] = v.x[2];
    }
    Vec3d(double x_, double y_, double z_){
        x[0] = x_; x[1] = y_; x[2] = z_;
    }

	// arithmetic
    Vec3d operator+(const Vec3d& v) const{
        return Vec3d(x[0]+v.x[0], x[1]+v.x[1], x[2]+v.x[2]);
    } 
    Vec3d operator-(const Vec3d& v) const{
        return Vec3d(x[0]-v.x[0], x[1]-v.x[1], x[2]-v.x[2]);
    } 
    Vec3d operator*(double d) {
        return Vec3d(x[0]*d, x[1]*d, x[2]*d);
    } 
    Vec3d operator/(double d) {
        return Vec3d(x[0]/d, x[1]/d, x[2]/d);
    } 
}; 

inline double dot(const Vec3d& v1, const Vec3d& v2) {return v1.x[0]*v2.x[0]+v1.x[1]*v2.x[1]+v1.x[2]*v2.x[2];}

inline double len(const Vec3d& v) {return sqrt(dot(v, v));}

inline Vec3d normalize(const Vec3d& v) {return Vec3d(v.x[0]/len(v), v.x[1]/len(v), v.x[2]/len(v));}

inline Vec3d cross(const Vec3d& v1, const Vec3d& v2) {
    return Vec3d(v1.x[1] * v2.x[2] - v1.x[2] * v2.x[1],
	             v1.x[2] * v2.x[0] - v1.x[0] * v2.x[2],
	             v1.x[0] * v2.x[1] - v1.x[1] * v2.x[0]);
}

struct Vertex {
	Vec3d p; // position
    Vertex() {}
    Vertex(const Vertex& v) {p=v.p;}
	Vertex(double x, double y, double z) : p(x, y, z) {}
}; 

struct Edge {
	int v_id[2];
    int f_id;
	Edge() {}
	Edge(int v1, int v2, int f) {v_id[0] = v1; v_id[1] = v2; f_id = f;}
}; 

struct Neighbour {
    int f_id;
	double Ang_Dist, Geod, Angle;
    double dist;

    Neighbour() {}
    Neighbour(const Neighbour& n) {f_id = n.f_id; Ang_Dist = n.Ang_Dist; Geod = n.Geod; Angle = n.Angle; dist = n.dist;}
	Neighbour(int f, double ang_dist, double geod, double angle) {f_id = f; Ang_Dist = ang_dist; Geod = geod; Angle = angle;}
}; 

struct Face {
    int label;
    int v_id[3];
    Vec3d center;
    Vec3d normal;
    std::vector<Neighbour> neighbours;

    Face() {}
    Face(const Face& f) {
        for (int i = 0; i < 3; ++i)
			v_id[i] = f.v_id[i];
        label = f.label;
        center = f.center;
        normal = f.normal;
        neighbours = f.neighbours;
    };
    Face(Vertex& v1, Vertex& v2, Vertex& v3, int id1, int id2, int id3) {
        v_id[0] = id1; v_id[1] = id2; v_id[2] = id3;
        center = (v1.p + v2.p + v3.p) / 3.0;
		normal = normalize(cross(v2.p - v1.p, v3.p - v1.p));
		label = 0;
    }
}; 

inline bool isConvex(const Face& f1, const Face& f2) {return dot(f1.normal, f2.center - f1.center) < EPS;}

inline double Normal_Angle(const Face& f1, const Face& f2) {return acos(dot(f1.normal, f2.normal));}

inline double Cal_Ang_Dist(const Face& f1, const Face& f2) {
    double eta = isConvex(f1, f2) ? ETA : 1.0;
	return eta * (1 - dot(f1.normal, f2.normal));
}

inline double Cal_Geo_Dist(const Face& f1, const Face& f2, const Vertex& v1, const Vertex& v2) {
    Vec3d axis = v2.p - v1.p;
    Vec3d f1v1 = f1.center - v1.p;
    Vec3d f2v1 = f2.center - v1.p;
    double angle1 = acos(dot(f1v1, axis) / (len(f1v1) * len(axis)));
	double angle2 = acos(dot(f2v1, axis) / (len(f2v1) * len(axis)));
	return len(f1v1) * len(f1v1) + len(f2v1) * len(f2v1) - 2 * len(f1v1) * len(f2v1) * cos(angle1 + angle2);
}

class Heap {
    public:
        std::vector<PII> h;
        std::vector<int> table;
        std::vector<bool> st;
        int tail;

        Heap() {}
        void Resize(int n) {table.resize(n); st.resize(n); h.resize(n+1);}
        void Swap(int k, int l);
        void Up(int k);
        void Down(int k);
        void Insert(PII val);
        void Modify(int k, PII val);
        int Pop();
};

class Graph {
    public:
        double avg_Ang_Dist, avg_Geod;
        std::vector<Vertex> vertices;
        std::vector<Edge> edges;
        std::vector<Face> faces;

        Graph() {}
};

class DualGraph {
    public:
        std::vector<std::vector<double> > w; // adjacency list
        std::vector<std::vector<double> > dis; // adjacency matrix
        std::vector<std::vector<int> > ne; // adjacency list
        Heap heap; // binary heap
        
        DualGraph() {}
        void Assign_Weight(Graph& graph);
        void Shortest_Path(Graph& graph);
};

struct FlowEdge {
	int id;
	double cap;
	double flow;
    FlowEdge() {}
    FlowEdge(int _id, double _cap, double _flow) {id = _id; cap = _cap; flow = _flow;}
};

struct FlowState
{
	int id;
	int front;
	double inc;
    FlowState() {}
    FlowState(int _id, int _front, double _inc) {id = _id; front = _front; inc = _inc;}
};

class GraphSolver {
    public:
        int REP_A, REP_B;
        int level, sub_n;
        double avg_dist;
        std::vector<double> P_A;
        std::vector<double> P_B;
        std::vector<int> sub_v;
        std::vector<std::vector<FlowEdge> > flow_net;
        int S, T;
        Graph* g;
        DualGraph* dg;

        GraphSolver() {}
        void Init(Graph* _g, DualGraph* _dg, int l, std::vector<int>& s);
        void Find_Max_Dist();
        void Solve();
        void Assign_Prob();
        void Assign_Label();
        void Update_Rep();
        void Graph_Cut();
        void Init_Flow(std::vector<int>& p);
        void Max_Flow(std::vector<int>& p);
        void DFS_Solve();
        double Cal_Avg_Dist();
};