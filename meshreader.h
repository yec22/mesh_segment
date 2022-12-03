#pragma once
#include <fstream>
#include "dsa.h"

class MeshReader{
    public:
        MeshReader() {}
        void read(std::string fname, Graph& graph){
            std::ifstream plyfile(fname);
            std::string tmp;
            
            // parse ply header
            int f_num, v_num;
            while(plyfile >> tmp) {
                if(tmp == std::string("end_header")) break;
                if(tmp == std::string("vertex")){
                    plyfile >> v_num;
                }
                if(tmp == std::string("face")){
                    plyfile >> f_num;
                }
            }

            // parse vertices
            double x, y, z;
            for (int i = 0; i < v_num; ++i){
                plyfile >> x >> y >> z;
                graph.vertices.push_back(Vertex(x, y, z));
            }

            // parse faces
            int n, v1, v2, v3;
            for (int i = 0; i < f_num; ++i){
                plyfile >> n >> v1 >> v2 >> v3;
                graph.edges.push_back(Edge(v1, v2, i));
                graph.edges.push_back(Edge(v2, v3, i));
                graph.edges.push_back(Edge(v3, v1, i));
                graph.faces.push_back(Face(graph.vertices[v1], graph.vertices[v2], graph.vertices[v3], v1, v2, v3));
            }

            std::cout << "load " << fname << " successfully! ";
            std::cout << "vertex_num: " << v_num << ", face_num: " << f_num << std::endl;
        }
};