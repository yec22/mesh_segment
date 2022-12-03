#pragma once
#include <fstream>
#include "dsa.h"

class MeshWriter {
    public:
        MeshWriter() {}
        void write(std::string fname, Graph& graph){
            std::ofstream plyfile(fname);

            // write ply header            
            int f_num, v_num;
            f_num = graph.faces.size();
            v_num = graph.vertices.size();
            plyfile << "ply" << std::endl;
	        plyfile << "format ascii 1.0" << std::endl;
	        plyfile << "element vertex " << v_num << std::endl;
	        plyfile << "property float x" << std::endl;
	        plyfile << "property float y" << std::endl;
	        plyfile << "property float z" << std::endl;
	        plyfile << "element face " << f_num << std::endl;
	        plyfile << "property list uchar int vertex_index" << std::endl;
	        plyfile << "property uint8 red" << std::endl;
	        plyfile << "property uint8 green" << std::endl;
	        plyfile << "property uint8 blue" << std::endl;
	        plyfile << "end_header" << std::endl;

            // write vertices
            for(int i = 0; i < v_num; ++i){
                plyfile << graph.vertices[i].p.x[0] << " ";
                plyfile << graph.vertices[i].p.x[1] << " ";
                plyfile << graph.vertices[i].p.x[2] << std::endl;
            }

            // write faces
            for(int i = 0; i < f_num; ++i){
                plyfile << "3 ";
                plyfile << graph.faces[i].v_id[0] << " ";
                plyfile << graph.faces[i].v_id[1] << " ";
                plyfile << graph.faces[i].v_id[2] << " ";
                // color
                int label = graph.faces[i].label;
		        plyfile << 60 * (label % 4 + 2) << " " << 80 * ((label + 1) % 3 + 2) << " " << 50 * ((label + 2) % 5 + 2) << std::endl;
            }

            plyfile.close();
            std::cout << "export " << fname << " successfully! ";
            std::cout << "vertex_num: " << v_num << ", face_num: " << f_num << std::endl;
        }
};