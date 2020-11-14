// CCSubMain.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "QuadMesh.h"
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <set>


using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	//read and write data files;
	//string path = "mesh.vtk";
    //ifstream is("D:\\0Research\\Loop with NURBS boundary\\subdivision\\data\\six.obj"); //
    //ifstream is("D:\\0Research\\Loop with NURBS boundary\\subdivision\\data\\torus.obj");
	ifstream is("D:\\0Research\\Loop with NURBS boundary\\subdivision\\data\\cube.obj");
    if (!is.is_open()) {
        cout << "fail to open the file" << endl;
        return -1;
    }

	QuadMesh quad_mesh;
	quad_mesh.ReadFile(is);
	quad_mesh.GeneratMesh();

	//Find adjacent faces for each vertex;
	for(int i = 0; i<quad_mesh.m_Vertices.size(); i++){
		cout<<"Vetex "<<quad_mesh.m_Vertices[i].m_nVerid<<", its faces: ";
		vector<CHalfEdge*>  ver_he;
		
		CHalfEdge  *he0, *he1, *he2; 
		he0 = quad_mesh.m_Vertices[i].pHalfEdge;
		
		//search to the starting half-edge by finding null dual half-edge;  anticlockwise; 
		he1 = he0;		
		while( he1->dual() != NULL)
		{		
			he1 = he1->dual();			
			he1 = he1->pPrevHalfedge; 
						
			if(he1 == he0)
				break;
		}
		

		// search all edges and faces in clockwise. If two faces are missed, this code will miss some edges and faces. 	
		he1 = he1->pNextHalfedge; 
		he2 = he1; 
		cout<<he1->pFace->m_nFaceid<<", ";
		while( he1->dual() != NULL)
		{
			he1 = he1->dual();
			he1 = he1->pNextHalfedge; 
			cout<<he1->pFace->m_nFaceid<<", ";

			if(he1 == he2)
				break;
		}
		cout<<endl;
	}


	//Find adjacent faces and vertices for each edge; 
	for(int i = 0; i<quad_mesh.m_Edges.size(); i++){
		cout<<"Edge: "<<i<<" "; 
		cout<<"vertices: "<<quad_mesh.m_Edges[i].m_nTwoVer[0]<<"<->"<<quad_mesh.m_Edges[i].m_nTwoVer[1]<<", Faces: ";
		cout<<quad_mesh.m_Edges[i].pLeftHalfEdge->pFace->m_nFaceid<<", ";
		if(quad_mesh.m_Edges[i].pRightHalfEdge != NULL)
		{
			cout<<quad_mesh.m_Edges[i].pRightHalfEdge->pFace->m_nFaceid;
		}
		cout<<endl;
	}




	
	cout<<"#Vertices:"<<(int)quad_mesh.m_Vertices.size()<<", "<<endl;
	cout<<"#Face:    "<<(int)quad_mesh.m_Faces.size()<<", "<<endl;
	cout<<"#Edges:   "<<(int)quad_mesh.m_Edges.size()<<", "<<endl;
	cout<<(int)quad_mesh.mEdges.size()<<endl;
	
	
	


	
	
	// 
    system("pause");

	int i; 
	cin>>i;
    return 0;
}

