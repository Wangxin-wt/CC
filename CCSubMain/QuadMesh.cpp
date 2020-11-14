#include "StdAfx.h"
#include "QuadMesh.h"
#include <cassert>
#include <set>
#include <iostream>
#include <algorithm>



typedef map<pair<index_t, index_t >, index_t > directedEdge2indexMap_t;
index_t directedEdge2faceIndex(const directedEdge2indexMap_t& de2fi, index_t vi, index_t vj) {
	if (de2fi.empty()) {
		cout << "fail to build mesh: de2fi == empty()" << endl;
		return -1;
	}

	auto it = de2fi.find(make_pair(vi, vj));
	//没有对应面的he返回-1

	if(it == de2fi.end()) {
		if (de2fi.find(make_pair(vj, vi)) == de2fi.end()) {
			cout << "fail to build mesh: de2fi.find(make_pair(vj, vi)) == de2fi.end()" << endl;
			return -1;
		}
		return -1;
	}
	return it->second;
}




int unorderedEdgesFromQuads(const vector<quad_t> &quads, vector<edge_t> &edgesOut) {
	typedef set<pair<index_t , index_t > > edgeSet_t;       //边集，防止重复建边
	edgeSet_t edges;
	for (index_t t = 0; t < quads.size(); ++t) { //ij, jk, kl, li
		edges.insert(make_pair(min(quads[t].i(), quads[t].j()), max(quads[t].i(), quads[t].j())));
		edges.insert(make_pair(min(quads[t].j(), quads[t].k()), max(quads[t].j(), quads[t].k())));
		edges.insert(make_pair(min(quads[t].k(), quads[t].l()), max(quads[t].k(), quads[t].l())));
		edges.insert(make_pair(min(quads[t].l(), quads[t].i()), max(quads[t].l(), quads[t].i())));
	}
	edgesOut.resize(edges.size());
	int e = 0;
	for (auto it = edges.begin(); it != edges.end(); it++, e++) {
		edgesOut.at(static_cast<unsigned long long int>(e)).start() = it->first;
		edgesOut.at(static_cast<unsigned long long int>(e)).end() = it->second;
		//cout<<"edge:" << it->first<<", "<<it->second<<endl;
	}
	return 0;
}


QuadMesh::QuadMesh(void)
{

}


QuadMesh::~QuadMesh(void)
{

}



void QuadMesh::GeneratMesh(void)
{
	int     id = 0, vid = -1;
	double  lengthedge;
	CVertex temp_ver;
	vector<vector<CHalfEdge*>>  point_to_ver;

	vector<CEdge> tedges;
	tedges.clear();
	m_Vertices.clear();	
	m_Edges.clear();
	m_Faces.clear();

	
	//push back vertices from gvertices into m_Vertices;
	for(int i=0; i<(int)mPoints.size(); i++)
	{
		temp_ver.m_nVerid =  mPoints[i].index;
		temp_ver.m_dVerx  =  mPoints[i].x;
		temp_ver.m_dVery  =  mPoints[i].y;
		temp_ver.m_dVerz  =  mPoints[i].z;
		temp_ver.pHalfEdge = NULL;

        m_Vertices.push_back(temp_ver);
	}

	point_to_ver.resize(m_Vertices.size()); 
	
	for(int faceid=0; faceid<(int)mQuads.size(); faceid++) //face by face;
	{
		CFace *face = new CFace;
		face->m_nFaceid = faceid;
		face->m_bMark = 0;
		
		for(int i=0; i<4; i++){
			face->m_nGver[i] = mQuads[faceid].v[i];			
		}
					
		CHalfEdge  *he[4]; 
		CVertex    *ve[4];
		for(int i=0; i<4; i++){
			ve[i] = new CVertex;
			he[i] = new CHalfEdge; 
		}
		
		////maintain half-edge structure;
		for(int i=0; i<4; i++){
			vid = mQuads[faceid].v[i]; //ids for four vertices;			
			ve[i]->m_dVerx = mPoints[vid].x;
			ve[i]->m_dVery = mPoints[vid].y;
			ve[i]->m_dVerz = mPoints[vid].z;
			
			ve[i]->m_nVerid = vid; 
			ve[i]->pHalfEdge = he[i]; 
			he[i]->pVertex = ve[i];
			he[i]->pFace = face;			
		}
		
		for(int i=0; i<4; i++){
			he[i]->pNextHalfedge = he[(i+1)%4]; //build connection;
			he[i]->pPrevHalfedge = he[(i+3)%4];
		}
		
		// attach half-edge to every vertex; Nice!
		for(int i=0; i<4; i++)
		{
			id = mQuads[faceid].v[i];
			point_to_ver[id].push_back(he[i]); 
		}
		
		face->pHalfedge = he[0];
			
		// attach half-edge to edge; to assign null half-edge
		bool flag = false;
		for(int i=0; i<4; i++)
		{
			int id_1 = he[i]->pVertex->m_nVerid;
			int id_2 = he[i]->pPrevHalfedge->pVertex->m_nVerid;
			int j = 0;
			id = mQuads[faceid].v[i]; 
			flag = false;

			//For every vertex, search from all its half-edges
			for(j=0; j < (int)point_to_ver[id].size(); j++) 
			{
				if(point_to_ver[id][j]->pNextHalfedge->source()->m_nVerid == id_1 && 
					point_to_ver[id][j]->pNextHalfedge->pVertex->m_nVerid == id_2 )
				{
					flag = true; break;
				}
				if(point_to_ver[id][j]->pNextHalfedge->source()->m_nVerid == id_2 && 
					point_to_ver[id][j]->pNextHalfedge->pVertex->m_nVerid == id_1 )
				{
					flag=true; break;
				}
			}

			CEdge *pe = new CEdge;
			if(flag == true)
			{					
				lengthedge = sqrt(pow(he[i]->source()->m_dVerx - he[i]->target()->m_dVerx,2)+
								pow(he[i]->source()->m_dVery - he[i]->target()->m_dVery,2));
				pe->m_dLength = lengthedge;
				pe->m_nTwoVer[0] = id_1; 
				pe->m_nTwoVer[1] = id_2;
				pe->pLeftHalfEdge = point_to_ver[id][j]->pNextHalfedge;
				pe->pRightHalfEdge = he[i];
				point_to_ver[id][j]->pNextHalfedge->pEdge = pe;
				he[i]->pEdge = pe;					
			}
			else
			{				
				lengthedge = sqrt(pow(he[i]->source()->m_dVerx - he[i]->target()->m_dVerx,2)+
								pow(he[i]->source()->m_dVery - he[i]->target()->m_dVery,2));
				pe->m_dLength = lengthedge;
				pe->m_nTwoVer[0] = id_1; 
				pe->m_nTwoVer[1] = id_2;
				pe->pLeftHalfEdge = he[i];
				pe->pRightHalfEdge = NULL;
				he[i]->pEdge = pe;  				
			}			
			tedges.push_back(*pe);
		}
		
		
		m_Faces.push_back(*face);
		

	}
	

	for(int i=0; i< (int)m_Vertices.size(); i++)
	{
		if(point_to_ver[i].size()!=0)
			m_Vertices[i].pHalfEdge = point_to_ver[i][0];
	}

	//delete repeated CEdges
	map <pair<int, int>, int> mpair;
	for(int ii=0; ii<tedges.size(); ii++)
	{
		mpair.insert(make_pair(make_pair(min(tedges[ii].m_nTwoVer[0],tedges[ii].m_nTwoVer[1]), max(tedges[ii].m_nTwoVer[0],tedges[ii].m_nTwoVer[1])),ii));		
	}
	cout<<"here"<<endl;
	for (auto it = mpair.begin(); it != mpair.end(); it++) {
		//cout<<it->second<<", ";
		m_Edges.push_back(tedges[it->second]);
	}

	//maintain RightHalfEdge from LeftHalfEdge for m_Edges;
	for(int ii=0; ii<m_Edges.size(); ii++)
	{
		m_Edges[ii].pRightHalfEdge = m_Edges[ii].pLeftHalfEdge->dual();
	}

}


void QuadMesh::ReadFile(ifstream& is) {
    int pcount = 0, tricount = 0, quadcount = 0; 
    string s;
    char c;

    cout<<"Read data file----begin"<<endl;
    while (is.get(c)) {
        switch (c) {
			// There is a bug if the number of vertices and faces is not indicated in obj files
		case '#':
			break;
        /*    
		case '#':
                getline(is, s);
                if (s.find("Vertices") != string::npos) {
                    mPoints.resize(size_t (stoll(s.substr(11))));
                }
                if (s.find("Faces") != string::npos) {
                    //triangles.resize(size_t (stoll(s.substr(8))));
					mQuads.resize(size_t (stoll(s.substr(8))));
                }
                break;
			*/
            case 'v':
                point_t p;
                p.index = pcount;
				is >> p.x >> p.y >> p.z;
                is.get();
                mPoints.push_back(p);
				pcount++;
                break;
			//Triangle case
			/*
            case 'f':
                triangle_t triangle;
                int x, y, z;
                is >> x >> y >> z;
                is.get();
                if (!(x <= mPoints.size() && y <= mPoints.size() && z <= mPoints.size())) {
                    cout << "x="  << x << "y=" << y << "z=" << z << endl;
                    cout << "maxSize=" << mPoints.size() << endl;
                    exit(-1);
                }
                if (tricount == 18429) {
                    cout << endl;
                }
                triangle.v[0] = mPoints[x - 1];
                triangle.v[1] = mPoints[y - 1];
                triangle.v[2] = mPoints[z - 1];
                triangles[tricount++] = triangle;
                break;
			*/
			case 'f':
                quad_t quad;
				quad.id = quadcount;
                int p00, p10, p11, p01;
                is >> p00 >> p10 >> p11 >> p01;
                is.get();
				/*
                if (!(p00 <= mPoints.size() && p01 <= mPoints.size() && p10 <= mPoints.size() && p11 <= mPoints.size())) {
                    cout << "x="  << p00 << "y=" << p10 << "z=" << p11 <<"w=" << p01 << endl;
                    cout << "maxSize=" << mPoints.size() << endl;
                    exit(-1);
                }
                if (quadcount == 18429) {
                    cout << endl;
                }
				*/
				// Shit, problem from obj file, face: v0,v1,v2,v3, 
				// vertices are indexed from 1, not from 0, need to -1; 
                quad.v[0] = p00-1;
                quad.v[1] = p10-1;
                quad.v[2] = p11-1;
				quad.v[3] = p01-1;

                //mQuads[quadcount++] = quad;
				mQuads.push_back(quad);
				quadcount++;
                break;
            default:
                break;
        }
    }    
	cout<<"read data file---ok!"<<endl; 
    //unorderedEdgesFromTriangles(triangles, edges);
	unorderedEdgesFromQuads(mQuads, mEdges);
    //mesh.build(mPoints, edges, triangles);	
	//mesh.build(mPoints, edges, quads);

}

bool IsSameEdge(CEdge e1, CEdge e2){
	bool flag = false; 
	if (e1.m_nTwoVer[0] == e2.m_nTwoVer[0] && e1.m_nTwoVer[1] == e2.m_nTwoVer[1])
		flag = true;
	if (e1.m_nTwoVer[1] == e2.m_nTwoVer[0] && e1.m_nTwoVer[0] == e2.m_nTwoVer[1])
		flag = true;
	return flag;
}
