#pragma once

#include <iostream>
//#include "mesh_types.h" //提供edge, triangle, quad;
#include "StructDefs.h"
#include <utility>
#include <vector>
#include <map>
#include <cassert>
#include <unordered_set>
#include <fstream>      // std::ifstream

using namespace std;

#define SAFEDELETE(p)    if(p) { delete p;    p=NULL;}
#define SAFEDELETES(p)   if(p) { delete []p;  p=NULL;}
#define TOLERANCE 1.0e-13

int unorderedEdgesFromQuads(const vector<quad_t>& quads, vector<edge_t>& edgesOut);
class CHalfEdge;
class CVertex;

//basic structure for quad_meshes
class CVertex{                     
public:    //attributes
	int        	m_nVerid;            //id
    int        	m_nvLevel;
	double     	m_dVerx, m_dVery, m_dVerz;    //coordinates
    CHalfEdge*  pHalfEdge;           //half edge
	bool       	m_bDealWith;	     //deal or not?
	int			m_nValency;
	
	CVertex()
	{
		m_nVerid = -1;
		m_nvLevel = 0;
		pHalfEdge = NULL;
		m_bDealWith = false;
		m_nValency = -1;
	}

	CVertex(int id, double vx, double vy, double vz, CHalfEdge *hep)
	{
		m_nVerid = id;
		m_dVerx = vx;
		m_dVery = vy;
		m_dVerz = vz;
		pHalfEdge = hep;
		m_bDealWith = false;
	
	}
	CVertex(const CVertex &v)
	{
		if(this != &v) 
		{
			m_nVerid = v.m_nVerid;
			m_dVerx = v.m_dVerx;
			m_dVery = v.m_dVery;
			m_dVerz = v.m_dVerz;
			pHalfEdge = v.pHalfEdge;
			m_bDealWith = v.m_bDealWith;
			m_nvLevel = v.m_nvLevel;
		}

	}
	CVertex & operator=(const CVertex &v)
	{
		m_nVerid = v.m_nVerid;
		m_dVerx = v.m_dVerx;
		m_dVery = v.m_dVery;
		m_dVerz = v.m_dVerz;
		pHalfEdge = v.pHalfEdge;
		m_bDealWith = v.m_bDealWith;
		m_nvLevel = v.m_nvLevel;
		m_nValency = v.m_nValency;
        return *this;
     } 
	
	~CVertex(){};
	

public:   //operation
     //bool IsBoundaryVertex(double x1,double x2,double y1,double y2); //lie on or not the boundary of rectangle of x1,x2,y1,y2
     //void VertexInHalfEdge(vector<CHalfEdge*> &lst_hes);	
	 //char VertexType();

	 //int LRHalfEdge(vector<CHalfEdge*> lst_hes,char type);
	 //int UDHalfEdge(vector<CHalfEdge*> lst_hes,char type);
	 //int LRUDHalfEdge(vector<CHalfEdge*> lst_hes,char type);
	 
};


class CFace{
public:
	int             m_nFaceid;			//face id;
	int				m_nGver[4];			//four vertices;
	int             m_nFaceLevel;		//face level, if needed;	
	bool          	m_bMark;            //sign to indicate refined or not
	bool          	IsSubdivid;    
	CHalfEdge  		*pHalfedge;			//any half-edge

public:
	CFace()
	{
	    m_nFaceid = -1;	    
		m_bMark = false;
		pHalfedge = NULL;
		m_nFaceLevel = 0;
		IsSubdivid = 0;
	}

	CFace & operator=(const CFace &cell1)
	{
		m_nFaceid = cell1.m_nFaceid;
		m_bMark   = cell1.m_bMark;
		pHalfedge = cell1.pHalfedge;
		m_nFaceLevel = cell1.m_nFaceLevel;
		IsSubdivid = cell1.IsSubdivid;

		 return *this;
     }
	~CFace(){
	/*SAFEDELETE(pHalfEdge);*/
	};
	
public://operation	
	bool IsBoundaryCell();  
	int  BoundaryCellType();
	char IsBoundaryCell_LRUD();
};


class CEdge{
public:
	int          m_nTwoVer[2];
	double       m_dLength;
    CHalfEdge    *pLeftHalfEdge;
	CHalfEdge    *pRightHalfEdge;
public:     
	 CEdge()
	{
       m_nTwoVer[0] = -1;
	   m_nTwoVer[1] = -1;
	   m_dLength = 0;
	   pLeftHalfEdge = NULL;
	   pRightHalfEdge = NULL;
	}
    
	CEdge & operator=(const CEdge &e1)
	{
		m_nTwoVer[0] = e1.m_nTwoVer[0];
		m_nTwoVer[1] = e1.m_nTwoVer[1];
		m_dLength    = e1.m_dLength;
		pLeftHalfEdge  = e1.pLeftHalfEdge;
		pRightHalfEdge = e1.pRightHalfEdge;
		return *this;
     }
	
	~CEdge(){
	
		/*SAFEDELETE(pLeftHalfEdge);
	SAFEDELETE(pRightHalfEdge);*/
	};
	
public:  //operation
	char edge_type();  
	bool IsBoundaryEdge(double x1, double x2, double y1, double y2);
	char EdgeTypeLRUD(CFace lyonface);  
	
};

bool IsSameEdge(CEdge e1, CEdge e2); 


class CHalfEdge{
public:
	CHalfEdge   *pPrevHalfedge;
	CHalfEdge   *pNextHalfedge;
	CEdge       *pEdge;
	CFace       *pFace;
	CVertex     *pVertex;
public:
	CHalfEdge()
	{
		pPrevHalfedge = NULL;
		pNextHalfedge = NULL;
		pEdge = NULL;
		pFace = NULL;
		pVertex = NULL;

	}
	~CHalfEdge()
	{
		SAFEDELETE(pPrevHalfedge);
		SAFEDELETE(pNextHalfedge);
		SAFEDELETE(pEdge);
		SAFEDELETE(pFace);
		SAFEDELETE(pVertex);
	}

public:   //operation
	inline CVertex *source() const
	{	
		return this->pPrevHalfedge->pVertex;
	}
	inline CVertex *target() const
	{	
		return this->pVertex;
	}
	CHalfEdge *dual() const
	{
			CEdge *pe=this->pEdge;
			assert(pe!=NULL);
			CHalfEdge *ph;
			if(pe->pLeftHalfEdge != this)
		    	ph = pe->pLeftHalfEdge;
			else
				ph = pe->pRightHalfEdge;
			return ph;
	}
   CHalfEdge& operator=(const CHalfEdge &hedge)
	{
		pPrevHalfedge=hedge.pPrevHalfedge;
		pNextHalfedge=hedge.pNextHalfedge;
		pEdge=hedge.pEdge;
		pFace=hedge.pFace;
		pVertex=hedge.pVertex;
      return *this;
	}

};

//define quad_meshes with half-edge structure;
class QuadMesh
{
public:
	QuadMesh(void);
	~QuadMesh(void);	

public:
	vector<point_t> mPoints;                                                    //点集
    vector<quad_t> mQuads;                                                      //四边形网格集合
	vector<edge_t> mEdges;      
    //typedef map<pair<index_t, index_t>, index_t > directedEdge2indexMap_t;
    //directedEdge2indexMap_t mDirectedEdge2heIndex;                             //根据有向边获取半边索引


	int              m_nMeshLevel;
	vector<CVertex>  m_Vertices;
	vector<CEdge>	 m_Edges;
	vector<CFace>    m_Faces;	

public:
	void ReadFile(ifstream& is);												//read data, maintain mPoints, mQuads, mEdgs; 
	void GeneratMesh(void);														//Todo

};


