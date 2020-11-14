#pragma once

#include <vector>
//#include "Numerical.h"

#define PI 3.14159265358979
class CPoint3D;

//ͨ�õ���ѧ��Ķ���

//mesh types; used for read data from obj files
typedef long index_t;

struct point_t {                    //����
	index_t index;                  //������
	double x, y, z;                 //��������
};

struct edge_t {                     //���������еı�
	index_t v[2];                   //�洢����������
	index_t& start() {
		return v[0];
	}

	index_t& end() {
		return v[1];
	}

	edge_t() {
		v[0] = v[1] = -1;           //��ʼ���������������
	}
};

struct triangle_t {                 //������
	index_t v[3];                   //����������������
	const index_t& i() const {
		return v[0];
	}
	const index_t& j() const {
		return v[1];
	}
	const index_t& k() const {
		return v[2];
	}
};

struct quad_t{						//�ı���
	index_t id;
	index_t v[4];
	const index_t& i() const {      //�ı����ĸ�������
		return v[0];
	}
	const index_t& j() const {
		return v[1];
	}
	const index_t& k() const {
		return v[2];
	}
	const index_t& l() const {
		return v[3];
	}
};


// Point3D and its operations
class CPoint3D
{
public:
	float	x, y, z;
	
	CPoint3D ( );					// normal constructor
	CPoint3D ( int, int, int );			// normal constructor
	CPoint3D ( float, float, float );		// normal constructor
	CPoint3D ( double, double, double );		// normal constructor
	virtual ~CPoint3D ();
	
	float	operator & ( CPoint3D& v );		// doc_product
	CPoint3D	operator + ( CPoint3D& v );
	CPoint3D	operator - ( CPoint3D& v );
	CPoint3D	operator * ( CPoint3D& v );		// cross_product
	CPoint3D	operator * ( int k );			// scale by int
	CPoint3D	operator * ( float k );			// scale by float
	CPoint3D	operator * ( double k );		// scale by double
	int operator ==(CPoint3D& v);
	//		int	operator = ( CPoint3D& v );
	void	unify();
	void	RangeUnify(float min,float max);
	float	length();
	//���Լ�ͶӰ���µ�����ϵ�м�
	CPoint3D ProjectToNewCoordinate(CPoint3D& coordX,CPoint3D& coordY,CPoint3D& coordZ);

};

void VectorToAngle(CPoint3D& vec,float& fPdir,float& fAdir);


class CAngle
{
	float m_fP;	//Polar angle
	float m_fA;	//azimuthal angle


	CPoint3D m_vPos;

	CAngle()
	{
		m_fP=0.0f;
		m_fA=0.0f;

		AngleToVector();
	}

	CAngle(float fP,float fA)
	{	
		m_fP=fP;
		m_fA=fA;
		
		AngleToVector();
	}

	~CAngle(){}

	CAngle operator +(CAngle& a)
	{
		return CAngle(m_fP+a.m_fP,m_fA+a.m_fA);
	}

	CAngle operator *(float a)
	{
		return CAngle(m_fP*a,m_fA*a);
	}


	CAngle operator -(CAngle& a)
	{
		return CAngle(m_fP-a.m_fP,m_fA-a.m_fA);
	}

	CAngle operator =(CAngle& a)
	{
		m_fP=a.m_fP;
		m_fA=a.m_fA;
		AngleToVector();
	}

	void CAngle::SetValue(float fP,float fA)
	{
		m_fP=fP;
		m_fA=fA;
		AngleToVector();
	}
	

	void AngleToVector()
	{		
		m_vPos.y=float(sin(PI/2-m_fP));
		m_vPos.x=float(cos(PI/2-m_fP)*cos(m_fA));
		m_vPos.z=float(-cos(PI/2-m_fP)*sin(m_fA));

		
	}


	

};

