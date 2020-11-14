
#ifndef _NUMERICAL_H_
#define _NUMERICAL_H_

// #include "stdafx.h"
#include <iostream>
#include <cmath>
// #include <vector>


#define PI 3.14159265358979
using namespace std;

// #pragma comment(lib,"CLAPACK.lib")
// #include"CLAPACK/f2c.h"
// extern "C"
// { 
// #include "CLAPACK/clapack.h" 
// }
// using namespace std;


// void SolveLinearEquation(const vector<vector<double>> &coeff_matrix,const vector<double> &load_vect,
// 						 vector<double> &m_solute)
// {
// 	int matrix_size = coeff_matrix.size();

// 	double *m_CoMatrix=new double[matrix_size*matrix_size];
// 	double *m_rhs=new double[matrix_size];
// 	int M = matrix_size;
// 	int nrhs = 1;
// 	int lda,ldb;
// 	int *ipiv=new int[matrix_size];
// 	int INFO;
// 	lda = M;
// 	ldb = M;

// 	for(int i=0;i<matrix_size;i++)
// 	{
// 		for(int j=0;j<matrix_size;j++)
// 		{	
// 			m_CoMatrix[j*matrix_size+i] = coeff_matrix[i][j];
// 		}	
               
// 	}
//     for ( int i=0;i<matrix_size;i++)
// 	    m_rhs[i] = load_vect[i];      //方程右端项
	
// 	//求解线性方程组
// 	dgesv_((integer *)&M, (integer *)&nrhs, m_CoMatrix, (integer *)&lda, (integer *)ipiv, m_rhs,(integer *)&ldb,(integer *) &INFO);	
// 	if (0 == INFO)  //为0表示求解正确，m_rhs存储的是求出的解
// 	{
// 		m_solute.resize(matrix_size);
// 		for (int i=0;i<matrix_size;i++)
// 		{
// 			m_solute[i] = m_rhs[i];
// 			/*cout<<m_solute[i]<<endl;*/
// 		}
// 		/*cout<<'\n';*/
// 	}
	
// 	delete []ipiv;
// 	delete []m_CoMatrix;
// 	delete []m_rhs;

// }

// void SolveOverEquation(int m,int n,const vector<vector<double>> &coeff_matrix,const vector<double> &rhs,
// 						 vector<double> &m_solute)
// {
// 	    int nrhs =1;
// 		int lda=m;   //max(1,m)
// 		int ldb=m;   //max(1,m,n) 
// 		int info, lwork, rank;
//         double rcond = -1.0;/* Negative rcond means using default (machine precision) value */
//         double wkopt;
//         double* work;
//         /* Local arrays */
//         /* iwork dimension should be at least 3*min(m,n)*nlvl + 11*min(m,n)=3*n*n1v1+11*n,
//                 where nlvl = max( 0, int( log_2( min(m,n)/(smlsiz+1) ) )+1 )
//                 and smlsiz = 25 */
//          int n1v1=(int)n/26;
// 		 int *iwork=new int[3*n*n1v1+11*n]; //int iwork[3*n*1+11*n];
//          double *s=new double[n];
//          double *a=new double[m*n];
// 		 for(int j=0;j<n;j++)
// 		 {
// 			 for(int i=0;i<m;i++)
// 				 a[j*m+i]=coeff_matrix[i][j];
// 		 }
// 		  /*ofstream coeff_file;
// 	     coeff_file.open("coeffA.txt");
// 		 for(int i=0;i<m;i++)
// 		 { 
// 			  for(int j=0;j<n;j++)
// 			     coeff_file<<a[j*m+i]<<" ";
// 			 coeff_file<<endl;
// 		 }
// 		 coeff_file.close();*/
// 		double *b=new double[m];
// 		for(int i=0;i<m;i++)
// 		    b[i]=rhs[i]; 
       
//         /* Query and allocate the optimal workspace */
//         lwork = -1;
// 		dgelsd_( (integer *)&m, (integer *)&n,(integer *)&nrhs, a,(integer *)&lda,
// 			     b, (integer *)&ldb, s, &rcond, (integer *)&rank, &wkopt,(integer *)&lwork,(integer *)iwork, (integer *)&info );
			
                        
//         //lwork = (int)wkopt;
//         //work = (double*)malloc( lwork*sizeof(double) );
//         ///* Solve the equations A*X = B */
//         //dgelsd_( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info );
                       
//         /* Check for convergence */
// 		if(info==0)
// 		{
// 			m_solute.clear();
// 			 for( int i = 0; i< n; i++ ) 
// 			 {
// 				 m_solute.push_back(b[i]); //cout<<coeff_basis[i]<<" ";
//              }
// 		}
// 		else 
// 		{
//                 printf( "The algorithm computing SVD failed to converge;\n" );
//                 printf( "the least squares solution could not be computed.\n" );
//                 exit( 1 );
//         }
//        delete[] iwork;
// 	   delete[] b;
// 	   delete[] s;
// 	   delete[] a;

// }

// void GaussIntegalPoint(int n, vector<double> &w,vector<double> &gausspoint)
// {
// 	int iter,m,el;
// 	int mm;
// 	double nn;
// 	for(int i = 0; i < n; i++)
// 	{
// 		w.push_back(0);	gausspoint.push_back(0);
// 	}
// 	iter = 2;
// 	m = (n + 1) / 2; el = n * (n + 1);
// 	mm = 4 * m -1;
// 	nn = (1 - (1 - 1.0 / n) / (8 * n * n));
// 	vector<double> t,x0;
// 	double tt;
// 	for(int i = 3; i <= mm; i+=4)
// 	{
// 		tt = (PI / (4 * n + 2)) * i;
// 		t.push_back(tt);
// 		x0.push_back(nn*cos(tt));
// 	}
// 	vector<double> pkm1,pkp1,t1,pk;
// 	int x0_size = x0.size();
// 	vector<double> den,d1,dpn,d2pn,d3pn,d4pn,u,v,h,p,dp,fx;
// 	t1.resize(x0_size);  pkm1.resize(x0_size);
// 	pk.resize(x0_size);  t1.resize(x0_size);
// 	pkp1.resize(x0_size);den.resize(x0_size);
// 	dpn.resize(x0_size); d2pn.resize(x0_size);
// 	d3pn.resize(x0_size);d4pn.resize(x0_size);
// 	d1.resize(x0_size);  u.resize(x0_size);
// 	v.resize(x0_size);   h.resize(x0_size);
// 	p.resize(x0_size);   fx.resize(x0_size);
// 	dp.resize(x0_size);

// 	for(int i = 0; i < iter; i++)
// 	{
// 		for(int j = 0; j < x0_size; j++)
// 		{
// 			pkm1[j] = 1;
// 			pk[j] = x0[j];
// 		}
// 		for(int k = 2; k <= n; k++)
// 			for(int r = 0; r < x0_size; r++)
// 			{
// 				t1[r] = x0[r] * pk[r];
// 				pkp1[r] = t1[r] - pkm1[r] - 1.0 * (t1[r] - pkm1[r]) / k + t1[r];
// 				pkm1[r] = pk[r]; pk[r] = pkp1[r];
// 			}
// 			for(int j = 0; j < x0_size; j++)
// 			{
// 				den[j] = 1.0 - x0[j] * x0[j];
// 				d1[j] = n * (pkm1[j] - x0[j] * pk[j]);
// 				dpn[j] = d1[j] / den[j];
// 				d2pn[j] = (2.0 * x0[j] * dpn[j] - el * pk[j]) / den[j];
// 				d3pn[j] = (4 * x0[j] * d2pn[j] + (2-el) * dpn[j]) / den[j];
// 				d4pn[j] = (6 * x0[j] * d3pn[j] + (6-el) * d2pn[j]) / den[j];
// 				u[j] = pk[j] / dpn[j]; v[j] = d2pn[j] / dpn[j];
// 				h[j] = -u[j] * (1 + 0.5 * u[j] * (v[j] + u[j] * (v[j] * v[j] - u[j] * d3pn[j] / (3 * dpn[j]))));
// 				p[j] = pk[j] + h[j] * (dpn[j] + 0.5 * h[j] * (d2pn[j] + 1.0 * h[j] / 3 * (d3pn[j] + 0.25 * h[j] * d4pn[j])));
// 				dp[j] = dpn[j] + h[j] * (d2pn[j] + 0.5 * h[j] * (d3pn[j] + h[j] * d4pn[j] / 3));
// 				h[j] = h[j] - p[j] / dp[j]; x0[j] = x0[j] + h[j];
// 			}
// 	}

// 		for(int j = 0; j < x0_size; j++)
// 		{
// 			gausspoint[j] = - x0[j] - h[j];
// 			fx[j] = d1[j] - h[j] * el * (pk[j] + h[j]/2 * (dpn[j] + h[j]/3 * (d2pn[j] + h[j] / 4 * (d3pn[j] + 0.2 * h[j] * d4pn[j]))));
// 			w[j] = 2 * (1.0-gausspoint[j] * gausspoint[j]) / (fx[j] * fx[j]);
// 		}
// 		if(m + m > n) gausspoint[m-1] = 0;
// 		if((m + m) != n) m = m - 1;
// 		vector<int> jj,n1j;
// 		jj.resize(m); n1j.resize(m);
// 		for(int j = 0; j < m; j++)
// 		{
// 			jj[j] = j; n1j[j] = n - jj[j] - 1;
// 			gausspoint[n1j[j]] = -gausspoint[jj[j]];
// 			w[n1j[j]] = w[jj[j]];
// 		}

// }

// int VectorMax(vector<double> Array)//冒泡排序
// {
// 	int maxindex = 0;
// 	int i =0;
//     int n = (int)Array.size();

// 	for(i=1;i<n;i++)
// 	{
// 		if(Array[i]>Array[maxindex])
// 			maxindex=i;
// 	}
// 	return maxindex;
// }

// int partition(vector<int> &RId,vector<double> &R,int low,int high)
// {
// 	double pivotkey;
// 	pivotkey=R[low];
// 	while(low<high)
// 	{
// 		while(low<high && R[high]>pivotkey)
// 			--high;
// 		//if(low<high)
// 		R[low]=R[high];
// 		while(low<high&&R[low]<pivotkey)
// 			++low;
// 		//if(low<high)
// 		R[high]=R[low];
		
// 	}
// 	R[low]=pivotkey;
// 	return low;
// }
/*
void Qsort(vector<int> &RId,vector<double> &R,int s,int t)
{
	int pivot;
	if(s<t)
	{
		pivot=partition(R,s,t); cout<<pivot<<endl;
		Qsort(R,s,pivot-1);
		Qsort(R,pivot+1,t);
	}

}
*/

#endif