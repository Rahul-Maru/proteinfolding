/*******************************************************************************
 * Copyright (c) 2005, Arun S Konagurthu, Monash University.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification are permitted provided that the following conditions are met:
 * 
 * (1) Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * (2) Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * (3) Neither the name of the University of Melbourne nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *  
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESSINTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY,OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/
#include <iostream>
using std::cout ;
using std::cerr ;
using std::endl ;
#include <iomanip>
#include <cmath>
#include "jacobi.h"
//#include "globals.h"

/* ***************************
 * GLOBAL FUNCTION DEFINITION
 * *************************** */
#define JACOBI_SIZE 4
#define ROTATE(Quat,i,j,k,l) g=Quat[i][j];h=Quat[k][l];Quat[i][j]=g-s*(h+g*tau);\
	Quat[k][l]=h+s*(g-h*tau);


//extern double Quat[][4] ;
//float Quat[3][3] ;
#define jac_verbose 0
void eigsrt( float d[] , float v[][4] , int n )
{
	int k,j,i;
	float p;

	for (i=0;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=0;j<n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}
void JACOBI_ROTATIONS( double Quat[][4] , float eigen_values[4] , float eigen_vectors[][4])
{
	// This program performs Jacobi rotations for calculation 
	// of eigenvalues and eigen vectors of a real symmetric 
	// quarternion matrix "quat[4][4]"
	
	//float d[4] , **v ;
	float d[4] , v[4][4] ;
	int  iq , ip ;
	float tresh , theta , tau , t , sm , s , h , g , c , *b , *z ;
	int NROT ;

	// Allocating space for b-single dimension
	b = new float [JACOBI_SIZE] ;
	
	// Allocating space for z-single dimension
	z = new float [JACOBI_SIZE] ;

	// Allocating space for v- double dimension
	//v = new (float *)[JACOBI_SIZE] ;
	//for( int i = 0 ; i < JACOBI_SIZE ; i++ )
	//	v[i] = new float [JACOBI_SIZE] ;
	
		
	
	// Initialize v to Identity matrix 
	for ( ip = 0 ; ip < JACOBI_SIZE ; ip++ )
	{
		for( iq = 0 ; iq < JACOBI_SIZE ; iq++ )
			v[ip][iq] = 0.0 ;
		v[ip][ip] = 1.0 ;		
	}
	
	
	// Initialize b and d to diagonal of a  
	for ( ip = 0 ; ip < JACOBI_SIZE ; ip++ )
	{
		b[ip] = d[ip] = Quat[ip][ip] ;
		z[ip] = 0.0 ;
	}
	NROT = 0 ;
	for( int i = 0 ; i <= 50 ; i++ )
	{
		sm = 0.0 ;
		// calculating sum of off-diagonal elements
		
		for ( ip = 0 ; ip < JACOBI_SIZE - 1 ; ip++ )
			for( iq = ip + 1 ; iq < JACOBI_SIZE ; iq++ )
				sm += fabs( Quat[ ip ][ iq ] ) ;
		if( sm == 0.0 )
		{
			delete[] z ;
			delete[] b ;
			for( int i = 0 ; i < JACOBI_SIZE ; i++ )
				eigen_values[i]= d[i] ;
			
			for ( int i = 0 ; i < JACOBI_SIZE ; i++ )
				for( int j = 0 ; j < JACOBI_SIZE ; j++ )
					eigen_vectors[i][j]=v[i][j] ;
			
			if(jac_verbose)
			{
				cout << std::setprecision(4) ;
				cout << "\t Eigenvalues of the Quaternion matrix:" ;
				for( int i = 0 ; i < JACOBI_SIZE ; i++ )
					cout << std::setw(17) << d[i] ;
				cout << endl ;
				cout << "\t Eigenvectors of the Quaternion Matrix:\n" ;
				for ( int i = 0 ; i < JACOBI_SIZE ; i++ )
				{
					cout<<"\t" ;
					for( int j = 0 ; j < JACOBI_SIZE ; j++ )
						cout << std::setw(12) << eigen_vectors[i][j] ;
					cout << "\n" ;
				}
						
				//for( int i = 0 ; i < JACOBI_SIZE ; i++ )
				//	cout << std::setw(7) << eigen_values[i] ;
				//cout << endl ;
				cout << "\tNumber of jacobi rotations required to converge: " ;
				cout << NROT << endl ;
			}
		 	eigsrt( eigen_values , eigen_vectors , JACOBI_SIZE ) ;
			if(jac_verbose)
			{
				cout << "\tSorted Eigenvalues of the Quaternion matrix:" ;
				for( int i = 0 ; i < JACOBI_SIZE ; i++ )
					cout << std::setw(17) << eigen_values[i] ;
				cout << endl ;
				cout << "\tSorted Eigenvectors of the Quaternion Matrix:\n" ;
				for ( int i = 0 ; i < JACOBI_SIZE ; i++ )
				{
					cout<<"\t" ;
					for( int j = 0 ; j < JACOBI_SIZE ; j++ )
						cout << std::setw(12) << eigen_vectors[i][j] ;
					cout << "\n" ;
				}
			}
			return ;
				
		}
		if( i < 4 )
			tresh = 0.2 * sm / ( JACOBI_SIZE * JACOBI_SIZE ) ;
		else
			tresh = 0.0 ;
		
		for ( ip = 0 ; ip < JACOBI_SIZE -1 ; ip++ )
		{
			for( iq = ip + 1 ; iq < JACOBI_SIZE ; iq++ )
			{
				g = 100.0 * fabs( Quat[ip][iq] ) ;
				// After four sweeps skip rotation if the off-diagonal element is small
				if ( i > 4 && (float)( fabs( d[ip] ) + g ) ==  (float)fabs( d[ip] ) &&
				     (float)( fabs( d[iq] ) + g ) == (float)fabs( d[iq] ) 
				   )
					Quat[ip][iq] = 0.0 ;
				else if ( fabs( Quat[ip][iq] ) > tresh )
				{
					h = d[iq] - d[ip] ;
					if( (float)( fabs(h) + g) == (float)fabs(h) )
					{
						t = ( Quat[ip][iq])/h ;      // t = 1/2theta
					}
					else
					{
						theta = 0.5 * h/(Quat[ip][iq]) ;
						t =  1.0 / ( fabs(theta) + sqrt( 1.0 + theta*theta ) ) ;
						if( theta  < 0.0 ) 
							t = -t ;
					}
					c = 1.0 / sqrt( 1+ t*t ) ;
					s = t* c ;
					tau = s / ( 1.0 + c ) ;
					h = t * Quat[ip][iq] ;
					z[ ip ] -= h ;
					z[ iq ] += h ;
					d[ ip ] -= h ;
					d[ iq ] += h ;
					Quat[ ip ][ iq ] = 0.0 ;
					for( int j = 0 ; j <= ip - 1 ; j++  )
					{
						ROTATE(Quat,j,ip,j,iq) 
					}
					for( int j = ip + 1 ; j <= iq -1 ; j++ )
					{
						ROTATE(Quat,ip,j,j,iq) 
					}
					for( int j =  iq + 1 ; j < JACOBI_SIZE ; j++ )
					{
						ROTATE(Quat,ip,j,iq,j) 
					}
					for( int j =  0 ; j < JACOBI_SIZE ; j++ )
					{
						ROTATE(v,j,ip,j,iq) 
					}
					++NROT ;
				}
			}
		}

		for( ip = 0 ; ip < JACOBI_SIZE ; ip++ )
		{
			b[ ip ] += z[ ip ] ;
			d[ ip ] = b[ ip ] ;
			z[ ip ] = 0.0 ; 
		}
	}
	cerr  << "too many iterations in jacobi" << endl ;
	
}
/*
int main()
{
	Quat[0][0] = 5 ; Quat[0][1] = -6 ; Quat[0][2] = -6 ;
	Quat[1][0] = -1 ; Quat[1][1] = 4 ; Quat[1][2] = 2 ;
	Quat[2][0] = 3 ; Quat[2][1] = -6 ; Quat[2][2] = -4 ;
	JACOBI_ROTATIONS();
	
	
}
*/
#undef JACOBI_SIZE
