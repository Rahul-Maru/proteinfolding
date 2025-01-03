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
using std::endl;
using std::flush ;
using std::ios ;
#include <iomanip>
using std::setprecision ;
using std::setw ;
#include <fstream>
using std::ifstream ;
using std::ofstream ;
#include <cmath>
#include "superpose_weightedRMS.h"
#include "jacobi.h"
#include "macros.h"
#define sup_verbose 0
#define sup_gibberish 0

//double CM_A[3] , CM_B[3] ; // Centers of Mass of vector sets A and B
//double **DELTA_Coords_PLUS , **DELTA_Coords_MINUS ;
//double Quat[4][4] ;
/* *******************************************
 * GLOBAL FUNCTION DEFINITION                *
 * ******************************************* */
void SUPERPOSE_WRMS(float **CoordsA , float **CoordsB , float *weights , int Size , float *RMSD , float ROTATION_MAT[][3] , double CM_A[] , double CM_B[] )
{
	//cout << "Performing Quaternion Superposition..." << std::flush ;
	//double CM_A[3] ,  CM_B[3] ;
	double **DELTA_Coords_PLUS ;
	double **DELTA_Coords_MINUS ;
	double Quat[4][4] ;
	float eigen_values[4] ;
	float eigen_vectors[4][4] ;
	float local_rmsd = -1 ;
	/* LOCAL FUNCTION PROTOTYPES */
	void Calculate_Center_of_Mass_WRMS( float ** , float ** , int , double [] , double []) ; 
	void Calculate_Coordinate_differentials_WRMS( float ** , float ** , int , double [] , double [] , double *** , double ***) ; 
	void Calculate_Quaternion_Matrix_WRMS( double ** , double ** , float * , int , double[][4] ) ;
	void Calculate_RMSD_WRMS(float [4] , int , float* ) ;
	void Calculate_Best_Rotation_Matrix_WRMS( float [][4] , float [][3]) ;
	void DEALLOC_WRMS(double ** , double ** , int ) ;
	
	Calculate_Center_of_Mass_WRMS(CoordsA , CoordsB , Size , CM_A , CM_B ); // for each of the two vector sets
	Calculate_Coordinate_differentials_WRMS(CoordsA , CoordsB , Size , CM_A , CM_B , &DELTA_Coords_PLUS , &DELTA_Coords_MINUS) ; 
	Calculate_Quaternion_Matrix_WRMS( DELTA_Coords_PLUS , DELTA_Coords_MINUS , weights,  Size ,  Quat ) ; // contructs a quaterion from the coordinate differences
	JACOBI_ROTATIONS( Quat , eigen_values , eigen_vectors ) ; // performs jacobi rotations on quaternion;  
								  // outputs sorted Eigen values & vectors
	Calculate_RMSD_WRMS( eigen_values , Size , &local_rmsd  ) ;
	*RMSD = local_rmsd ;
	
	Calculate_Best_Rotation_Matrix_WRMS( eigen_vectors , ROTATION_MAT ) ;
	DEALLOC_WRMS(DELTA_Coords_PLUS , DELTA_Coords_MINUS , Size ) ;

}
void Calculate_Center_of_Mass_WRMS(float **CoordsA , float **CoordsB , int Size , double CM_A[] , double CM_B[] )
{
	
	CM_A[ 0 ] = CM_A[ 1 ] = CM_A[ 2 ] = 0.0 ;
	CM_B[ 0 ] = CM_B[ 1 ] = CM_B[ 2 ] = 0.0 ;
	double temp = 0 ;
	for( int i = 0 ; i < 3 ; i++ )
	{
		temp = 0 ;
		for( int j = 0 ; j < Size ; j++ )
			temp+=CoordsA[j][i] ;
		CM_A[ i ] = temp / Size ;
	}
	for( int i = 0 ; i < 3 ; i++ )
	{
		temp = 0 ;
		for( int j = 0 ; j < Size ; j++ )
			temp+=CoordsB[j][i] ;
		CM_B[ i ] = temp / Size ;
	}
	if( sup_verbose )
	{
		cout << setprecision(3) << ios::fixed;
		cout << "\n\tCenter of Mass of vectors set A : " ;
		for( int i = 0 ; i < 3 ; i++ )
			cout << std::setw(9) << CM_A[ i ] ;
		cout << endl ;
		cout << "\tCenter of Mass of vectors set B : " ;
		for( int i = 0 ; i < 3 ; i++ )
			cout << std::setw(9) << CM_B[ i ] ;
		cout << endl ;
		cout << "\tTarget vector       : " ;
		for( int i = 0 ; i < 3 ; i++ )
			cout << std::setw(9) << CM_B[i] - CM_A[ i ] ;
		cout << endl ;
	}
}

void Calculate_Coordinate_differentials_WRMS(float **CoordsA , float **CoordsB , int Size , double CM_A[] , double CM_B[] , double ***DELTA_Coords_PLUS , double ***DELTA_Coords_MINUS ) 
{
	//Allocating space for DELTA_X_PLUS and DELTA_X_MINUS.
	(*DELTA_Coords_PLUS) = new double* [Size] ;
	for( int i = 0 ; i < Size ; i ++ )
		(*DELTA_Coords_PLUS)[i] = new double [3] ;
	(*DELTA_Coords_MINUS) = new double* [Size] ;
	for( int i = 0 ; i < Size ; i ++ )
		(*DELTA_Coords_MINUS)[i] = new double [3] ;
	
	// Calculating coordinate differences, DELTA_X_PLUS and DELTA_X_MINUS. Assume A to be Movable and B fixed
	for( int i = 0 ; i < Size  ; i++ )
		for( int j = 0 ; j < 3 ; j ++ )
		{
			(*DELTA_Coords_PLUS)[i][j] =  ( CoordsA[i][j] - CM_A[j] ) + ( CoordsB[i][j] - CM_B[j] ) ;
			(*DELTA_Coords_MINUS)[i][j] = ( CoordsA[i][j] - CM_A[j] ) - ( CoordsB[i][j] - CM_B[j] ) ;
		}
}

void Calculate_Quaternion_Matrix_WRMS( double **DELTA_Coords_PLUS , double **DELTA_Coords_MINUS , float *weights , int Size, double Quat[][4] )
{
	//Calculate Square-Symmetric Quaternion Matrix
	//initialize
	for( int i = 0 ; i < 4 ; i++ )
		for( int j = 0 ; j < 4 ; j++ )
			Quat[i][j] = 0.0 ;
	
	// Filling upper triangle of the quaternion matrix
	for( int i = 0 ; i < Size  ; i++ )
	{
		//Diags =Sum of squared cyclic coordinate differences
		Quat[0][0] += pow( weights[i] , 2 ) *
			      (pow( DELTA_Coords_MINUS[i][0] , 2 ) +
			       pow( DELTA_Coords_MINUS[i][1] , 2 ) +
			       pow( DELTA_Coords_MINUS[i][2] , 2 )) ;
		
		Quat[1][1] += pow( weights[i] , 2 ) *
			     (pow( DELTA_Coords_PLUS[i][1]  , 2 ) +
			      pow( DELTA_Coords_PLUS[i][2]  , 2 ) +
			      pow( DELTA_Coords_MINUS[i][0] , 2 )) ;
		
		Quat[2][2] += pow( weights[i] , 2 ) *
			     (pow( DELTA_Coords_PLUS[i][0]  , 2 ) +
			      pow( DELTA_Coords_PLUS[i][2]  , 2 ) +
			      pow( DELTA_Coords_MINUS[i][1] , 2 )) ;
		
		Quat[3][3] += pow( weights[i] , 2 ) *
			     (pow( DELTA_Coords_PLUS[i][0]  , 2 ) +
			      pow( DELTA_Coords_PLUS[i][1]  , 2 ) +
			      pow( DELTA_Coords_MINUS[i][2] , 2 )) ;
		// Cross differences
		Quat[0][1] += pow( weights[i] , 2 ) *
			     (DELTA_Coords_PLUS[i][1] * DELTA_Coords_MINUS[i][2] -
			      DELTA_Coords_MINUS[i][1] * DELTA_Coords_PLUS[i][2]) ;

		Quat[0][2] += pow( weights[i] , 2 ) *
			     (DELTA_Coords_MINUS[i][0] * DELTA_Coords_PLUS[i][2] -
			      DELTA_Coords_PLUS[i][0] * DELTA_Coords_MINUS[i][2]) ;

		Quat[0][3] += pow( weights[i] , 2 ) *
			     (DELTA_Coords_PLUS[i][0] * DELTA_Coords_MINUS[i][1] -
			      DELTA_Coords_MINUS[i][0] * DELTA_Coords_PLUS[i][1]) ;

		Quat[1][2] += pow( weights[i] , 2 ) *
			     (DELTA_Coords_MINUS[i][0] * DELTA_Coords_MINUS[i][1] -
			      DELTA_Coords_PLUS[i][0] * DELTA_Coords_PLUS[i][1]) ;

		Quat[1][3] += pow( weights[i] , 2 ) *
			     (DELTA_Coords_MINUS[i][0] * DELTA_Coords_MINUS[i][2] -
			      DELTA_Coords_PLUS[i][0] * DELTA_Coords_PLUS[i][2]) ;

		Quat[2][3] += pow( weights[i] , 2 ) *
			     (DELTA_Coords_MINUS[i][1] * DELTA_Coords_MINUS[i][2] -
			      DELTA_Coords_PLUS[i][1] * DELTA_Coords_PLUS[i][2]) ;
	}
	// Fill the rest by transposing it onto itself
	Quat[1][0] = Quat[0][1] ; Quat[2][0] = Quat[0][2] ; Quat[2][1] = Quat[1][2] ;  
	Quat[3][0] = Quat[0][3] ; Quat[3][1] = Quat[1][3] ; Quat[3][2] = Quat[2][3] ;
	if( sup_verbose )
	{
		cout << "\tQuaternion_Matrix:\n" ;
		for( int i = 0 ; i < 4 ; i++ )
		{
			cout << "\t" ;
			for( int j = 0 ; j < 4 ; j++ )
				cout << setw(12) << Quat[i][j] ;
			cout << endl;
		}
		
	}
}
void Calculate_RMSD_WRMS( float eigen_values[] , int Size , float *RMSD ) 
{
	*RMSD = sqrt(fabs( eigen_values[3] / Size )) ;
	
}
void Calculate_Best_Rotation_Matrix_WRMS( float eigen_vectors[][4] , float ROTATION_MAT[][3]) 
{
	ROTATION_MAT[0][0] = pow( eigen_vectors[0][3] , 2 ) + pow( eigen_vectors[1][3] , 2 ) -  
			     pow( eigen_vectors[2][3] , 2 ) - pow( eigen_vectors[3][3] , 2 ) ;
	ROTATION_MAT[1][0] = 2 *( 
			          eigen_vectors[1][3] * eigen_vectors[2][3] +
			          eigen_vectors[0][3] * eigen_vectors[3][3] 
			        ) ;
	ROTATION_MAT[2][0] = 2 *( 
			          eigen_vectors[1][3] * eigen_vectors[3][3] -
			          eigen_vectors[0][3] * eigen_vectors[2][3] 
			        ) ;
	ROTATION_MAT[0][1] = 2 *( 
			          eigen_vectors[1][3] * eigen_vectors[2][3] -
			          eigen_vectors[0][3] * eigen_vectors[3][3] 
			        ) ;
	ROTATION_MAT[1][1] = pow( eigen_vectors[0][3] , 2 ) + pow( eigen_vectors[2][3] , 2 ) -  
			     pow( eigen_vectors[1][3] , 2 ) - pow( eigen_vectors[3][3] , 2 ) ;
	ROTATION_MAT[2][1] = 2 *( 
			          eigen_vectors[2][3] * eigen_vectors[3][3] +
			          eigen_vectors[0][3] * eigen_vectors[1][3] 
			        ) ;
	ROTATION_MAT[0][2] = 2 *( 
			          eigen_vectors[1][3] * eigen_vectors[3][3] +
			          eigen_vectors[0][3] * eigen_vectors[2][3] 
			        ) ;
	ROTATION_MAT[1][2] = 2 *( 
			          eigen_vectors[2][3] * eigen_vectors[3][3] -
			          eigen_vectors[0][3] * eigen_vectors[1][3] 
			        ) ;
	ROTATION_MAT[2][2] = pow( eigen_vectors[0][3] , 2 ) + pow( eigen_vectors[3][3] , 2 ) -  
			     pow( eigen_vectors[1][3] , 2 ) - pow( eigen_vectors[2][3] , 2 ) ;
	if(gibberish)
	{
		cout << "\n\t Rotation Matrix corresponding to the best superposition :\n" ;
		for( int i = 0 ; i < 3 ; i ++ )
		{
			cout << "\t" ;
			for( int j = 0 ; j < 3 ; j ++ )
				cout << setw(12) <<ROTATION_MAT[i][j] ;
			cout << endl ;
		}
	}
}

void DEALLOC_WRMS( double **a , double **b , int Size )
{
	for ( int i = 0 ; i < Size ; i++ )
		delete[]  a[i]  ;
	delete[] a ;
	for ( int i = 0 ; i < Size ; i++ )
		delete[]  b[i]  ;
	delete[] b ;
}
