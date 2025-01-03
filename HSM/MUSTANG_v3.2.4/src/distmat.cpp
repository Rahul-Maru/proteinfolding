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
using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::flush;
#include<iomanip>
using std::setprecision;
using std::setw;

#include <fstream>
using std::ifstream;

#include<math.h>
#include "macros.h"
#include "globals.h"
#include "distmat.h"
/* ************************************************
 * GLOBAL FUNCTION DEFINITION  
 * ************************************************/
void CALCULATE_DISTANCE_MATRICES( )
{
	if(!meditate)
		cout << "Calculating Distance Matrices..." << flush ;
	
	void ALLOCATE_Distance_Matrices( ) ;
	ALLOCATE_Distance_Matrices() ;

	float Euclidean_Distance( float [] , float [] ) ; // local function prototype
	// Calculating Intra CA-CA distances for every structure
	for( int K = 0 ; K < NSTRUCTS ; K++ )
	{
		for( int i = 0 ; i < PROT_SIZES[K] ; i++ )
			for( int j = 0 ; j <= i ; j++ )
			{
				distance_matrices[K][i][j] = Euclidean_Distance( PROT[K][i].CA_coords , PROT[K][j].CA_coords ) ;
				if( i != j )
					distance_matrices[K][j][i] = distance_matrices[K][i][j] ;
			}
	}
	
	if(gibberish)
	{
		//cout << std::fixed ;
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::showpoint);

		for( int K = 0 ; K < NSTRUCTS ; K++ )
		{
			cout <<endl <<endl ;
			for( int i = 0 ; i < PROT_SIZES[K] ; i++ )
			{
				for( int j = 0 ; j < PROT_SIZES[K] ; j++ )
					cout << setw(6) << setprecision(1)  << distance_matrices[ K ][ i ][ j ] << " " ;
				cout << endl ;
			}
			cout <<endl <<endl ;
		}
		//exit(0);
	}
	if(!meditate)
	{
		cout << setw(37) << " " ;
		cout << "[ \033[32;4mOK\033[0m ]\n" ;
	}

}
/* *****************************************
   * LOCAL FUNCTION DEFINITIONS
   ****************************************/
float Euclidean_Distance( float A[] , float B[] )
{
	return ( sqrt( pow( (B[0] - A[0]) , 2) + pow( (B[1] - A[1]) , 2) + pow( (B[2] - A[2]) , 2)  ) ) ;
}	


void ALLOCATE_Distance_Matrices( ) 
{
	distance_matrices = new float** [NSTRUCTS] ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		distance_matrices[i] = new float* [ PROT_SIZES[i] ] ;
		for( int j = 0 ; j < PROT_SIZES[i] ; j++ )
			distance_matrices[i][j] = new float [ PROT_SIZES[i] ] ;
	}
}
