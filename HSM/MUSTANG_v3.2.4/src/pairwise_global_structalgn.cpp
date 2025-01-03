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
using std::cin; 	using std::cout; 	using std::cerr; 
using std::flush; 	using std::endl; 	using std::ios;
using std::fixed;


#include<iomanip>
using std::setprecision ; using std::setw ; 

#include<fstream>
using std::ifstream ;

#include <ctime>
#include <cstdio>
#include <cstdlib>

#include "pairwise_global_structalgn.h"
#include "globals.h"
#include "macros.h"
#include "ew.h"
#include "refine_pairalgn.h"
#include "de_alloc_routines.h"
#define LINE_WIDTH 15 
template <class Type> void DE_ALLOC_2D( Type ** , int ) ;
void PAIRWISE_GLOBAL_STRUCTURAL_ALIGNMENT( int a , int b )
{
	//cout << std::fixed ;
	float **edge_weights ;
	//clock_t start, end ;
	//start = clock();
	CALCULATE_EDGE_WEIGHTS( a , b , &edge_weights ) ;
	//end = clock();
	//cout << "Exe_Time:" <<  (end -start) << endl ;
	//exit(0);
	void COPY_EDGE_WEIGHTS_TO_GLOBAL_DEFINITION( int , int , float ** ) ;
	//COPY_EDGE_WEIGHTS_TO_GLOBAL_DEFINITION( a , b,  edge_weights ) ;
	//if( gibberish )
	if( 0 )
	{
		//exit(0) ;
		cerr << "Pair: " << a << " " << b << endl ;
	
		cerr << "EdgeWeightsi_2: " << endl ;
		cerr.setf(ios::fixed, ios::floatfield);
		cerr.setf(ios::showpoint);

		cerr << setw(15) << " " ;
		for( int j = 0 ; j < PROT_SIZES[b] ; j++ )
			cerr << setw(10) << PROT[b][j].res_name << "(" << setw(3) << j << ")";
		cerr << "\n" ;
		for(int i = 0 ; i < PROT_SIZES[a] ; i++ )
		{
			cerr << setw(10) << PROT[a][i].res_name << "(" << setw(3) << i << ")" ;
			for( int j = 0 ; j < PROT_SIZES[b] ; j++ )
				cerr << setprecision(1) <<setw(15) << edge_weights[i][j] ;
			cerr << endl ;
		}
		exit(0);
	}

	void DYNAMIC_PROGRAMMING_METHOD_FOR_SEQUENTIAL_WEIGHTED_MATCHING(  int , int , float ** ) ;
	int NCYCLES = 1 ;
	for( int i = 0 ; i < NCYCLES ; i++ )
	{
		DYNAMIC_PROGRAMMING_METHOD_FOR_SEQUENTIAL_WEIGHTED_MATCHING( a , b , edge_weights  ) ;

		//DELETE_MISALIGNED_EDGE_WEIGHTS( a , b , edge_weights ) ;
		REFINE_PAIRWISE_ALIGNMENT( a , b , edge_weights ) ;
		COPY_EDGE_WEIGHTS_TO_GLOBAL_DEFINITION( a , b,  edge_weights ) ;
	}

	void DE_ALLOC_float2D( float ** , int ) ;
	DE_ALLOC_float2D( edge_weights , PROT_SIZES[a] ) ;
}
void  DYNAMIC_PROGRAMMING_METHOD_FOR_SEQUENTIAL_WEIGHTED_MATCHING( int a , int b , float **edge_weights )
{
	//if(!meditate)
	//	cout << "Performing Dynamic programming method....." << flush ;
	float arr[3] ;
	
	float **DP_MATRIX ;
	char **derivation;
	
	float max( float [] , char * ) ;
	
	// Allocate sufficient space for the DP matrix
	DP_MATRIX = new float* [ (PROT_SIZES[a]+1) ] ;
	derivation = new char* [ (PROT_SIZES[a]+1) ] ;

	for( int i = 0 ; i <  PROT_SIZES[a]+1 ; i ++ )
	{
		DP_MATRIX[i]  = new float [ ( PROT_SIZES[b] + 1 ) ] ;
		derivation[i] = new char  [ ( PROT_SIZES[b] + 1 ) ] ;
	}
	
	//boundary conditions
	DP_MATRIX[0][0]  = 0 ; 
	derivation[0][0] = '0' ;
	
	for( int x = 1 , y = 0 ; x < PROT_SIZES[b]+1 ; x++ )
	{
		DP_MATRIX[y][x] = 0 ;
		derivation[y][x] = '1' ;
	}
	
	for( int x = 0 , y = 1 ; y < PROT_SIZES[a]+1 ; y++ )
	{
		DP_MATRIX[y][x] = 0 ;
		derivation[y][x] = '2' ;
	}

	//rest of the matrix
	for( int y = 1 ; y < PROT_SIZES[a]+1 ; y++ )
		for( int x = 1 ; x < PROT_SIZES[b]+1 ; x++ )
		{
			//diagonal derivation
			arr[0] = DP_MATRIX[y-1][x-1] + edge_weights[y-1][x-1] ;
			
			//horizontal derivation
			arr[1] = DP_MATRIX[y][x-1]  ;
			
			//vertical
			arr[2] = DP_MATRIX[y-1][x]  ;
			
			DP_MATRIX[y][x] = max( arr , &derivation[y][x] ) ;
		}
	//cout << "out of rest-DP\n" ;
	//if( gibberish )
	if( 0 )
	{
		cout << "Pair: " << a << " " << b << endl ;
	
		cout << "DP_SCORES: " << endl ;
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::showpoint);

		cout << setw(10) << " " ;
		cout << setw(10) << "-" ;
		for( int j = 0 ; j < PROT_SIZES[b] ; j++ )
			cout << setw(10) << PROT[b][j].res_name ;
		cout << "\n" ;
		for(int i = 0 ; i < PROT_SIZES[a]+1 ; i++ )
		{
			if( i == 0 ) cout << setw(10) << "-" ;
			else	cout << setw(10) << PROT[a][i-1].res_name ;
			for( int j = 0 ; j < PROT_SIZES[b] ; j++ )
			{
				cout << setprecision(1) <<setw(7) << DP_MATRIX[i][j] ;
				cout << "(" ;
				if( derivation[i][j] == '0' ) cout<< "D)"  ;
				else if( derivation[i][j] == '1' ) cout<< "H)"  ;
				else if( derivation[i][j] == '2' ) cout<< "V)" ;
			}
			cout << endl ;
		}
	}
	
	//find size of alignment
	int y = PROT_SIZES[a] ;
	int x = PROT_SIZES[b] ;
	int algn_size = 0 ;
	//cout << y << " " << x << "\n" ;
	while( y >= 0 && x >=0 )
	{
		if( y == 0 && x == 0 )
			break;
		if( derivation[y][x] == '0' )
		{
			y--; x-- ;
		}
		else if( derivation[y][x] == '1' )
		{
			 x-- ;
		}
		else if( derivation[y][x] == '2' )
		{
			 y-- ;
		}
		algn_size++ ;

	}
	
	//backtrack to find the optimal path
	int *A ,  *B;
	//alloc
	A = new int [algn_size + 10 ] ; // plus some extra just in case
	B = new int [algn_size + 10 ] ;
	y = PROT_SIZES[a] ;
	x = PROT_SIZES[b] ;
	int ind = 0 ;
	//cout << y << " " << x << "\n" ;
	while( y >= 0 && x >=0 )
	{
		if( y == 0 && x == 0 )
			break;
		if( derivation[y][x] == '0' )
		{
			A[ind] = y-1 ;
			B[ind] = x-1 ;
			ind++;y--; x-- ;
		}
		else if( derivation[y][x] == '1' )
		{
			A[ind] = -99 ;
			B[ind] = x-1 ;
			ind++; x-- ;
		}
		else if( derivation[y][x] == '2' )
		{
			A[ind] = y-1 ;
			B[ind] = -99 ;
			ind++; y-- ;
		}
		//cout << A[ind-1] << " " << B[ind-1] << " " << ind << " " << (int)derivation[y][x] << "|\t" ;
		//exit(0);
	}

	int temp;
	//reverse A
	for( int i = 0 , j = ind - 1 ; i < j ; i++, j-- )
	{
		temp = A[i] ;
		A[i] = A[j] ;
		A[j] = temp;
	}
	//reverse B
	for( int i = 0 , j = ind - 1 ; i < j ; i++, j-- )
	{
		temp = B[i] ;
		B[i] = B[j] ;
		B[j] = temp;
	}

	int ALIGNMENT_SIZE = ind ;
  	//if (!meditate)
	//	cout << "........[ OK ]" ;
	
	void STORE_RESULT_TO_GLOBAL_LIBRARY( int * , int * , int , int , int ) ;
	STORE_RESULT_TO_GLOBAL_LIBRARY( A , B , a , b , ALIGNMENT_SIZE ) ;

	//cout << "in1..." << flush ;
	DE_ALLOC_2D( DP_MATRIX  , (PROT_SIZES[a]+1) ) ;
	//cout << "in2..." << flush ;
	DE_ALLOC_2D( derivation , (PROT_SIZES[a]+1) ) ;
	DE_ALLOC_1D(A) ;
	DE_ALLOC_1D(B) ;
	//cout << "...out" << flush ;
	
	//print alignment
	//if(1 )//&& a == 0 && b == 4 )
	if(0)
	{
		cout << endl ;
		cout << setw(10) << struct_names[a] << ": ";
		for( int i = 0 ; i < ALIGNMENT_SIZE ; i++ )
			if(A[i] != -99 )
				cout  << PROT[a][A[i]].res_name ;
			else
				cout << "-" ;
		cout << endl ;
		cout << setw(10) << struct_names[b] << ": ";
		for( int i = 0 ; i < ALIGNMENT_SIZE ; i++ )
			if(B[i] != -99 )
				cout  << PROT[b][B[i]].res_name ;
			else
				cout << "-" ;
		cout << endl ;
		exit(0);
		
	}
		
	
	
}
float max( float arr[] , char *ind )
{
	float MAX = -9999.99 ;
	for( int i = 0 ; i < 3 ; i++ )
	{
		if( arr[i] > MAX )
		{
			MAX = arr[i] ;
			*ind = 48 + i ;
		}
	}
	return MAX;
}

void STORE_RESULT_TO_GLOBAL_LIBRARY( int *AA , int *BB , int a , int b  , int ALIGNMENT_SIZE ) 
{
	// Adding the result into the Global_Library 
	//int z =0 ;
	int g_index = 0 ;
	//int x , y ;

	g_index = (NSTRUCTS-1)*a + b - ( (a*(a+1)/2) + 1 ) ;
	//cout <<"\n" << g_index << "\n" ;
	Global_Library[g_index].mates = new int* [2] ;
	Global_Library[g_index].mates[0] = new int[ PROT_SIZES[a] ] ;
	Global_Library[g_index].mates[1] = new int[ PROT_SIZES[b] ] ;
	Global_Library[g_index].sizes[0] = PROT_SIZES[a]  ;
	Global_Library[g_index].sizes[1] = PROT_SIZES[b]  ;
	
	//assign GlobalA and B to -99 initially
	for (int i = 0 ; i < PROT_SIZES[a] ; ++i )
		Global_Library[g_index].mates[0][i] = -99 ;
	for (int i = 0 ; i < PROT_SIZES[b] ; ++i )
		Global_Library[g_index].mates[1][i] = -99 ; 
	
	for (int i = 0 ; i < ALIGNMENT_SIZE ; ++i )
	{
		if( AA[i] != -99 || BB[i]!= -99 )
		{
			if( AA[i] != -99 && AA[i] < PROT_SIZES[a] )
				Global_Library[g_index].mates[0][AA[i]] = BB[i] ;
			if( BB[i] != -99  && BB[i] < PROT_SIZES[b] )
				Global_Library[g_index].mates[1][BB[i]] = AA[i] ;
		}
	}
	if(show_result)
	//if(1)
	{
		cout << "\n\n******************\nPairwise structural alignment Result:\n" ;
		cout << "Size of Structure-" << a << "(" << struct_names[a] << ")" << " is " << PROT_SIZES[a] << endl ;
		cout << "Size of Structure-" << b << "(" << struct_names[b] << ")" << " is " << PROT_SIZES[b] << endl ;
		/*
		int loop_cntr = ALIGNMENT_SIZE / LINE_WIDTH ;
		int res_cntr = 0 ;
		int DIFF = 99999 ;
		DIFF = PROT_SIZES[a] > PROT_SIZES[b] ? (PROT_SIZES[a] - PROT_SIZES[b]) : (PROT_SIZES[b] - PROT_SIZES[a]) ;
		for( int i = 0 ; i < ALIGNMENT_SIZE ; i += LINE_WIDTH )
		{
			for( int j = i ; j < i+ LINE_WIDTH && j < ALIGNMENT_SIZE ; j++ )
				if( AA[y] != -99 )
		}
		for( int i = 0 ; i < PROT_SIZES[a] ; i += LINE_WIDTH )
		{
			for ( y = i ; y < i + LINE_WIDTH && y< PROT_SIZES[a] ; ++y )
			{
				//cout<<y << endl ; 
				if( A[ y  ] != -99 )
					cout << setw(4) <<  y  << "(" << prot_residue_chain[a][ y ] << ")";
				else
					cout << setw(7) << "----" ;

			}
			cout <<"\n" ;
			for ( y = i ; y < i + LINE_WIDTH && y < PROT_SIZES[a] ; ++y )
				cout << setw(4) << "|" << "   " ;
			cout <<"\n" ;
			for ( y = i ; y < i + LINE_WIDTH && y < PROT_SIZES[a] ; ++y )
			{
				if( A[ y ] < PROT_SIZES[b]  && A[y] != -99 )
					cout << setw(4) << A[ y  ] << "(" << prot_residue_chain[b][ A[y] ] << ")";
				else
					cout << setw(7) << "----" ;
					
			}
			cout << endl << endl ;
		}
		*/
		cout << "\n\t\t<====================>\nMates:\n" ;
		for( int i = 0 ; i < PROT_SIZES[a] ; i++ )
		{
			if( i % 9 == 0 )
				cout <<"\n" ;
			cout << setw(3)<< i << "(" << setw(3) << Global_Library[g_index].mates[0][i] << ",";
			if( Global_Library[g_index].mates[0][i] == -99 )
				cout << PROT[a][i].res_name << " " << "-" ;
			else
				cout << PROT[a][i].res_name << " "<< PROT[b][Global_Library[g_index].mates[0][i]].res_name; 

			cout << ") " ;
		}
		cout << endl ;
		for( int i = 0 ; i < PROT_SIZES[b] ; i++ )
		{
			if( i % 9 == 0 )
				cout <<"\n" ;
			cout << setw(3)<< i << "(" << setw(3) << Global_Library[g_index].mates[1][i]  ;
			if( Global_Library[g_index].mates[1][i] == -99 )
				cout << PROT[b][i].res_name << " " << "-" ;
			else
				cout << PROT[b][i].res_name << " "<< PROT[a][Global_Library[g_index].mates[1][i]].res_name; 

			cout << ")  " ;
		}
		cout << endl<< "**********************" << endl ;

	}
	


	
}
void DE_ALLOC_int2D( int **array , int size )
{
	for( int i = 0 ; i < size ; i++ )
		delete[](array[i]) ;
	delete[](array) ;
}



void COPY_EDGE_WEIGHTS_TO_GLOBAL_DEFINITION( int a , int b , float **edge_weights ) 
{
	int ind =  (NSTRUCTS-1)*a + b - ( a*(a+1)/2 + 1 ) ;
	// Allocate space for edge_weights
	Edge_Weights[ind] = new float* [PROT_SIZES[a]] ;
	for(int i = 0 ; i < PROT_SIZES[a] ; i++ )
		Edge_Weights[ind][ i ] = new float[PROT_SIZES[b]] ;
	// initialize to zero
	for(int i = 0 ; i < PROT_SIZES[a] ; i++ )
		for( int j = 0 ; j < PROT_SIZES[b] ; j++ )
			Edge_Weights[ind][ i ][ j ] = 0 ;
	// copy local edge_weights into global Edge_Weights
	for(int i = 0 ; i < PROT_SIZES[a] ; i++ )
		for( int j = 0 ; j < PROT_SIZES[b] ; j++ )
			Edge_Weights[ind][ i ][ j ] = edge_weights[i][j] ;
}
void DE_ALLOC_float2D( float **array , int size ) 
{
	for( int i = 0 ; i < size ; i++ )
		delete[]( array[i] ) ;
	delete[]( array );
}


/*
template <class Type>
void DE_ALLOC_2D( Type **array , int size )
{
	cout << "here1..." << size << " " << array << flush ;
	for( int i = 0 ; i < size ; i++ )
		delete[]( array[i] ) ;
	cout << "here2..." << flush ;
	delete[]( array );
	cout << "...here3" << flush ;
}
*/

