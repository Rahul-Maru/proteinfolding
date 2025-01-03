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
using std::endl ;
using std::flush ;
#include <iomanip>
using std::setw ;
using std::setprecision ;

#include "macros.h"
#include "globals.h"
#include "primary_lib_gen.h"
#include "pairwise_global_structalgn.h"
//#include "pairwise_local_structalgns.h"
#include "merge_global_local_libs.h"
/* *******************************************
 * GLOBAL FUNCTION DEFINITION                *
 * ******************************************* */
void PRIMARY_LIBRARY_GENERATION() 
{
	// Local function prototypes
	void GLOBAL_LIBRARY_GENERATION( );
	void LOCAL_LIBRARY_GENERATION( char [][15] , int ) ;
	void PRINT_GLOBAL_LIBRARY( int );
	void PRINT_LOCAL_LIBRARY( int );
	void PRINT_MERGE_LIBRARY(int) ;
	void PRINT_GLOBAL_EDGE_WEIGHTS(int) ;

	
	//alloc space to Global_Library datastructure
	int a_size = ((NSTRUCTS*(NSTRUCTS-1))/2) ;
	Global_Library = new struct glob_lib [a_size] ;
	
	//alloc space to Max_Frags_Library datastructure
	Max_Frags_Library = new struct maximal_fragments* [a_size] ;
	//init	
	for( int i = 0 ; i <  a_size; i++ )  Max_Frags_Library[i] = NULL ; 

	GLOBAL_LIBRARY_GENERATION();
	//LOCAL_LIBRARY_GENERATION();
	MERGE_GLOBAL_AND_LOCAL_LIBRARIES() ;
	/*
	if(gibberish)
//	if(1)
	{
		PRINT_GLOBAL_LIBRARY(NSTRUCTS) ;
		PRINT_LOCAL_LIBRARY(NSTRUCTS) ;
		PRINT_MERGE_LIBRARY(NSTRUCTS) ;
		PRINT_GLOBAL_EDGE_WEIGHTS(NSTRUCTS) ;
	}
	*/
}
void GLOBAL_LIBRARY_GENERATION()
{
	 void ALLOC_MEM_EDGE_WEIGHTS() ; // This func will only alloc the required num of sentinel float** ptrs
	 ALLOC_MEM_EDGE_WEIGHTS() ;   // to the global definition of float ***Edge_Weights .
	 int tot_pairs = ((NSTRUCTS * (NSTRUCTS - 1))/2) ;
	 int pair_cntr = 1 ;
	 if(!meditate)
	 {
		 cout << "Pairwise Comparisons:     "  ;
	 }
	 for( int I = 0 ; I < NSTRUCTS ; I++ )
		 for( int J = I + 1 ; J < NSTRUCTS ; J++ )
		 {
			 if(!meditate)
			 {
				 cout << setw(8) << pair_cntr << " of " << setw(8) << tot_pairs << flush ;

				 
			 }
			 PAIRWISE_GLOBAL_STRUCTURAL_ALIGNMENT( I , J );
			 if(!meditate)
			 {
				 for( int i = 0 ; i < 20 ; i++ ) cout << "\b" ;
				 cout << flush ;
			 }
			 
			 pair_cntr++ ;
			//exit(0);
		 }
	 
	 if(!meditate)
	 {
	 	cout << setw(43) << " ";
		cout << "[ \033[32;4mOK\033[0m ]\n" ;
	 }


}
void ALLOC_MEM_EDGE_WEIGHTS() 
{
	const int size =  (NSTRUCTS*(NSTRUCTS-1))/2 ;
	Edge_Weights = new float** [size] ;
	
}
/*
void LOCAL_LIBRARY_GENERATION()
{
	 for( int I = 0 ; I < NSTRUCTS ; I++ )
		 for( int J = I + 1 ; J < NSTRUCTS ; J++ )
		 {
			 if(verbose)
				 cout << "Dealing with Pair : (" << I << " , " << J << ")" << std::endl ;
			 PAIRWISE_LOCAL_STRUCTURAL_ALIGNMENTS( I , J );
			 //exit(0);
		 }
}
void PRINT_GLOBAL_LIBRARY( int NSTRUCTS )
{
	int temp_limit = NSTRUCTS*(NSTRUCTS-1)/2 ;
	 for( int I = 0 ; I < NSTRUCTS ; I++ )
		 for( int J = I + 1 ; J < NSTRUCTS ; J++ )
		 {
			int a= I,b=J ;
			int g_index = (NSTRUCTS-1)*a + b - ( a*(a+1)/2 + 1 ) ;
			cout << endl<< I <<" " << J << "**********************" << endl ;
			cout << Global_Library[g_index].sizes[0] << " ";
			cout << Global_Library[g_index].sizes[1] << "\n";
			for( int i = 0 ; i < prot_sizes[a] ; i++ )
			{
				if( i % 10 == 0 )
					cout <<"\n" ;
				cout << setw(3)<< i << "(" << setw(3) << Global_Library[g_index].mates[0][i] << ")  " ;
			}
			cout << endl ;
			for( int i = 0 ; i < prot_sizes[b] ; i++ )
			{
				if( i % 10 == 0 )
					cout <<"\n" ;
				cout << setw(3)<< i << "(" << setw(3) << Global_Library[g_index].mates[1][i] << ")  " ;
			}
			
			cout << endl<< "**********************" << endl ;
		}
}
void PRINT_LOCAL_LIBRARY( int NSTRUCTS )
{
		int temp_limit = NSTRUCTS*(NSTRUCTS-1)/2 ;
		cout << "\n-------------------------------------------------------\n" << endl;
		for( int i = 0 ; i < temp_limit ; i++ )
		{
		int temp_cntr = 0 ;
			//cout << i << "\n" ;
		struct loc_lib *traveler = Local_Library[i] ;
		while ( traveler != NULL )
		{
			cout << "{(" << setw(3) <<  traveler -> frag_start[0] << "--" << setw(3) << traveler -> frag_start[1] << ")" ;
			cout << "," << setw(3) << traveler -> frag_size  ;
			cout << "," << setw(4) << setprecision(2) << traveler -> frag_rmsd <<"}  " ;
			traveler = traveler -> link ;
			if( ++temp_cntr % 5 == 0 )
				cout << endl ;
		}
				
		cout << "\n-------------------------------------------------------" << endl;
		}

}
void PRINT_MERGE_LIBRARY( int NSTRUCTS) 
{
	const int lim= NSTRUCTS*(NSTRUCTS-1)/2 ;
	struct merge_lib *mtrav ;
	for( int j = 0 ; j < lim ; j++ )
	{
		mtrav =  Merged_Library[j] ;
		while( mtrav != NULL )
		{
			cout << endl << "++********************" << endl ;
			cout << mtrav->seq_index[0]<< "(" << mtrav->sizes[0] << ") ";
			cout << mtrav->seq_index[1]<< "(" << mtrav->sizes[1] << ")\n";
			for( int i = 0 ; i < mtrav->sizes[0] ; i++ )
			{
				if( i % 10 == 0 )
					cout <<"\n" ;
				cout << setw(3)<< mtrav->res_indices[0][i] << "(" << setw(3) << mtrav->mates[0][i] << ")  " ;
			}
			cout << endl ;
			for( int i = 0 ; i < mtrav->sizes[1] ; i++ )
			{
				if( i % 10 == 0 )
					cout <<"\n" ;
				cout << setw(3)<< mtrav->res_indices[1][i] << "(" << setw(3) << mtrav->mates[1][i] << ")  " ;
			}
			
			cout << endl<< "++********************" << endl ;
			mtrav = mtrav->link ;
		}
		
	//	sleep(2) ;
	}
	cout << "TOTAL NUMBER OF RECORDS IN MERGE-LIB:" << MERGE_SIZE << endl ;

}
void PRINT_GLOBAL_EDGE_WEIGHTS( int NSTRUCTS ) 
{
	 for( int I = 0 ; I < NSTRUCTS ; I++ )
		 for( int J = I + 1 ; J < NSTRUCTS ; J++ )
		 {
			cout << "-------\n" ;
			int a= I,b=J ;
			int g_index = (NSTRUCTS-1)*a + b - ( a*(a+1)/2 + 1 ) ;
			cout << setprecision(1) ;
			for( int i = 0 ; i < prot_sizes[a] ; i++ )
			{
				for( int j = 0 ; j < prot_sizes[b] ; j ++ )
					cout << setw(6) << Edge_Weights[g_index][i][j] ;
				cout << endl ;
			}
			cout << "-------\n" ;
		 }
}
*/
