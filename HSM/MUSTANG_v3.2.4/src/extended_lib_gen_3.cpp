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
#include <iomanip>
using std::setw ;
using std::setprecision ;

#include "globals.h"
#include "macros.h"
#include "extended_lib_gen.h"
#include "de_alloc_routines.h"

int Aminus1, Aminus2, Aminus3, Aminus4, Aplus1, Aplus2, Aplus3, Aplus4 ;
int Bmprime1, Bmprime2, Bmprime3, Bmprime4, Bpprime1, Bpprime2, Bpprime3, Bpprime4 ;
/* *******************************************
 * GLOBAL FUNCTION DEFINITION                *
 * ******************************************* */
void EXTENDED_LIBRARY_GENERATION() 
{
	if(!meditate)
		cout << "Extending Edge-Weights thro' intermediate structures..." ;
	void ALLOC_AND_INITIALIZE_EXTENDED_EDGE_WEIGHTS(); 
	ALLOC_AND_INITIALIZE_EXTENDED_EDGE_WEIGHTS();

	for( int I = 0 ; I < NSTRUCTS ; I++ )
		for( int J = I+1 ;  J <  NSTRUCTS ; J++ )
		{
			int ind = (NSTRUCTS-1)*I + J - ( I*(I+1)/2 + 1 ) ;

			//first add original arc weight into extended_edge_weights Library
			struct merge_lib *merge_trav  ;
			merge_trav = Merged_Library[ind] ;
			while( merge_trav != NULL )
			{
				for( int i = 0 ; i < merge_trav->sizes[0] ; i++ )
				{
					int y = merge_trav->res_indices[0][i] ;
					int z = merge_trav->mates[0][y] ;
					if(z == -99) continue; // DUMMY mate 
					
					Extended_Edge_Weights[ind][y][z] = Edge_Weights[ind][y][z] ;
				}
				merge_trav = merge_trav -> link ;
			}


			for( int intermediate = 0  ; intermediate < NSTRUCTS ; intermediate++ )
			{

				if ( intermediate == I || intermediate == J )
					continue ;

				merge_trav = Merged_Library[ind] ;
				while( merge_trav != NULL )
				{
					//adding transitive-arcs derived through extension via intermediate structures
					// ie if A_x-C_z and C_z-B_y only then add arc A_x-B_y
					for( int i = 0 ; i < merge_trav->sizes[0] ; i++ )
					{
						int A =-9, B = -9, C =-9 ;
					       	A = merge_trav->res_indices[0][i] ;
						int B_orig = merge_trav->mates[0][A] ;
						
						int temp_ind_AC ;
						if (intermediate > I )
							 temp_ind_AC = (NSTRUCTS-1)*I + intermediate - ( I*(I+1)/2 + 1 ) ;
						else
							temp_ind_AC =  (NSTRUCTS-1)*intermediate + I - ( intermediate*(intermediate+1)/2 + 1 ) ;
						
						struct merge_lib *temp_merge_trav_AC  ;
						temp_merge_trav_AC = Merged_Library[temp_ind_AC] ;
						while(temp_merge_trav_AC!=NULL)
						{
							if( intermediate > I )
							{
							// to determine if A is a part of this record
							 if( 
							   temp_merge_trav_AC->res_indices[0][0] <= A &&
							   temp_merge_trav_AC->res_indices[0][0]+
							             temp_merge_trav_AC->sizes[0] > A
							   ) 
							 {
							    int tmp_index = A - temp_merge_trav_AC->res_indices[0][0] ;	 
						  	    C = temp_merge_trav_AC->mates[0][tmp_index] ;
							 }
							 else
							 {
								 temp_merge_trav_AC = temp_merge_trav_AC->link ;
								 continue ;
							 }
							}
							else
							{
							// to determine if A is a part of this record
							 if( 
							   temp_merge_trav_AC->res_indices[1][0] <= A &&
							   temp_merge_trav_AC->res_indices[1][0]+
							             temp_merge_trav_AC->sizes[1] > A
							   ) 
							 {
							    int tmp_index = A - temp_merge_trav_AC->res_indices[1][0] ;	 
						  	    C = temp_merge_trav_AC->mates[1][tmp_index] ;
							 }
							 else
							 {
								 temp_merge_trav_AC = temp_merge_trav_AC->link ;
								 continue ;
							 }
							}


							int temp_ind_CB ;
							if (intermediate > J )
								 temp_ind_CB = (NSTRUCTS-1)*J + intermediate - ( J*(J+1)/2 + 1 ) ;
							else
								temp_ind_CB =  (NSTRUCTS-1)*intermediate + J - ( intermediate*(intermediate+1)/2 + 1 ) ;
							struct merge_lib *temp_merge_trav_CB ;
							temp_merge_trav_CB = Merged_Library[temp_ind_CB] ;
							while(temp_merge_trav_CB!=NULL)
							{
							  if( intermediate > J )
							  {
							   // to determine if C is a part of this record
							   if( 
							     temp_merge_trav_CB->res_indices[1][0] <= C &&
							     temp_merge_trav_CB->res_indices[1][0]+
							             temp_merge_trav_CB->sizes[1] > C
							     )  
							    {
							      int tmp_index = C - temp_merge_trav_CB->res_indices[1][0] ; 
						  	      B = temp_merge_trav_CB->mates[1][tmp_index] ;
							    }
							    else
							    {
							  	 temp_merge_trav_CB = temp_merge_trav_CB->link ;
								 continue ;
							    }
							  }
							  else
							  {
							   // to determine if C is a part of this record
							   if( 
							      temp_merge_trav_CB->res_indices[0][0] <= C &&
							      temp_merge_trav_CB->res_indices[0][0]+
							             temp_merge_trav_CB->sizes[0] > C
							     ) 
							     {
							      int tmp_index = C - temp_merge_trav_CB->res_indices[0][0] ;
						  	      B = temp_merge_trav_CB->mates[0][tmp_index] ;
							     }
							    else
							    {
							  	 temp_merge_trav_CB = temp_merge_trav_CB->link ;
								 continue ;
							    }
							   }
							  
							  if( B!= -99 && B_orig != -99 && merge_trav->mates[1][B] != -99 )
							  {
								  float AAA, BBB ;
							      	  AAA = Edge_Weights[ind][A][B_orig] ;
							      	  BBB = Edge_Weights[ind][ merge_trav->mates[1][B] ][B] ;
								  float max = AAA > BBB ? AAA : BBB ;
								  if( Extended_Edge_Weights[ind][A][B] >=  max && max >= 0 )
										Extended_Edge_Weights[ind][A][B] +=1000 ;
								  else if( max >=0 )  Extended_Edge_Weights[ind][A][B] = max ;
							  }
							  ///*
							  else if( B!= -99 && B_orig != -99 && merge_trav->mates[1][B] == -99)
							  {

								  float AAA ;
							      	  AAA = Edge_Weights[ind][A][B_orig] ;
								  float max = AAA ;
								  if( Extended_Edge_Weights[ind][A][B] >=  max-1 && max-1 >= 0 )
										Extended_Edge_Weights[ind][A][B] +=1000 ;
								  else if( max-1 >=0 )  Extended_Edge_Weights[ind][A][B] = max-1 ;
							  }
							  else if( B!= -99 && B_orig == -99 && merge_trav->mates[1][B] != -99)
							  {
								  float  BBB ;
							      	  BBB = Edge_Weights[ind][ merge_trav->mates[1][B] ][B] ;
								  float max = BBB ;
								  if( Extended_Edge_Weights[ind][A][B] >=  max-1 && max-1 >= 0 )
										Extended_Edge_Weights[ind][A][B] +=1000 ;
								  else if( max >=0 )  Extended_Edge_Weights[ind][A][B] = max-1 ;
							  }
							  //*/
							  else if( B!= -99 && B_orig == -99 && merge_trav->mates[1][B] == -99 )
							  {
							      /*
								 The paln here is not to extended an arc between A and B thro' C if
								 A and B are both unmatched in the pairwise alignment. However if 
								 A-1 implies Bdoubleprime from pairwise library which in turn is 
								 equivalent to B-1 then allow the extension.

								 This helps to stop random drift aways in the final multiple alignment.
								 at the same time allow those that are continuous.
							      */
							      Aminus1 = A-1 ;  Aminus2 = A-2 ; Aminus3 = A-3 ; Aminus4 = A-4 ;
							      if(Aminus1>=0) Bmprime1 = merge_trav->mates[0][Aminus1] ; else Bmprime1 = -1 ;
							      if(Aminus2>=0) Bmprime2 = merge_trav->mates[0][Aminus2] ; else Bmprime2 = -1 ;
							      if(Aminus3>=0) Bmprime3 = merge_trav->mates[0][Aminus3] ; else Bmprime3 = -1 ;
							      if(Aminus4>=0) Bmprime4 = merge_trav->mates[0][Aminus4] ; else Bmprime4 = -1 ;
							      
							      int cntm = 0;
							      if( Bmprime1 == B-1 ) cntm++ ;
							      if( Bmprime2 == B-2 ) cntm++ ;
							      if( Bmprime3 == B-3 ) cntm++ ;
							      if( Bmprime4 == B-4 ) cntm++ ;
							      Aplus1 = A+1 ; Aplus2 = A+2 ; Aplus3 = A+3 ;  Aplus4 = A+4 ;
							      if(Aplus1<merge_trav->sizes[0]) Bpprime1 = merge_trav->mates[0][Aplus1] ;  else Bpprime1 = -1 ;
							      if(Aplus2<merge_trav->sizes[0]) Bpprime2 = merge_trav->mates[0][Aplus2] ;  else Bpprime2 = -1 ;
							      if(Aplus3<merge_trav->sizes[0]) Bpprime3 = merge_trav->mates[0][Aplus3] ;  else Bpprime3 = -1 ;
							      if(Aplus4<merge_trav->sizes[0]) Bpprime4 = merge_trav->mates[0][Aplus4] ;  else Bpprime4 = -1 ;
							      int cntp = 0;
							      if( Bpprime1 == B+1 ) cntp++ ;
							      if( Bpprime2 == B+2 ) cntp++ ;
							      if( Bpprime3 == B+3 ) cntp++ ;
							      if( Bpprime4 == B+4 ) cntp++ ;

							      int extend_flag = YES ;
							      int Aprime = merge_trav->mates[1][B] ;
							      if( B_orig == -99 && merge_trav->mates[1][B] == -99 && cntm < 3 && cntp < 3) 
								      extend_flag = NO ;
							      else if( B_orig == -99 && Aprime != -99 ) 
							      {
								      int delta = (A-Aprime) >= 0 ? (A-Aprime) : (Aprime-A) ;
								      if( delta >= 10 ) extend_flag = NO ;
							      }
							      else if( B_orig != -99 && Aprime == -99 ) 
							      {
								      int delta = (B_orig-B) >= 0 ? (B_orig-B) : (B-B_orig) ;
								      if( delta >= 10 ) extend_flag = NO ;
							      }

							      //if(  B_orig != -99 || merge_trav->mates[1][B] != -99 || cntm >= 3 ||cntp >=3 ) 
							      if( extend_flag == YES )
							      {

									if( Extended_Edge_Weights[ind][A][B] >= 0 )
										Extended_Edge_Weights[ind][A][B] ++ ;
									else
										Extended_Edge_Weights[ind][A][B] = 1 ;
							      }
							  }
						       	  temp_merge_trav_CB= temp_merge_trav_CB->link ;       
							}
						       	temp_merge_trav_AC= temp_merge_trav_AC->link ;       
						}

					}
					merge_trav = merge_trav -> link ;
				}
				

			}// end of intermediate-for loop
		}// end of J-for loop
	// Convert all 0.0 edge_weights in the extended library to -1000
	// this will force to induce a gap rather than align residues where there is no signal
	for( int a = 0 ; a < NSTRUCTS ; a ++ )
		for( int b = a+1 ; b < NSTRUCTS ; b++ )
		{
			// alloc 2D float array
			int ind =  (NSTRUCTS-1)*a + b - ( a*(a+1)/2 + 1 ) ;
			for(int k = 0 ; k < PROT_SIZES[a] ; k++ )
				for( int l = 0 ; l < PROT_SIZES[b] ; l++ )
					if( Extended_Edge_Weights[ind][ k ][ l ]  == 0 )
						Extended_Edge_Weights[ind][ k ][ l ] = -1000 ;

		}

	void PRINT_EXTENDED_EDGE_WEIGHTS() ;
	//DE_ALLOC_1D( max_edge_weights) ;
	if(!meditate)
	{
		cout << setw(14) << " ";
		cout << "[ \033[32;4mOK\033[0m ]\n" ;
	}

	
}
void ALLOC_AND_INITIALIZE_EXTENDED_EDGE_WEIGHTS()
{
	const int size =  NSTRUCTS*(NSTRUCTS-1)/2 ;
	Extended_Edge_Weights = new float** [size] ;
	for( int a = 0 ; a < NSTRUCTS ; a ++ )
		for( int b = a+1 ; b < NSTRUCTS ; b++ )
		{
			// alloc 2D float array
			int ind =  (NSTRUCTS-1)*a + b - ( a*(a+1)/2 + 1 ) ;
			Extended_Edge_Weights[ind] = new float* [PROT_SIZES[a]] ;
			for(int k = 0 ; k < PROT_SIZES[a] ; k++ )
				Extended_Edge_Weights[ind][ k ] = new float [PROT_SIZES[b]] ;
			//initialize 2D array to 0
			for(int k = 0 ; k < PROT_SIZES[a] ; k++ )
				for( int l = 0 ; l < PROT_SIZES[b] ; l++ )
					Extended_Edge_Weights[ind][ k ][ l ] = 0 ;

		}
}
void PRINT_EXTENDED_EDGE_WEIGHTS( int NSTRUCTS ) 
{
	 for( int I = 0 ; I < NSTRUCTS ; I++ )
		 for( int J = I + 1 ; J < NSTRUCTS ; J++ )
		 {
			cout << "-------\n" ;
			int a= I,b=J ;
			int g_index = (NSTRUCTS-1)*a + b - ( a*(a+1)/2 + 1 ) ;
			cout << setprecision(1) ;
			for( int i = 0 ; i < PROT_SIZES[a] ; i++ )
			{
				for( int j = 0 ; j < PROT_SIZES[b] ; j ++ )
					cout << setw(7) << Extended_Edge_Weights[g_index][i][j] ;
				cout << endl ;
			}
			cout << "-------\n" ;
		 }
}
