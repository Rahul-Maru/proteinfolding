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
using std::fixed;

#include <fstream>
using std::ifstream;

#include <cstdlib>

#include<math.h>
#include "macros.h"
#include "globals.h"
#include "ew.h"
#include "superpose_2.h"
#include "superpose_weightedRMS.h"
#include "de_alloc_routines.h"
#include "alloc_routines.h"

#define THETA_elastic  0.2 
#define ALPHA  20 
void CHECK_FOR_WRONG_TERMINAL_MATCHING_NEW( int , int, int, int , int , int* , int*) ;
void CHECK_FOR_MISALIGNED_TERMINAL_SECONDARY_STRUCTURAL_ELEMENTS( int, int, int, int , int* ) ;
void  tentative_PairAlign( int , int , float ***edge_weights , int *[2] , int *);
int determineDegreeOfPruning( int , int , int *[2] );
void determinePruningCodes( int  , int  , int *[2] , int *[2] );

struct maxfrag_pointer_linkedlist
{
	struct maximal_fragments *mflib_ptr ;
	struct maxfrag_pointer_linkedlist *link;
};

struct superposable_fragment_group
{
	struct maxfrag_pointer_linkedlist *mf_hdr_ll ;	
	struct maxfrag_pointer_linkedlist *mf_ll_insert_pos ;	
	int nfrags;
	struct superposable_fragment_group *link;
} *cluster_hdr = NULL , *clust_insert_pos = NULL ;

int NCLUSTERS = 0 ;

float CPSim_MAXIMAL_FRAGMENTS( int a,int b, int I, int J , int F_LEN , float *arr ) // quadratic
{
	//reset arr to 0
	for( int i = 0 ; i < F_LEN ; i++ )
		arr[i] = 0.0 ;

	float SIM = 0 ;
	float temp1 = 0 , temp2 = 0 , temp3 = 0 , env_function = 0 ;
	for( int m = 0 ; m < F_LEN ; m++ )
	{
		//for( int n = m + 1 ; n < CP_WINDOW_SIZE ; n++ )
		for( int n = 0 ; n < F_LEN ; n++ )
		{
			if( m == n) continue;
			
			temp1 = 0 , temp2 = 0 , temp3 = 0 , env_function = 0 ;
			temp1 = fabs( distance_matrices[a][I+m][I+n] - distance_matrices[b][J+m][J+n] ) ;
			temp2 = ( distance_matrices[a][I+m][I+n] + distance_matrices[b][J+m][J+n] ) / 2 ;
			if( temp2 == 0 ) temp1 = 0 ;
			else temp1/= temp2 ;
			temp3 = THETA_elastic - temp1 ;
			//env_function = exp(-(pow( (temp2/ALPHA) , 2 ))) ;
			env_function = 1 ;
			if( n > m )
				SIM += (temp3 * env_function ) ;
			arr[m] += ( temp3 * env_function )  ;
		}
	}
	return(SIM) ;
}


float CPSim_MAXIMAL_FRAGMENTS( int a, int b, int I,int J, int K, int L, int R_LEN, int C_LEN, float *R , float *C )// quartic
{
	//reset R and C to 0
	for( int i = 0 ; i < R_LEN ; i++ ) R[i] = 0.0 ;
	for( int i = 0 ; i < C_LEN ; i++ ) C[i] = 0.0 ;

	float SIM = 0 ;
	float temp1 = 0 , temp2 = 0 , temp3 = 0 , env_function = 0 ;
	for( int m = 0 ; m < R_LEN ; m++ )
		for( int n = 0 ; n < C_LEN ; n++)
		{
			temp1 = 0 , temp2 = 0 , temp3 = 0 , env_function = 0 ;
			temp1 = fabs( distance_matrices[a][I+m][K+n] - distance_matrices[b][J+m][L+n] ) ;
			temp2 = ( distance_matrices[a][I+m][K+n] + distance_matrices[b][J+m][L+n] ) / 2 ;
			if( temp2 == 0 ) temp1 = 0 ;
			else temp1/= temp2 ;
			temp3 = THETA_elastic - temp1 ;
			env_function = exp(-(pow( (temp2/ALPHA) , 2 ))) ;
			temp1 = (temp3 * env_function ) ;
			SIM += temp1 ;
			R[m] += temp1 ; C[n] += temp1 ;
		}
	return SIM ;

}
struct overlap_detector
{
	char ind ;
	struct maximal_fragments *ptr ;
};
float Examine_LR_DistProfiles( int , int , float***, float***, float*** , int , int,  int , float* ,float[][3] , double[] , double[] ) ;


void DE_ALLOC_MAXIMAL_FRAGMENT_LIB_indv( int a , int b )
{
	struct  maximal_fragments *ptr1, *ptr2 ;
	//int Size =  ( (NSTRUCTS * (NSTRUCTS-1))/2 ) ;
	int g_index = (NSTRUCTS-1)*a + b - ( (a*(a+1)/2) + 1 ) ;
	//for( int i = 0 ; i < Size ; i++ )
	//{
		ptr1 = Max_Frags_Library[g_index] ;
		while( ptr1 != NULL )
		{
			ptr2 = ptr1 -> link ;
			delete ptr1 ;
			ptr1 = ptr2 ;
		}
//	}
	//delete[] Max_Frags_Library ;
}

/* ***************************************************
   * GLOBAL FUNCTION DEFINITION                      *
   *************************************************** */
void CALCULATE_EDGE_WEIGHTS( int a , int b , float ***edge_weights )
{
	void PRINT_PROT_INFO( int a, int b ) ;
	// PRINT_PROT_INFO( a, b ) ;
	struct maximal_fragments *mf_insrt_point = NULL ;
	
	int g_index = (NSTRUCTS-1)*a + b - ( (a*(a+1)/2) + 1 ) ;
	struct overlap_detector **Overlap_ind ;
	//alloc mem
	Overlap_ind = new struct overlap_detector* [ PROT_SIZES[a] ] ;
	for( int i = 0 ; i <  PROT_SIZES[a] ; i ++ )
		Overlap_ind[i] = new struct overlap_detector [ PROT_SIZES[b] ] ;
	//init
	for( int i = 0 ; i <  PROT_SIZES[a] ; i ++ )
		for( int j = 0 ; j < PROT_SIZES[b] ; j++ )
		{
			Overlap_ind[i][j].ind = '0' ;
			Overlap_ind[i][j].ptr = NULL ;
		}
	
	int T_NUM = 0 ;
	//if( !meditate )
	//	cout << "Calculating Edge Weights....." << flush ;
	//int tmp_cntr = 0 ;
	//float Phi_elastic_function( int , int , int , int , int ,  int , int , int  ); // local function prototype
	// Allocate space for edge_weights
	(*edge_weights) = new float* [PROT_SIZES[a]] ;
	for(int i = 0 ; i < PROT_SIZES[a] ; i++ )
		(*edge_weights)[ i ] = new float [PROT_SIZES[b]] ;
	// initialize to zero
	for(int i = 0 ; i < PROT_SIZES[a] ; i++ )
		for( int j = 0 ; j < PROT_SIZES[b] ; j++ )
			(*edge_weights)[ i ][ j ] = 0 ;


	//find all maximal congruent fragments then calculating weights using DALI SIM-FUNC using 
	// joint superposition method
	

	float **Coords_Set_stationary;
	float **Coords_Set_moving ;
	//alloc sufficient memory to above structures
	int tmp_size= -1;
	if( PROT_SIZES[a] >= PROT_SIZES[b] ) tmp_size = 2*PROT_SIZES[a] ;
	else tmp_size = 2*PROT_SIZES[b] ;
	Coords_Set_stationary = new float* [tmp_size] ;
	Coords_Set_moving = new float* [tmp_size] ;
	for( int i = 0 ; i < tmp_size ; i++ ) 
	{
		Coords_Set_stationary[i] = new float [3] ;
		Coords_Set_moving[i] = new float [3] ;
	}
	float RMSD;
	float ROTATION_MAT[3][3];
	double CM_stationary[3] ,  CM_moving[3] ;

	for( int I = 0 ; I < PROT_SIZES[a] - MIN_CP_WINDOW_SIZE + 1; I++ )
		for( int J = 0 ; J < PROT_SIZES[b] - MIN_CP_WINDOW_SIZE + 1 ; J++ )
		{
			int Extn_Len = 0 ;
			int Overlap_flag = OFF ;
			float prev_rmsd = -1 ;
			float prev_rmsds[4] = {-1,-1,-1,-1} ; int tmp_ind = 0 ;
			int temp_size = 0  ;
			int seq_iden_flag = OFF ;
			
			if( Overlap_ind[I][J].ind == '1' ) 
			{
				Overlap_flag = ON ;
				temp_size = 1 ;
				for( int k = 1 ; I+k < PROT_SIZES[a] && J+k < PROT_SIZES[b] ; k++ )
					if( Overlap_ind[I+k][J+k].ind == '1' ) temp_size++ ;
					else break;
				if( temp_size <= MIN_CP_WINDOW_SIZE ) Extn_Len = 0 ;
				else Extn_Len = temp_size - MIN_CP_WINDOW_SIZE ;
			}
				
			//cout << PROT_SIZES[a] << " " << PROT_SIZES[b] << " " << I << " " << J << "|\n";
			for(  ;	I + MIN_CP_WINDOW_SIZE + Extn_Len -1 < PROT_SIZES[a] &&
				J + MIN_CP_WINDOW_SIZE + Extn_Len -1 < PROT_SIZES[b] ;
				Extn_Len++ )
			{
				//cout << "\t\t" << Extn_Len << std::flush ;
				int NEQUIV = MIN_CP_WINDOW_SIZE + Extn_Len ;
				//copy info into Coords_Set_stationary and Coords_Set_moving
				for( int K = 0 ; K < NEQUIV ; K++ )
				{
					for( int L = 0 ; L < 3 ; L++ )
					{
						Coords_Set_stationary[K][L] =  PROT[a][I+K].CA_coords[L] ;
						Coords_Set_moving[K][L] =  PROT[b][J+K].CA_coords[L] ;
					}
				}
				
				SUPERPOSE_2( 	
						Coords_Set_moving , Coords_Set_stationary , NEQUIV , 
						&RMSD , ROTATION_MAT , CM_moving , CM_stationary 
					    ) ;
				//check if RMSD of this superposition is within the Thershold
				//cout << " " << RMSD <<"\n" ;
				if( RMSD > RMSD_THRESH_MAXIMAL ) break ;
				//special cases
				if( NEQUIV < 5 ) {
					if( NEQUIV == 4 && RMSD > 0.7 )  break ;
					else if( NEQUIV == 5 && RMSD > 0.9 )  break ;
				}
				
				prev_rmsd = RMSD ;
				if( tmp_ind < 4 ) prev_rmsds[tmp_ind++] = RMSD ;
			}
			//cout << I << " " << J << " " << RMSD << " " <<  prev_rmsd << endl ; 
			// find similarity score of this maximal fragment
			int FRAG_LEN = MIN_CP_WINDOW_SIZE + Extn_Len -1 ;
			if( Extn_Len == 0 ) 
			{
				int cntr = 0 ;
				for( int i = 0 ; I+i < PROT_SIZES[a] && J+i < PROT_SIZES[b] && i < MIN_CP_WINDOW_SIZE ; i++ )
					if( PROT[a][I+i].res_name == PROT[b][J+i].res_name ) cntr++ ;
					else break ;
				if( cntr >= 4 )
				{
					FRAG_LEN = cntr ;
					seq_iden_flag = ON ;
				}	
				else	continue ;
			}
			
			if( Overlap_flag == ON ) if( temp_size + (MIN_CP_WINDOW_SIZE/2) >= FRAG_LEN ) continue;
			//if( Overlap_flag == ON ) if( temp_size <= MIN_CP_WINDOW_SIZE/2 ) continue;

			int ss_flag = OFF ;
			if( seq_iden_flag == OFF && ( (FRAG_LEN <= 7 ) ||  (FRAG_LEN < 10 && prev_rmsd > 1.5) ) )
			{
				// check if atleast some continuous stretch of atleast 3 residue matches share same secondary structure.
				// if yes then continue as is else try and prune it
				int cnt  = 0 ;
				int max_cnt = 0 ;
				for( int i = 0 ; i < FRAG_LEN ; i++ )
				{
					int A = SS_identifier[a][I+i].SS_code ;
					int B = SS_identifier[b][J+i].SS_code ;
					if( (A == 0 &&  B == 0) || ( (A == 1 ||A == 2) && (B == 1 || B == 2 ) ) ) 
						cnt++ ;
					else
					{
						if( cnt > max_cnt )
							max_cnt = cnt ;
						cnt = 0 ;
					}
				}	
				max_cnt = cnt ;
				if( FRAG_LEN >= 6 && FRAG_LEN <= 7 && max_cnt >= 3 ) 
				{	
					ss_flag = ON ;
				}
				if( FRAG_LEN > 7 && max_cnt < 3 ) 
				{
					int flag = ON ;
					for( int p = FRAG_LEN - MIN_CP_WINDOW_SIZE ; p >= 0 ; p-- )
						if( prev_rmsds[p] <= 1.0 && prev_rmsds[p] >= 0  )
						{
							if( p == 0  && prev_rmsds[p+1] < 0.9 ) p++ ;
	
							flag = OFF ;
							FRAG_LEN -= ((FRAG_LEN - MIN_CP_WINDOW_SIZE) - p ) ;
							prev_rmsd = prev_rmsds[p] ;
							break;
						}
					if(flag == ON ) continue;
				}
			}
			if( seq_iden_flag == OFF && ss_flag == OFF && FRAG_LEN == 6 && prev_rmsd > 0.6) continue ;
			//if(FRAG_LEN == 7 && prev_rmsd > 1.5) continue ;
			//if(FRAG_LEN == 8 && prev_rmsd > 0.9) continue ;
			
			int adjustment_on_head = 0 ;
			int adjustment_on_tail = FRAG_LEN ;
			//this func tries to eliminate the obvious terminal misalignments inthe maximal fragment pairs.
			//if( I + FRAG_LEN < PROT_SIZES[a] - 10 || J + FRAG_LEN < PROT_SIZES[b] - 10 ) 
			if( seq_iden_flag == OFF )
			{
				CHECK_FOR_MISALIGNED_TERMINAL_SECONDARY_STRUCTURAL_ELEMENTS( I, J , a, b , &FRAG_LEN ) ;
				if( FRAG_LEN < MIN_CP_WINDOW_SIZE ) continue ; 
				if( Overlap_flag == ON ) if( temp_size + (MIN_CP_WINDOW_SIZE/2) >= FRAG_LEN ) continue;
				
				adjustment_on_head = 0 ;
				adjustment_on_tail = FRAG_LEN ;
				CHECK_FOR_WRONG_TERMINAL_MATCHING_NEW( I , J, a, b , FRAG_LEN, 
						                   &adjustment_on_tail , &adjustment_on_head) ;
				FRAG_LEN = adjustment_on_tail ;
			
				if( FRAG_LEN < MIN_CP_WINDOW_SIZE ) continue ; 
				if( Overlap_flag == ON ) if( temp_size + (MIN_CP_WINDOW_SIZE/2) >= FRAG_LEN ) continue;
				//if( Overlap_flag == ON ) if( temp_size <= MIN_CP_WINDOW_SIZE/2 ) continue;
			}
			if( adjustment_on_head > 0  ) 
			{
				if( adjustment_on_head < 5 ) continue;
				else
				{
					FRAG_LEN = adjustment_on_head + 1 ;
					if( Overlap_flag == ON ) if( temp_size + (MIN_CP_WINDOW_SIZE/2) >= FRAG_LEN ) continue;
					adjustment_on_head = 0; 
				}
			}
			
			


	
			float Similarity_Score = 0 ;	
			float *tmp_arr = new float [ FRAG_LEN ] ;

			Similarity_Score = CPSim_MAXIMAL_FRAGMENTS( a,b,I+adjustment_on_head,J+adjustment_on_head,FRAG_LEN,tmp_arr );
			///*
			// overweight maximal fragment pairs that almost cover the length of 
			// of atlest one of the structures. This is necessary as this fragment 
			// will not gain from the ensuing quartic approach while there is a 
			// danger of wrong fragment pairs overwhelming the correct ones.
			if( (float)FRAG_LEN/(float)PROT_SIZES[a] >= 0.8 || 
			    (float)FRAG_LEN/(float)PROT_SIZES[b] >= 0.8 ) // if the FRAG_LEN >= 80% of STRUC(*) LEN
				//Similarity_Score *= (4*pow( (float)(JOINT_RMSD_THRESH_MAXIMAL+1-prev_rmsd),2)) ;
				Similarity_Score *= 7 ;
			if( // to avoid errors in very short fragments 
				( PROT_SIZES[a] < 100 && (float)FRAG_LEN/(float)PROT_SIZES[a] >= 0.45 )  ||
				( PROT_SIZES[b] < 60 && (float)FRAG_LEN/(float)PROT_SIZES[b] >= 0.45 )  
			   )		
				//Similarity_Score *= (pow( (float)(JOINT_RMSD_THRESH_MAXIMAL+1-prev_rmsd),2)) ;
				Similarity_Score *= 7 ;
				
			//	*/
			if( Similarity_Score >= 0 )
			{
				for ( int m  = 0 ; m < FRAG_LEN ; m++ )
				{
					/*
					if( 
						(I+adjustment_on_head+m==13 && J+adjustment_on_head+m==8)||
						(I+adjustment_on_head+m==5 && J+adjustment_on_head+m==1)
					  )
					{
						cout << endl << I << " " << J << " " << FRAG_LEN ;
						cout << " " << PROT[a][I+adjustment_on_head+m].res_name << " " ;
						cout << " " << PROT[b][J+adjustment_on_head+m].res_name << " " ;
						cout << " " << adjustment_on_head << " " << adjustment_on_tail << endl ;
					}
					*/
					//avoid weights to the overlapping part of the maximal fragmentpair; 
					//weight only the extension
					if( Overlap_flag == ON && m < temp_size - adjustment_on_head ) continue ;
					else 
					(*edge_weights)[I+adjustment_on_head+m][J+adjustment_on_head+m] += ( ( Similarity_Score  + tmp_arr[m])) ;
				}

			}				
			DE_ALLOC_1D(tmp_arr) ;
			
			// store or update (if there is an overlap) fragment into 
			// its corresponding maximal_fragments data structure
			//if( Overlap_flag == OFF )
			{
				int flag = OFF ;
				//cout << I << " " << J << " " <<  prev_rmsd << " "<< FRAG_LEN << endl ; 
				if( Max_Frags_Library[g_index] == NULL )
				{
					Max_Frags_Library[g_index] = new struct maximal_fragments ;
					mf_insrt_point = Max_Frags_Library[g_index] ;
					mf_insrt_point -> link = NULL;
				}
				else
				{
					mf_insrt_point->link = new struct maximal_fragments ;
					mf_insrt_point = mf_insrt_point -> link ;
					mf_insrt_point -> link = NULL;
				}
				if( Overlap_flag == OFF )
				{
					mf_insrt_point -> start[0] = I+adjustment_on_head ; 
					mf_insrt_point -> start[1] = J+adjustment_on_head ;
				
					mf_insrt_point -> size = FRAG_LEN ;
					mf_insrt_point -> rmsd = prev_rmsd ;
				}
				else
				{
					if( Overlap_ind[I][J].ptr->size <= FRAG_LEN )
					{
						mf_insrt_point -> start[0] = I+adjustment_on_head ; 
						mf_insrt_point -> start[1] = J+adjustment_on_head ;
					
						mf_insrt_point -> size = FRAG_LEN ;
						mf_insrt_point -> rmsd = prev_rmsd ;

						//remove overlap part from the old fragment
						Overlap_ind[I][J].ptr->size -= temp_size ;
					}
					else
					{
						//remove overlap part from the new fragment
						if( adjustment_on_head > temp_size ) adjustment_on_head -= temp_size ;
						else adjustment_on_head = 0 ;
						
						mf_insrt_point -> start[0] = I + temp_size + adjustment_on_head ; 
						mf_insrt_point -> start[1] = J + temp_size + adjustment_on_head ;
					
						mf_insrt_point -> size = ( FRAG_LEN = FRAG_LEN - temp_size ) ;
						mf_insrt_point -> rmsd = prev_rmsd ;
						flag = ON ;
					}
				}
				
				// fill Overlap_inds corresponding to this fragment
				for( int m = 0 ; m < FRAG_LEN ; m++ )
				{
					//if( I+m == 14 && J+m == 10 )
					//	cout << "HERE\n"  ;
						
					//Overlap_ind[I+adjustment_on_head+m][J+adjustment_on_head+m].ind = '1' ;
					//Overlap_ind[I+adjustment_on_head+m][J+adjustment_on_head+m].ptr =  mf_insrt_point ;
					if( flag == ON )
					{
						Overlap_ind[I+temp_size+m][J+temp_size+m].ind = '1' ;
						Overlap_ind[I+temp_size+m][J+temp_size+m].ptr =  mf_insrt_point ;
					}
					else
					{
						Overlap_ind[I+m][J+m].ind = '1' ;
						Overlap_ind[I+m][J+m].ptr =  mf_insrt_point ;
					}
				}
				T_NUM++ ;
			}
			/*
			else
			{
				//cout << I << " " << J << " " <<  prev_rmsd << " " << FRAG_LEN << "********" << endl ; 
				Overlap_ind[I][J].ptr -> size = Overlap_ind[I][J].ptr -> size + FRAG_LEN - temp_size ;
				Overlap_ind[I][J].ptr -> rmsd = (Overlap_ind[I][J].ptr -> rmsd + prev_rmsd)/2 ;
				// fill Overlap_inds corresponding to this fragment
				for( int m = 0 ; m < FRAG_LEN ; m++ )
				{
					Overlap_ind[I+m][J+m].ind = '1' ;
					Overlap_ind[I+m][J+m].ptr = Overlap_ind[I][J].ptr ;
				}
			}
			*/

		}
	//cout << endl << T_NUM << endl ;
	//exit(0) ;
	int *Tentative_Alignment[2] ;
	int *TAPruningCodes[2] ; // these codes will be used to vary degree of pruning
	int gap_flag = OFF ;
	Tentative_Alignment[0] = new int [PROT_SIZES[a]] ; TAPruningCodes[0] = new int [PROT_SIZES[a]] ;
	Tentative_Alignment[1] = new int [PROT_SIZES[b]] ; TAPruningCodes[1] = new int [PROT_SIZES[b]] ;
	tentative_PairAlign( a , b , edge_weights , Tentative_Alignment ,  &gap_flag) ;

	int PRUNING_WINDOW_MAX = determineDegreeOfPruning( a , b , Tentative_Alignment) ;
	//int PRUNING_WINDOW_MIN = PRUNING_WINDOW_MAX;
	int PRUNING_WINDOW = PRUNING_WINDOW_MAX;
	determinePruningCodes( a , b , Tentative_Alignment , TAPruningCodes ) ;
	//gap_flag = OFF ;

	if( gap_flag == ON )
	{
		int weed_cntr = 0 ;
		int cntr = 1 ;
		/* 
		   Prune the maximal fragment pairs list using the tentative alignment from
		   the quadratic procedure. This helps in (1) speed up the method, (2) reduce
		   the noise from the presence (and quartic weighting of) impractical pairs.
		*/
		
		struct maximal_fragments *mf_ptr1 = NULL , *mf_ptr2 = NULL , *mf_ptr_prev = NULL;
		mf_ptr1 =  Max_Frags_Library[g_index] ;
		if( mf_ptr1 != NULL ){ mf_ptr_prev = mf_ptr1 ; mf_ptr1 = mf_ptr1->link ; }
		
		while( mf_ptr1 != NULL )
		{
			int weedout_flag = NO ;
			for( int i = 0 ; i < mf_ptr1->size ; i++ )
			{
				int X = mf_ptr1->start[0] + i ;
				int Y = mf_ptr1->start[1] + i ;
				int Y_prime = Tentative_Alignment[0][X] ;
				int X_prime = Tentative_Alignment[1][Y] ;
				if( Y_prime == -99 && X_prime == -99 ) continue ;
				
				///*
				else if( Y_prime != -99 && Y_prime < Y && Y - PRUNING_WINDOW > Y_prime ){ 
					weedout_flag = YES ; break; }
				else if( Y_prime != -99 && Y_prime > Y && Y + PRUNING_WINDOW < Y_prime ){ 
					weedout_flag = YES ; break; }
				else if( X_prime != -99 && X_prime < X && X - PRUNING_WINDOW > X_prime ){ 
					weedout_flag = YES ; break; }
				else if( X_prime != -99 && X_prime > X && X + PRUNING_WINDOW < X_prime ){ 
					weedout_flag = YES ; break; }
			}

			if( weedout_flag == YES )
			{
				mf_ptr2 = mf_ptr1->link ;
				delete mf_ptr1 ;
				mf_ptr_prev->link = mf_ptr2 ;
				mf_ptr1 = mf_ptr2 ;
				weed_cntr++ ;
			}
			else
			{
				if(gibberish)
				{
					cout << "Fragment " <<  setw(4) << cntr << ":" ;
					cout << setw(3) << mf_ptr1->start[0] << "-" ;
					cout << setw(3) << mf_ptr1->start[1] << " , ";
					cout << setw(3) << mf_ptr1->size ;
					cout << endl ;
					
				}
				mf_ptr_prev = mf_ptr1 ;
				mf_ptr1 = mf_ptr1->link ;
			}
		}
		//cout << endl << T_NUM - weed_cntr << endl ;

		/*
		   Accumulate over the weights from quadratic procedure the weigths using 
		   the quartic method
		*/
		float overlapX1  = 0 , overlapX2 =0, overlapY1 = 0, overlapY2 = 0 , ovlp = 0 ;
		mf_ptr1 =  Max_Frags_Library[g_index] ;
		while( mf_ptr1 != NULL )
		{
			if(mf_ptr1->size <= 0) {
				mf_ptr1 =mf_ptr1->link;
				continue;
			}
			mf_ptr2 = mf_ptr1 -> link ;
			while( mf_ptr2 != NULL )
			{
				if(mf_ptr2->size <= 0) {
					mf_ptr2 =mf_ptr2->link;
					continue;
				}
				overlapX1  = overlapX2 = overlapY1 = overlapY2 = ovlp = 0 ;
				if( mf_ptr1 -> start[0] <= mf_ptr2->start[0] && 
				    mf_ptr1 -> start[0] +  mf_ptr1->size > mf_ptr2->start[0] 
				  )
				{
					ovlp = mf_ptr1 -> start[0] +  mf_ptr1->size - mf_ptr2->start[0] ;
					overlapX1 = (ovlp*100)/mf_ptr1->size ;
					overlapX2 = (ovlp*100)/mf_ptr2->size ;
				}
				else if( mf_ptr1 -> start[0] > mf_ptr2->start[0] && 
				         mf_ptr2 -> start[0] + mf_ptr2->size > mf_ptr1->start[0] 
				       )
				{
					ovlp = mf_ptr2 -> start[0] +  mf_ptr2->size - mf_ptr1->start[0] ;
					overlapX1 = (ovlp*100)/mf_ptr2->size ;
					overlapX2 = (ovlp*100)/mf_ptr1->size ;
				}

				if( mf_ptr1 -> start[1] <= mf_ptr2->start[1] && 
				    mf_ptr1 -> start[1] +  mf_ptr1->size > mf_ptr2->start[1] 
				  )
				{
					ovlp = mf_ptr1 -> start[1] +  mf_ptr1->size - mf_ptr2->start[1] ;
					overlapY1 = (ovlp*100)/mf_ptr1->size ;
					overlapY2 = (ovlp*100)/mf_ptr2->size ;
				}
				else if( mf_ptr1 -> start[1] > mf_ptr2->start[1] && 
				         mf_ptr2 -> start[1] + mf_ptr2->size > mf_ptr1->start[1] 
				       )
				{
					ovlp = mf_ptr2 -> start[1] +  mf_ptr2->size - mf_ptr1->start[1] ;
					overlapY1 = (ovlp*100)/mf_ptr2->size ;
					overlapY2 = (ovlp*100)/mf_ptr1->size ;
				}
				if( overlapX1 > 10 || overlapX2 > 10 || overlapY1 > 10 || overlapY2 > 10 )
				{
					mf_ptr2 =mf_ptr2->link;
				    	continue;
				}
					

				int NEQUIV = mf_ptr1->size + mf_ptr2->size ;
				//copy info from mf_ptr1 into Coords_Set_stationary and Coords_Set_moving
				for( int K = 0 ; K < mf_ptr1->size ; K++ )
				{
					for( int L = 0 ; L < 3 ; L++ )
					{
						Coords_Set_stationary[K][L] =  PROT[a][(mf_ptr1->start[0]+K)].CA_coords[L] ;
						Coords_Set_moving[K][L] =  PROT[b][(mf_ptr1->start[1]+K)].CA_coords[L] ;
						
					}
				}			
				
				//append info from mf_ptr2 into Coords_Set_stationary and Coords_Set_moving
				for( int K = mf_ptr1->size , l = 0 ; l < mf_ptr2->size ; K++ , l++ )
				{
					for( int L = 0 ; L < 3 ; L++ )
					{
						Coords_Set_stationary[K][L] =  PROT[a][(mf_ptr2->start[0]+l)].CA_coords[L] ;
						Coords_Set_moving[K][L] =  PROT[b][(mf_ptr2->start[1]+l)].CA_coords[L] ;
					}
				}
				
				SUPERPOSE_2( 	
						Coords_Set_moving , Coords_Set_stationary , NEQUIV , 
						&RMSD , ROTATION_MAT , CM_moving , CM_stationary 
					   ) ;

				//check if RMSD of this superposition is within the Thershold
				if( RMSD <= JOINT_RMSD_THRESH_MAXIMAL ) 
				{
					float Similarity_Score = 0 ;	
					float *tmp_arrR = new float [ mf_ptr1->size ] ;
					float *tmp_arrC = new float [ mf_ptr2->size ] ;

					Similarity_Score = 
						CPSim_MAXIMAL_FRAGMENTS( a,b, mf_ptr1->start[0],mf_ptr1->start[1],
									 mf_ptr2->start[0], mf_ptr2->start[1],
									 mf_ptr1->size, mf_ptr2->size,
									 tmp_arrR , tmp_arrC 
								       );
					if( Similarity_Score < 0 )
					{
						mf_ptr2 = mf_ptr2 -> link ;
						DE_ALLOC_1D(tmp_arrR);
						DE_ALLOC_1D(tmp_arrC);

						continue;
					}

					//assign edge weights to the matched residues 
					for( int m = 0 ; m < mf_ptr1->size ; m++ )
					{
						(*edge_weights)[mf_ptr1->start[0]+m][mf_ptr1->start[1]+m] +=( 
							pow(
								(float)(JOINT_RMSD_THRESH_MAXIMAL+1-RMSD),2
							   ) *
							( Similarity_Score  + tmp_arrR[m])) ;
					}
					
					for( int m = 0 ; m < mf_ptr2->size ; m++ )
					{
						(*edge_weights)[mf_ptr2->start[0]+m][mf_ptr2->start[1]+m] += (
							pow(
								(float)(JOINT_RMSD_THRESH_MAXIMAL+1-RMSD),2
							   ) *
							( Similarity_Score  + tmp_arrC[m])) ;
										}
					DE_ALLOC_1D(tmp_arrR);
					DE_ALLOC_1D(tmp_arrC);

				}

				mf_ptr2 = mf_ptr2 -> link ;
			}
			mf_ptr1 = mf_ptr1 -> link ;
		}
	}
	DE_ALLOC_1D( Tentative_Alignment[0] ) ;
	DE_ALLOC_1D( Tentative_Alignment[1] ) ;
	DE_ALLOC_1D( TAPruningCodes[0] ) ;
	DE_ALLOC_1D( TAPruningCodes[1] ) ;

	
	//Penalize edge-weights = 0.0
	for(int i = 0 ; i < PROT_SIZES[a] ; i++ )
		for( int j = 0 ; j < PROT_SIZES[b] ; j++ )
			if( (*edge_weights)[ i ][ j ] == 0 ) 
				(*edge_weights)[ i ][ j ] = -1000 ;
	
	if( gibberish )
	{
	
		cout << "EdgeWeights: " << std::fixed << endl ;
		
		cout << setw(7) << " " ;
		for(int j = 0 ; j < PROT_SIZES[b] ; j++ )
			cout << setw(8) << PROT[b][j].res_name ;
		cout << "\n" ;
		
		for(int i = 0 ; i < PROT_SIZES[a] ; i++ )
		{
			cout << setw(7) << PROT[a][i].res_name ;
			for( int j = 0 ; j < PROT_SIZES[b] ; j++ )
				cout << setprecision(1) <<setw(7) << (*edge_weights)[i][j] << " " ;
			cout << endl ;
		}
	}
	DE_ALLOC_2D( Coords_Set_stationary , tmp_size );
	DE_ALLOC_2D( Coords_Set_moving , tmp_size );
	DE_ALLOC_2D( Overlap_ind , PROT_SIZES[a] );
	DE_ALLOC_MAXIMAL_FRAGMENT_LIB_indv( a, b ) ;
	
}
void rotate_vector_frag( float *v , const float R[][3] , const int S )
{
	float s[3] ;
	for( int i = 0 ; i < S ; i++ )
	{
		s[i]=v[i] ;
		v[i]=0.0 ;
	}
	for( int i = 0 ; i < S ; i++ )
		for( int j = 0 ; j < S ; j++ )
			v[i] = v[i] + s[j] * R[i][j] ;
}

float distance_frag( float *A , float *B ) 
{
	float D = 0.0 ;
	for( int i = 0 ; i <  3 ; i++ )
		D += pow( (A[i] - B[i]) , 2 ) ;
	return( sqrt(D) );
}
void CHECK_FOR_WRONG_TERMINAL_MATCHING_NEW( int I , int J , int a, int b, int NEQUIV , int *adj_T , int *adj_H) 
{

	int EXAMINATION_SIZE = (int)((float)((NEQUIV*30)/100)+0.5) ;
	int FRAGSIZE = 6  ;
	int FRAGSIZE_alloc =  FRAGSIZE ;
	int CORESIZE = 4 ;
	const float DIST_CUTOFF = 2.75 ;
	const float DIST_LLIMIT = 1.5 ;
	const float DIST_ULIMIT = 3.2 ;
	float RMSD;
	float ROTATION_MAT[3][3] ;
	double CM_stationary[3] ,  CM_moving[3] ; 
	float **Coords_Set_stationary , **Coords_Set_moving ;
	float *weights ;
	ALLOC_2D( &Coords_Set_stationary, FRAGSIZE_alloc , 3 ) ; 
	ALLOC_2D( &Coords_Set_moving, FRAGSIZE_alloc , 3 ) ;
	ALLOC_1D( &weights , FRAGSIZE_alloc ) ;
	float temp_coords[3] ;
	float dists[6] ;
	//int recalc_size ;
	int H_adj_flag = OFF ;
	int T_adj_flag = OFF ;
	int start_pos ;

	//dealing with left-hand terminal
	if( I > 5 || J > 5 ) // Nterminal- NO cut
	{
		start_pos = ( NEQUIV-(2*EXAMINATION_SIZE) ) >= 4 ? 4 : ( NEQUIV-(2*EXAMINATION_SIZE) ) ;
		CORESIZE = start_pos ;
		start_pos += (EXAMINATION_SIZE-1) ;
		FRAGSIZE = ( CORESIZE+EXAMINATION_SIZE ) >= 6 ? 6 : ( CORESIZE+EXAMINATION_SIZE ) ;
		for( int ii = start_pos ; ii-FRAGSIZE+1 >= 0 ; ii-- )
		{
			H_adj_flag = OFF ;
			//copy info into Coords_Set_stationary and Coords_Set_moving
			for( int K = 0 ; K < FRAGSIZE ; K++ )
			{
				for( int L = 0 ; L < 3 ; L++ )
				{
					Coords_Set_stationary[K][L] =  PROT[a][I+ii-K].CA_coords[L] ;
					Coords_Set_moving[K][L]     =  PROT[b][J+ii-K].CA_coords[L] ;
				}
			}
		
			for( int j = 0 ; j < FRAGSIZE ; j++ )
				if(j < CORESIZE ) weights[j] = 100 ;
				else weights[j] = 1 ;
			
			SUPERPOSE_WRMS( 	
					Coords_Set_moving , Coords_Set_stationary , weights , FRAGSIZE , 
					&RMSD , ROTATION_MAT , CM_moving , CM_stationary 
					);
		
			for( int i = 0  ; i < FRAGSIZE ; i++ )
			{
				for( int j = 0 ; j < 3 ; j++ )
					temp_coords[j] = Coords_Set_moving[i][j] ;
			
				for( int j = 0 ; j < 3 ; j++ )
					temp_coords[j] = temp_coords[j] - CM_moving[j] ;

				rotate_vector_frag( temp_coords , ROTATION_MAT , 3 ) ;
				for( int j = 0 ; j < 3 ; j++ )
					temp_coords[j] = temp_coords[j] + CM_stationary[j] ;

				for( int j = 0 ; j < 3 ; j++ )
					Coords_Set_moving[i][j] = temp_coords[j] ;

			}

			//calculating frag dists
			for( int i = 0  ; i < FRAGSIZE ; i++  )
				dists[i] = distance_frag( Coords_Set_stationary[i] , Coords_Set_moving[i] ) ;
		
			//Fix Head
			for( int  j = 1  ; j < FRAGSIZE ; j++  )
				if( dists[j]  >= DIST_CUTOFF ) 
				{
					if( dists[j-1] <= DIST_LLIMIT && dists[j] < DIST_ULIMIT)
					{
						*adj_H = ii-j ; // allow the deviating residue pair and cut after that
						H_adj_flag = ON ;
						break;
					}
					else if( dists[j] >= DIST_ULIMIT )
					{
						*adj_H = ii-j+1 ; // cut exactly before the deviation
						H_adj_flag = ON ;
						break;
		
					}
				}
			if( H_adj_flag == ON) break ;
		}
	}
	if( I + NEQUIV < PROT_SIZES[a] - 10 || J + NEQUIV < PROT_SIZES[b] - 10 ) 
	{
		//dealing with right-hand terminal
		start_pos = ( NEQUIV-(2*EXAMINATION_SIZE) ) >= 4 ? 4 : ( NEQUIV-(2*EXAMINATION_SIZE) ) ;
		CORESIZE = start_pos ;
		start_pos = (NEQUIV-EXAMINATION_SIZE - start_pos) ;
		FRAGSIZE = ( CORESIZE+EXAMINATION_SIZE ) >= 6 ? 6 : ( CORESIZE+EXAMINATION_SIZE ) ;
		for( int ii = start_pos ; ii+FRAGSIZE <= NEQUIV ; ii++ )
		{
			T_adj_flag = OFF ;
			//copy info into Coords_Set_stationary and Coords_Set_moving
			for( int K = 0 ; K < FRAGSIZE ; K++ )
			{
				for( int L = 0 ; L < 3 ; L++ )
				{
					Coords_Set_stationary[K][L] =  PROT[a][I+ii+K].CA_coords[L] ;
					Coords_Set_moving[K][L]     =  PROT[b][J+ii+K].CA_coords[L] ;
				}
			}
		
			for( int j = 0 ; j < FRAGSIZE ; j++ )
				if(j < CORESIZE ) weights[j] = 100 ;
				else weights[j] = 1 ;
			
			SUPERPOSE_WRMS( 	
					Coords_Set_moving , Coords_Set_stationary , weights , FRAGSIZE , 
					&RMSD , ROTATION_MAT , CM_moving , CM_stationary 
					);
		
			for( int i = 0  ; i < FRAGSIZE ; i++ )
			{
				for( int j = 0 ; j < 3 ; j++ )
					temp_coords[j] = Coords_Set_moving[i][j] ;
			
				for( int j = 0 ; j < 3 ; j++ )
					temp_coords[j] = temp_coords[j] - CM_moving[j] ;

				rotate_vector_frag( temp_coords , ROTATION_MAT , 3 ) ;
				for( int j = 0 ; j < 3 ; j++ )
					temp_coords[j] = temp_coords[j] + CM_stationary[j] ;

				for( int j = 0 ; j < 3 ; j++ )
					Coords_Set_moving[i][j] = temp_coords[j] ;

			}

			//calculating frag dists
			for( int i = 0  ; i < FRAGSIZE ; i++  )
				dists[i] = distance_frag( Coords_Set_stationary[i] , Coords_Set_moving[i] ) ;
		
			//Fix Tail 
			for( int  j = 1  ; j < FRAGSIZE ; j++  )
				if( dists[j]  >= DIST_CUTOFF ) 
				{
					if( dists[j-1] <= DIST_LLIMIT && dists[j] < DIST_ULIMIT)
					{
						*adj_T =  ii+j+1-*adj_H ; // allow the deviating residue pair and cut after that.
						T_adj_flag = ON ;
						break;
					}
					else
					{
						*adj_T =  ii+j-*adj_H ; // cut exactly beforethe deviating residue pair.
						T_adj_flag = ON ;
						break;

					}
				}
			if( T_adj_flag == ON) break ;
		}
	}

	if( T_adj_flag == OFF ) *adj_T -= *adj_H ;

	DE_ALLOC_1D(weights) ;
	DE_ALLOC_2D( Coords_Set_stationary , FRAGSIZE_alloc ) ;
	DE_ALLOC_2D( Coords_Set_moving , FRAGSIZE_alloc ) ;
}

void CHECK_FOR_MISALIGNED_TERMINAL_SECONDARY_STRUCTURAL_ELEMENTS( int I, int J , int a, int b , int *LEN )
{
	int FRAG_LEN = *LEN ;
	int A = SS_identifier[a][I+FRAG_LEN-1].SS_code ;
	int B = SS_identifier[b][J+FRAG_LEN-1].SS_code ;
	int Astart = SS_identifier[a][I+FRAG_LEN-1].start, Aend = SS_identifier[a][I+FRAG_LEN-1].end ; 
	int Bstart = SS_identifier[b][J+FRAG_LEN-1].start, Bend = SS_identifier[b][J+FRAG_LEN-1].end ; 
	int repair_flag = OFF ;

	//helix
	if( A == 0 &&  A == B )
	{
		if( Aend - (I+FRAG_LEN-1) >= 2 && Bend - (J+FRAG_LEN-1) >= 2 )
		{
			*LEN =  (Astart - I) >= (Bstart - J) ? (Astart-I): (Bstart-J) ;
			repair_flag = ON ;
		}
	}
	//sheet
	else if( ( A == 1 || A == 2 ) && ( B == 1 || B == 2 ) )
	{
		if( Aend - (I+FRAG_LEN-1) >= 2 && Bend - (J+FRAG_LEN-1) >= 2 )
		{
			*LEN =  (Astart - I) >= (Bstart - J) ? (Astart-I): (Bstart-J) ;
			repair_flag = ON ;
		}
	}


	//if very long MFP(>=25), and terminal part ends in helix or sheet, of length >= 6 but <= 10
	//cut the MFP before the terminal secondary structural element.
	for( int i = 1 ; i <= 3 ; i++  )
	{
		A = SS_identifier[a][I+FRAG_LEN-i].SS_code ; B = SS_identifier[b][J+FRAG_LEN-i].SS_code ;
		Astart = SS_identifier[a][I+FRAG_LEN-i].start; Aend = SS_identifier[a][I+FRAG_LEN-i].end ; 
		Bstart = SS_identifier[b][J+FRAG_LEN-i].start; Bend = SS_identifier[b][J+FRAG_LEN-i].end ; 
		if( ( A == 0  && B == 0)  ||  (( A == 1 || A == 2 ) && ( B == 1 || B == 2 )) ) break ;	
	}
	
	//long MFP, short terminla helix
	if( repair_flag == OFF && FRAG_LEN >=25 && 
	    ( ( A == 0  && B == 0) )	&&
	    ( ( Aend -Astart+1 >= MIN_CP_WINDOW_SIZE) || (Bend-Bstart+1 >= MIN_CP_WINDOW_SIZE ) ) &&
	    ( ( Aend -Astart+1 <= 6) || (Bend-Bstart+1 <= 6 ) )
	  )
	{
		//cout << I << " " << J << " " << *LEN  << "| " ;
		*LEN =  (Astart - I) <= (Bstart - J) ? (Astart-I): (Bstart-J) ;
		//cout << *LEN << endl ;
		repair_flag = ON ;
	}
	//long MFP, short terminal sheet
	else if( repair_flag == OFF && FRAG_LEN >=25 && 
	    ( (( A == 1 || A == 2 ) && ( B == 1 || B == 2 ))  )	&&
	    ( ( Aend -Astart+1 >= MIN_CP_WINDOW_SIZE) || (Bend-Bstart+1 >= MIN_CP_WINDOW_SIZE ) ) &&
	    ( ( Aend -Astart+1 <=10) || (Bend-Bstart+1 <= 10 ) )
	  )
	{
		//cout << I << " " << J << " " << *LEN  << "| " ;
		*LEN =  (Astart - I) <= (Bstart - J) ? (Astart-I): (Bstart-J) ;
		//cout << *LEN << endl ;
		repair_flag = ON ;
	}
}
void  tentative_PairAlign( int a , int b , float ***edge_weights , int *TALGN[2] , int *gflag )
{
	float arr[3] ;
	
	float **DP_MATRIX ;
	char **derivation;
	
	float max_talgn( float [] , char * ) ;
	
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
			arr[0] = DP_MATRIX[y-1][x-1] + (*edge_weights)[y-1][x-1] ;
			
			//horizontal derivation
			arr[1] = DP_MATRIX[y][x-1]  ;
			
			//vertical
			arr[2] = DP_MATRIX[y-1][x]  ;
			
			DP_MATRIX[y][x] = max_talgn( arr , &derivation[y][x] ) ;
		}
	//cout << "out of rest-DP\n" ;
	if( gibberish )
	//if( 1 )
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
	
	// convert A-B alignment into mates format
	for( int i = 0 ; i < ALIGNMENT_SIZE ; i++ )
		if( A[i] == -99 ) 
		{
			if( *gflag == OFF ) *gflag = ON ;
			continue;
		}
		else TALGN[0][ A[i] ]  = B[i] ;
	
	for( int i = 0 ; i < ALIGNMENT_SIZE ; i++ )
		if( B[i] == -99 ) 
		{	
			if( *gflag == OFF ) *gflag = ON ;
			continue;
		}
		else TALGN[1][ B[i] ]  = A[i] ;

	
	
	if(gibberish)
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
		//exit(0);
	}
		
	DE_ALLOC_2D( DP_MATRIX  , (PROT_SIZES[a]+1) ) ;
	DE_ALLOC_2D( derivation , (PROT_SIZES[a]+1) ) ;
	DE_ALLOC_1D( A) ;
	DE_ALLOC_1D( B) ;
	
	
}
float max_talgn( float arr[] , char *ind )
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


int  determineDegreeOfPruning( int a , int b , int *TALGN[2] )
{
	int tot_gaps = 0  ;
	int tot_matches = 0 ;
	int longest_gap_stretch  = 0;
	int gap_start_flag = OFF ;

	int temp_cntr = 0 ;
	for( int i = 0 ; i < PROT_SIZES[a] ; i++ )
	{
		if( TALGN[0][i] == -99 ) 
		{
			if( gap_start_flag == OFF )
			{
				temp_cntr = 1 ;
				gap_start_flag = ON ;
			}
			else
			{
				temp_cntr++ ;
			}

			tot_gaps++ ;
		}
		else 
		{
			if( gap_start_flag == ON && temp_cntr > longest_gap_stretch  )
				longest_gap_stretch = temp_cntr ;
			temp_cntr = 0 ;
			gap_start_flag = OFF ;
			
			tot_matches++ ;
		}
	}


	for( int i = 0 ; i < PROT_SIZES[b] ; i++ )
	{
		if( TALGN[1][i] == -99 ) 
		{
			if( gap_start_flag == OFF )
			{
				temp_cntr = 1 ;
				gap_start_flag = ON ;
			}
			else
				temp_cntr++ ;

			tot_gaps++ ;
		}
		else 
		{
			if( gap_start_flag == ON && temp_cntr > longest_gap_stretch  )
				longest_gap_stretch = temp_cntr ;
			temp_cntr = 0 ;
			gap_start_flag = OFF ;
		}
	}
	
	int alen = tot_gaps + tot_matches ;

	float perc_gaps	 = (float) tot_gaps / (float)alen ;

	if(gibberish)
	{
		cout << endl << endl ;
		cout << "TG:" << tot_gaps << endl ;
		cout << "AL:" << alen << endl ;
		cout << "P:" << perc_gaps << endl ;
		cout << "LL:" << longest_gap_stretch << endl ;
		cout << endl << endl ;
	}

	//if( perc_gaps < 0.05 ) return 5 ;
	if( perc_gaps < 0.1 ) return 15 ;
	else if( perc_gaps < 0.15 ) return 15 ;
	else if( perc_gaps < 0.2 ) return 20 ;
	else if( perc_gaps < 0.25 ) return 25 ;
	else 
	{
		return (longest_gap_stretch > 30 ? longest_gap_stretch+10 : 30) ;
	}
	
}

void determinePruningCodes( int a , int b , int *TALGN[2] , int *TCodes[2] ) 
{
	//codes:
	//0: allow full checks (+ or - PRUNING_WINDOW_MAX)
	//1: belongs to LEFT terminal part of a large contiguous algn; allow MFPs (>-PRUNING_WINDOW_MIN)
	//2: belongs to RIGHT terminal part of a large contiguous algn; allow MFPs (<+PRUNING_WINDOW_MIN)
	//3: belongs to MIDDLE part of a large contiguous algn; ignore MFPs
	for( int i = 0 ; i < PROT_SIZES[a] ; i++ )
		TCodes[0][i] = 0 ;
	for( int i = 0 ; i < PROT_SIZES[b] ; i++ )
		TCodes[1][i] = 0 ;
	
	int x , y , len ;
	int TERMINAL_SIZE = 10 ;
	
	for( int i = 0 ; i < PROT_SIZES[a] ; i++ )
	{
		if( TALGN[0][i] == -99 ) continue ;

		x = i ;
		y = TALGN[0][i] ;
		len = 0 ;
		
		for( int j = i+1 ; j < PROT_SIZES[a] ; j++ )
		{
			if( TALGN[0][j]-y != j-x ) break;
			len++;
		}
		if(len > 25 )
		{
			int cnt = 0;
			for( int j = i ; j < i+len+1 ; j++ )
			{
				if( cnt < TERMINAL_SIZE)
				{
					TCodes[0][j] =  1 ;
					TCodes[1][ TALGN[0][j] ] =  1 ;
				}
				else if ( len-cnt < TERMINAL_SIZE )
				{
					TCodes[0][j] =  2 ;
					TCodes[1][ TALGN[0][j] ] =  2 ;
				}
				else
				{
					TCodes[0][j] =  3 ;
					TCodes[1][ TALGN[0][j] ] =  3 ;
				}
				cnt++ ;
			}
		}
		i += len ;
	}

	if(gibberish)
	{
		for( int i = 0 ; i < PROT_SIZES[a] ; i++ )
			cout << i << "(" << TALGN[0][i] << "- " << TCodes[0][i] << " )  " ;
		cout << endl ;
		exit(0);
	}

}

void PRINT_PROT_INFO( int a, int b ) 
{
	int sizes[2] ; 
	int ind[2] ;
	sizes[0] = PROT_SIZES[a] ; ind[0] = a;
	sizes[1] = PROT_SIZES[b] ; ind[1] = b ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		cout << " ---------------- \n" ;
		for( int j = 0 ; j < sizes[i] ; j++ )
		{
			cout << setw(6) << ios::left << j ;
			for( int L = 0 ; L < 3; L++ )
			{
				cout << setw(7) << setprecision(1) << ios::fixed ;
				cout <<  PROT[ ind[i] ][j].CA_coords[L];
			}
			cout << "\n" ;
		}
		cout << " ---------------- \n\n\n" ;
	}
	exit(0) ;
		
}
