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
#include <cmath>

#include "globals.h"
#include "macros.h"
#include "refine_pairalgn.h"
#include "superpose_2.h"
#include "de_alloc_routines.h"

//#define CA_CA_CUTOFF 6.0
//void rotate_vector(float * , const float [][3] , int ) ;
void REFINE_PAIRWISE_ALIGNMENT( int a , int b , float **EW );


float Calc_distance( float *A , float *B ) 
{
	float d = 0 ;
	for( int i = 0 ; i <  3 ; i++ )
		d += pow( (A[i] -B[i] ), 2 ) ;
	
	return( sqrt(d) ) ;
}

void rotate_vector( float *v , const float R[][3] , const int S )
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

void REFINE_PAIRWISE_ALIGNMENT( int a , int b , float **edge_weights )
{
	int g_index = (NSTRUCTS-1)*a + b - ( (a*(a+1)/2) + 1 ) ;
	int dist1 = 0 , dist2 = 0 , dist3 = 0;
	// convert from pairwise matches from "mates"-format to "alignment"-format!
	int ALEN = 0 ;
	int prev_seq2_indx = -1 ;

	for( int i = 0 ; i < PROT_SIZES[a] ; i++ )
	{
		int Y = Global_Library[g_index].mates[0][i] ;
		
		ALEN++ ;
		//may be gaps in the 2nd seq/struct
		if(Y!= -99)
		{
			for( int j = 0 ; j < Y - prev_seq2_indx -1 ; j++ ) 
				ALEN++ ;
		        prev_seq2_indx = Y ;
		}
	}
	//check for C-terminal gaps in 2nd seq/struct
	for( int j = 0 ; j < PROT_SIZES[b] - prev_seq2_indx -1 ; j++ )  ALEN++ ;
	

	int **A2I_hash ; // a convenient hash from alignment-to-(residue)-indices 
	//alloc A2I_hash
	A2I_hash = new int* [2] ;
	for( int i = 0 ; i < 2 ; i++ )
		A2I_hash[i] = new int [ALEN] ;

	 prev_seq2_indx = -1 ;
	 int l = 0 ;
	for( int i = 0 ; i < PROT_SIZES[a] ; i++ )
	{
		int X = i ;
		int Y = Global_Library[g_index].mates[0][i] ;
		
		
		//may be gaps in the 2nd seq/struct
		if(Y!= -99)
		{
			for( int j = 0 ; j < Y - prev_seq2_indx -1 ; j++ ) 
			{
				A2I_hash[0][l] = -99 ;
				A2I_hash[1][l] = prev_seq2_indx + j  + 1 ;
				l++ ;
			}
		        prev_seq2_indx = Y ;
		}
		A2I_hash[0][l] = X ;
		A2I_hash[1][l] = Y ;
		l++ ;

	}
	//check for C-terminal gaps in 2nd seq/struct
	for( int j = 0 ; j < PROT_SIZES[b] - prev_seq2_indx - 1 ; j++ )  
	{
		A2I_hash[0][l] = -99 ;
		A2I_hash[1][l] = prev_seq2_indx + j  + 1 ;
		l++ ;
	}

	if( 0 )
	{
		//print'em
		cout << "\n" ;
		for( int i = 0 ; i < 2 ; i++ )
		{
			for( int j = 0 ; j < ALEN ; j++ )
				cout << setw(4) << A2I_hash[i][j] ;
	
			cout << "\n" ;
		}
		//exit(0);
	}

	const int DIST = 2 ;
	dist1  = dist2 = dist3 = 0 ;
	for( int i = 0 ; i < ALEN ; i++ )
	{
		if( A2I_hash[0][i] == -99 || A2I_hash[1][i] == -99 ) dist1++ ;
		else
		{
			if( dist1 < DIST && A2I_hash[0][i] != 0 && A2I_hash[1][i] != 0  )  {  dist1 = 0 ; continue ; }


			int j , k ;
			dist2 = 0 ;
			for( j = i+1 ; j < ALEN ; j++ )
			{
				if( A2I_hash[0][j] != -99 && A2I_hash[1][j] != -99 ) break ;
				else dist2++ ;
			}
			
			dist3 = 0 ;
			for( k = j+1 ; k < ALEN ; k++ )
			{
				if( A2I_hash[0][k] != -99 && A2I_hash[1][k] != -99 ) break ;
				else dist3++ ;
			}

			//check for runaway matches now!
			
			  //before that a bit of preprocessing
			  //N-terminus
			  if( A2I_hash[0][i] == 0 || A2I_hash[1][i] == 0 ) dist1 =DIST ;
			  //C-terminus
			  if( A2I_hash[0][i] == PROT_SIZES[a]-1 || A2I_hash[1][i] == PROT_SIZES[b]-1 ) dist2 =DIST ;


			  
			if( dist1 >= DIST && dist2 >= DIST ) 
			{
				//remove this edge and penalize the match
				int X = A2I_hash[0][i] ;
				int Y = A2I_hash[1][i] ;
				edge_weights[ X ][ Y ] = -1 ;

				Global_Library[g_index].mates[0][X] = -99 ; 
				Global_Library[g_index].mates[1][Y] = -99 ;

				dist1 += ( dist2 + 1 ) ;
			}

			//couplet case
			if( dist1 >= DIST && ( dist3 >=DIST || k >= ALEN ) )
			{
				//remove this two edges and penalize the match
				int X = A2I_hash[0][i] ;
				int Y = A2I_hash[1][i] ;
				edge_weights[ X ][ Y ] = -1 ;

				Global_Library[g_index].mates[0][X] = -99 ;
				Global_Library[g_index].mates[1][Y] = -99 ;

				if( j < PROT_SIZES[a] )
				{
					int Xprime = A2I_hash[0][j] ;
					int Yprime = A2I_hash[1][j] ;
					edge_weights[ Xprime ][ Yprime ] = -1 ;

					Global_Library[g_index].mates[0][Xprime] = -99 ;
					Global_Library[g_index].mates[1][Yprime] = -99 ;
				}
			}

			dist1 = 0 ;	
		}
	}

	
	// fix terminal mismatch(es) that gravitate towards a longer maximal fragment pair.
	int Lsize = 1 , Rsize = 1 ; 
	int Lstart[2] , Rstart[2] ;
	int MinRepairSize = 6 ;
	for( int I = 0 ; I < ALEN  ; )
	{
		if( A2I_hash[0][I] == -99 || A2I_hash[1][I] == -99 ) 
		{
			I++ ;
			continue ;
		}

		//Left block
		Lstart[0] = A2I_hash[0][I] ;
		Lstart[1] = A2I_hash[1][I] ;
		Lsize = 1 ;
		for( int J = I+1 ; J < ALEN ; J++ ) 
			if( A2I_hash[0][J] == -99 || A2I_hash[1][J] == -99 ) break ;
			else Lsize++ ;

		if( I+Lsize >= ALEN ) break ;
		if( Lsize < MinRepairSize ) 
		{
			I += Lsize ; 
			continue ;
		}

		//if what follows is a string of gaps,
		// if yes, check if they are from any one of the two seqs/structs
		int gap_sign = ( A2I_hash[0][I+Lsize] == -99 ? 0 : 1 )	;
		int gap_size = 1 ;
		int isPlainGap = YES ;
		for( int J = I+Lsize  ; J < ALEN ; J++ )
		{
			if( A2I_hash[0][J] == -99 || A2I_hash[1][J] == -99 )
			{
				if( J == I+Lsize ) gap_sign = ( A2I_hash[0][I+Lsize] == -99 ? 0 : 1 )	;
				else 
				{
					int thisGapSign =( A2I_hash[0][J] == -99 ? 0 : 1 ) ;
					if( gap_sign != thisGapSign ) 
					{
						isPlainGap = NO ;
						break ;
					}
					gap_size++ ;
				}
			}
			else break ;
		}
		if( isPlainGap == NO )
		{
			I += ( Lsize + gap_size ) ;
			continue;
		}
		if( I+Lsize+gap_size >= ALEN) break ;

		//Right block
		Rstart[0] = A2I_hash[0][I+Lsize+gap_size] ;
		Rstart[1] = A2I_hash[1][I+Lsize+gap_size] ;
		Rsize = 1 ;
		for( int J = I+Lsize+gap_size+1 ; J < ALEN ; J++ )
			if( A2I_hash[0][J] == -99 || A2I_hash[1][J] == -99 ) break ;
			else Rsize++ ;
		
		if( Rsize < MinRepairSize )
		{
			I += (Lsize+gap_size) ;
			continue;
		}

		
		//now superpose 4 residue-matches in the end region of Lblock and 4 residue-matches in the start region of the Rblock
		//note, 4 residues barring the terminal two residue-matches in both the blocks. These terminal two are candidates for
		//refinement.
		float **Coords_Set_stationary;    //| these are used to obtain the Rotation matrix
		float **Coords_Set_moving ;       //|

		float **longer_coords_list ;      //| these are used to repair terminal mismatches
		float **shorter_coords_list;      //|
		int *longerlist_indx  , *shorterlist_indx ;
		
		const int NCORE = 4 ;
		const int TERSIZE = 2 ;
		
		int Size = NCORE*2 ;
		//alloc memory
		Coords_Set_stationary = new float* [Size] ;
		Coords_Set_moving = new float* [Size] ;
		for( int i = 0 ; i < Size ; i++ ) 
		{
			Coords_Set_stationary[i] = new float [3] ;
			Coords_Set_moving[i] = new float [3] ;
		}

		int temp_size1 = TERSIZE*2 + gap_size ;
		longer_coords_list = new float* [temp_size1] ;
		longerlist_indx = new int [temp_size1] ;
		for( int i = 0 ; i < temp_size1 ; i++ )
			longer_coords_list[i] = new float [3] ;

		
		int temp_size2 = TERSIZE*2 ;
		shorter_coords_list= new float* [ temp_size2 ] ;
		shorterlist_indx = new int [temp_size2] ;
		for( int i = 0 ; i < temp_size2 ; i++ )
			shorter_coords_list[i] = new float [3] ;
		
		
		//prepare coordinates to superpose
		int s = 0 ;
		for( int x = Lsize - NCORE - TERSIZE ; x < Lsize - TERSIZE ; x++ , s++ )
			for( int y = 0 ; y < 3 ; y ++ )
			{
				Coords_Set_stationary[s][y] = PROT[a][ Lstart[0]+x ].CA_coords[y] ;
				Coords_Set_moving[s][y]     = PROT[b][ Lstart[1]+x ].CA_coords[y] ;
			}
		
		for( int x = TERSIZE ; x < TERSIZE+NCORE ; x++ , s++ )
			for( int y = 0 ; y < 3 ; y ++ )
			{
				Coords_Set_stationary[s][y] = PROT[a][ Rstart[0]+x ].CA_coords[y] ;
				    Coords_Set_moving[s][y] = PROT[b][ Rstart[1]+x ].CA_coords[y] ;
			}

		if( s != Size )
		{
			cerr << "Error! Impossible...this just can't be!\n" ;
			exit(0) ;	
		}
		
		float RMSD;
		float ROTATION_MAT[3][3] ;
		double CM_stationary[3] ,  CM_moving[3] ;
		SUPERPOSE_2( 	
				Coords_Set_moving , Coords_Set_stationary , Size , 
				&RMSD , ROTATION_MAT , CM_moving , CM_stationary 
			   ) ;

		const float RMSD_REFINE_LIMIT  =  1.5 ;
		const float REFINE_DISTANCE  =  3.0 ;
		if( RMSD > RMSD_REFINE_LIMIT ) 
		{
			DE_ALLOC_2D( Coords_Set_stationary , Size ) ;
			DE_ALLOC_2D( Coords_Set_moving , Size ) ;
			DE_ALLOC_2D( longer_coords_list , temp_size1 ) ;
			DE_ALLOC_2D( shorter_coords_list , temp_size2 ) ;
			DE_ALLOC_1D( longerlist_indx ) ;
			DE_ALLOC_1D( shorterlist_indx ) ;

			I+= (Lsize+gap_size) ; 
			continue ;
		}

		//prepare longer_coords_list and shorter_coords_list , longerlist_indx , shorterlist_indx
		s = 0 ;
		//find which one of the 2 structs is longer in this region
		int longind  = ( Rstart[0] - Lstart[0] > Rstart[1] - Lstart[1] ) ? 0 : 1 ;
		int shortind = (longind == 0 ) ? 1 : 0 ;
		int longlist_size , shortlist_size ;
			
		for( int x = Lstart[longind]+Lsize-TERSIZE ; x < Rstart[longind]+TERSIZE ; x++ , s++ )
		{
			if( longind == 0 )
				for( int y = 0 ; y < 3 ; y ++ )
					longer_coords_list[s][y] = PROT[a][x].CA_coords[y] ;
			if( longind == 1 )
				for( int y = 0 ; y < 3 ; y ++ )
					longer_coords_list[s][y] = PROT[b][x].CA_coords[y] ;


			longerlist_indx[s] = x ;
		}
		longlist_size = s ;

		
		s = 0 ;
		for( int x = Lstart[shortind]+Lsize-TERSIZE ; x < Rstart[shortind]+TERSIZE ; x++ , s++ )
		{
			if( shortind == 0 )
				for( int y = 0 ; y < 3 ; y ++ )
					shorter_coords_list[s][y] = PROT[a][x].CA_coords[y] ;
			if( shortind == 1 )
				for( int y = 0 ; y < 3 ; y ++ )
					shorter_coords_list[s][y] = PROT[b][x].CA_coords[y] ;


			shorterlist_indx[s] = x ;
		}
		shortlist_size = s ;

		if( longind == 0 )
		{
			float  temp_coords[3] ;
			for( int  x = 0 ; x < shortlist_size ; x++  )
			{
				//rotate shorter_coords_list
				for( int l = 0 ; l < 3 ; l++ )
					temp_coords[l] = shorter_coords_list[x][l] ;
						
				for( int l = 0 ; l < 3 ; l++ )
					temp_coords[l] = temp_coords[l] - CM_moving[l] ;
				rotate_vector( temp_coords , ROTATION_MAT , 3 ) ;
				
				for( int l = 0 ; l < 3 ; l++ )
					temp_coords[l] = temp_coords[l] + CM_stationary[l] ;
				for( int l = 0 ; l < 3 ; l++ )
					shorter_coords_list[x][l] = temp_coords[l] ;
			}
		}
		else
		{
			float  temp_coords[3] ;
			for( int  x = 0 ; x < longlist_size ; x++  )
			{
				//rotate longer_coords_list
				for( int l = 0 ; l < 3 ; l++ )
					temp_coords[l] = longer_coords_list[x][l] ;
						
				for( int l = 0 ; l < 3 ; l++ )
					temp_coords[l] = temp_coords[l] - CM_moving[l] ;
				rotate_vector( temp_coords , ROTATION_MAT , 3 ) ;
				
				for( int l = 0 ; l < 3 ; l++ )
					temp_coords[l] = temp_coords[l] + CM_stationary[l] ;
				for( int l = 0 ; l < 3 ; l++ )
					longer_coords_list[x][l] = temp_coords[l] ;
			}
		}

		//now, check terminal residue distances and remember nearest residue;
		int indx  = -1 ;
		for( int x = 0 ; x < shortlist_size ; x++ )
		{
			float mindev = 99999 ;
			float deviation ;
			for( int y = indx+1 ; y < longlist_size ; y++ )
			{
				deviation = Calc_distance( longer_coords_list[y] , shorter_coords_list[x] ) ;
				if( deviation < mindev ) 
				{
					indx = y ;
					mindev = deviation ;
				}
			}
			
			//change the matching
			//if( A2I_hash[longind][ longerlist_indx[indx] ] != A2I_hash[shortind][ shorterlist_indx[x] ]  && mindev < REFINE_DISTANCE )
			if( Global_Library[g_index].mates[longind][ longerlist_indx[indx] ] !=  shorterlist_indx[x]   && mindev < REFINE_DISTANCE )
			{
				int X =   longerlist_indx[indx]  ;
				int Y =  shorterlist_indx[x] ;
				//int oldY = Global_Library[g_index].mates[longind][X] ;
				
				int oldXmate = Global_Library[g_index].mates[longind][X] ;
				int oldYmate = Global_Library[g_index].mates[shortind][Y] ;
				/*
				if( oldY != -99 ) 
				{
					Global_Library[g_index].mates[longind][oldY] = -99 ;
					if(longind == 0 ) edge_weights[ X ][ oldY ] = -1 ;
					else              edge_weights[ oldY ][ X ] = -1 ;
				}
				*/
				Global_Library[g_index].mates[longind][X] = Y ;
				Global_Library[g_index].mates[shortind][Y] = X ;
				if( longind )	
				{
					
					float temp = 0 ;
				       if( oldXmate != -99  && oldYmate  != -99 )
				       {
				       		temp = edge_weights[oldXmate][X] >= edge_weights[Y][oldYmate] ? 
						                                    edge_weights[oldXmate][X] : edge_weights[Y][oldYmate] ;
				       }
				       else
				       {
					       if( oldXmate != -99 )
						       temp =  edge_weights[oldXmate][X] ;
					       else if( oldYmate != -99 )
						       temp =  edge_weights[Y][oldYmate] ;
				       }
				       edge_weights[Y][X] = temp+1 ;
				}
				else
				{
					float temp = 0 ;
				       if( oldXmate != -99  && oldYmate  != -99 )
				       {
						temp =  edge_weights[X][oldXmate] >= edge_weights[oldYmate][Y] ? 
						                                     edge_weights[X][oldXmate] : edge_weights[oldYmate][Y] ;
				       }
				       else
				       {
					       if( oldXmate != -99 )
						       temp =  edge_weights[X][oldXmate] ;
					       else if( oldYmate != -99 )
						       temp =  edge_weights[oldYmate][Y] ;
				       }
					edge_weights[X][Y] = temp+1 ;
				}


				if( oldXmate != -99 ) Global_Library[g_index].mates[shortind][oldXmate] = -99 ;
				if( oldYmate != -99 ) Global_Library[g_index].mates[longind][oldYmate] = -99 ;
				
				/*
				if( oldXmate != -99 && oldYmate != -99 )
				{
					cerr << "longind" << longind << ":: " ;
					cerr << "[" << X << "," << oldXmate << "]" ;
					cerr << "[" << Y << "," << oldYmate << "]" ;
					cerr << " replaced with " << X << "," << Y << " in " ;
					cerr << struct_names[a] << " vs. " << struct_names[b] << "\n" ;
				}
				*/
			}

		}
		

		

		DE_ALLOC_2D( Coords_Set_stationary , Size ) ;
		DE_ALLOC_2D( Coords_Set_moving , Size ) ;
		DE_ALLOC_2D( longer_coords_list , temp_size1 ) ;
		DE_ALLOC_2D( shorter_coords_list , temp_size2 ) ;
		DE_ALLOC_1D( longerlist_indx ) ;
		DE_ALLOC_1D( shorterlist_indx ) ;

		I += Lsize + gap_size ;
	}


	//finally check runaways in the N and C terminus.
	// eliminate 

	//Nterminus ;
	int pos1 = 0 , pos2 = 0 ;
	int struc_ind = 0 ;
	int some_flag = ON ;
	for( int j = 0 ; j < 4 ; j++ )
	{
		some_flag = ON ;
		for( int i = pos1 ; i < ALEN ; i++ )
		{
			if( ( A2I_hash[0][i] != -99 && A2I_hash[1][i] != -99)) 
			{
				if ( A2I_hash[0][i] == j || A2I_hash[1][i] == j )   
				{
					if( A2I_hash[0][i] == j ) struc_ind = 0 ;
					else struc_ind = 1 ;
					pos1 = i ;
					break ;
				}
				else if( A2I_hash[0][i] >= 4 && A2I_hash[1][i] >= 4 )
				{
					some_flag = OFF ;
					 break ;
				}
			}
		}
		if( some_flag == OFF ) break ;
		for( int i = pos1+1 ; i < ALEN ; i++ )
		{
			if( ( A2I_hash[0][i] != -99 && A2I_hash[1][i] != -99) && ( A2I_hash[struc_ind][i] >= 4  )  ) 
			{
				pos2 = i ;
				break ;
			}
		}
		if( pos2 -pos1 >=5+3 )
		{

			int X = A2I_hash[0][pos1] ;
			int Y = A2I_hash[1][pos1] ;
			if(X!= -99 && Y!= -99) edge_weights[ X ][ Y ] = -1 ;

			Global_Library[g_index].mates[0][X] = -99 ;
			Global_Library[g_index].mates[1][Y] = -99 ;
		}
		pos1++ ;
	}

	//Cterminus ;
	pos1 = 0 ; pos2 = 0 ;
	struc_ind = 0 ;
	for( int j = 0 ; j < 4 ; j++ )
	{
		some_flag = ON ;
		for( int i = ALEN-pos1-1 ; i >= 0  ; i-- )
		{
			if( ( A2I_hash[0][i] != -99 && A2I_hash[1][i] != -99) )
			{
				if( A2I_hash[0][i] == PROT_SIZES[a]-j-1 || A2I_hash[1][i] == PROT_SIZES[b]-j-1 )   
				{
					if( A2I_hash[0][i] == PROT_SIZES[a]-j-1 ) struc_ind = 0 ;
					else struc_ind = 1 ;
					pos1 = i ;
					break ;
				}
				else if( A2I_hash[0][i] <= PROT_SIZES[a]-5 && A2I_hash[1][i] <= PROT_SIZES[b]-5 )
				{
					some_flag = OFF ;
					 break ;
				}
					
			}
		}
		if( some_flag == OFF ) break ;
		for( int i = pos1-1 ; i >= 0 ; i-- )
		{
			if( ( A2I_hash[0][i] != -99 && A2I_hash[1][i] != -99) )
			{
				if( (struc_ind = 0 && A2I_hash[0][i] <= PROT_SIZES[a] - 5)  ||
				    (struc_ind = 1 && A2I_hash[1][i] <= PROT_SIZES[b] - 5)  
				  ) 
				{
					pos2 = i ;
					break ;
				}
			}
		}
		if( pos1 -pos2 >=5+3 )
		{

			int X = A2I_hash[0][pos1] ;
			int Y = A2I_hash[1][pos1] ;
			if(X!= -99 && Y!= -99)	edge_weights[ X ][ Y ] = -1 ;

			if (X != -99) Global_Library[g_index].mates[0][X] = -99 ;
			if (Y != -99) Global_Library[g_index].mates[1][Y] = -99 ;
		}
		pos1-- ;
	}

	DE_ALLOC_2D( A2I_hash , 2 ) ;

}
