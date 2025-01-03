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
using std::endl ;using std::ios; using std::cerr ;

#include <iomanip>
using std::setw ;
using std::setprecision ;
#include<fstream>
using std::ifstream ;
using std::ofstream ;
#include<cmath>
#include<cstring>
#include <cstdlib>


#include "macros.h"
#include "globals.h"
#include "progress_align.h"
#include "neighbour_joining.h"
#include "upgma.h"
#include "de_alloc_routines.h"
#define LINE_WIDTH 12
float SUBMAT[20][20] ;
char aa_str[] = "ACDEFGHIKLMNPQRSTVWY" ;
float **ProgAlgn_distmat ;
void PROGRESSIVE_ALIGNMENT_USING_EXTENDED_EDGE_WEIGHTS() 
{
	if( !meditate)
		cout << "Final Alignment using the extended weights..." ;
	struct tree_info *Tree = new struct tree_info [2*NSTRUCTS] ;
	//Alloc mem
	ProgAlgn_distmat = new float* [NSTRUCTS] ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		ProgAlgn_distmat[i] = new float [NSTRUCTS] ;
	//Initialize
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		for( int j = 0 ; j < NSTRUCTS ; j++ )
			ProgAlgn_distmat[i][j] = 0 ;
	
	void COMPUTE_DISTANCE_MATRIX_FOR_PROGRESSIVE_ALIGNMENT(float ** );
	void GENERATE_GUIDE_TREE(float ** , struct tree_info*  );
	void PRINT_NODE_INFO(struct tree_info* ) ;
	void PROG_ALGN_along_GuideTree(struct tree_info* ) ;
	
	COMPUTE_DISTANCE_MATRIX_FOR_PROGRESSIVE_ALIGNMENT( ProgAlgn_distmat ) ;
	GENERATE_GUIDE_TREE(  ProgAlgn_distmat , Tree );
	//PRINT_NODE_INFO(Tree ) ;
	PROG_ALGN_along_GuideTree(Tree ) ;
	DE_ALLOC_2D( ProgAlgn_distmat , NSTRUCTS ) ;
	//releasing  memory associated to the Tree
	for( int i = 0 ; i < 2*NSTRUCTS-1 ; i++ )
		delete[] Tree[i].node_list ;
	delete[] Tree ;

	if(!meditate)
	{
		cout << setw(24) << " ";
		cout << "[ \033[32;4mOK\033[0m ]\n" ;
		
	}
	//exit(0);
}
void COMPUTE_DISTANCE_MATRIX_FOR_PROGRESSIVE_ALIGNMENT( float **distmat )
{
	//local func prototype
	float SEQ_ALGN_SCORE( int , int  ) ; // this uses sequence alignments
	float SEQ_ALGN_SCORE_2( int , int  ) ;// this uses pairwise structural alignments over BLOSUM
	

	float *indv_maxs = new float [ (NSTRUCTS * (NSTRUCTS-1))/2 ] ;
	float tot_max = -99999 ; 
	void STR_ALGN_INIT( float * , float* ) ; // this uses calcs the max weights indv'ly and totally so that
	STR_ALGN_INIT( indv_maxs , &tot_max ) ; 
	float STR_ALGN_SCORE( int , int , float * , float  ) ; // this uses MUSTANG's scoring matrix

	
	void PRINT_DISTMAT( int , float ** , int) ;
	void CHANGE_SIMILARITY_TO_DISTANCE( int , float ** ) ;

	for( int i = 0 ; i < NSTRUCTS ; i++ )
		for( int j = i + 1 ; j < NSTRUCTS ; j++ )
		{
			//int a=i , b=j ;
			//int ind = (NSTRUCTS-1)*a + b - ( a*(a+1)/2 + 1 ) ;
			//distmat[i][j] = SEQ_ALGN_SCORE( i , j ) ;
			//distmat[i][j] = SEQ_ALGN_SCORE_2( i , j ) ;
			distmat[i][j] = STR_ALGN_SCORE( i , j , indv_maxs , tot_max ) ;
			
			distmat[j][i] = distmat[i][j] ;
			//cout << "\n" << distmat[j][i] << "\n";
		}
	//PRINT_DISTMAT( NSTRUCTS , distmat , 0 ) ;
	CHANGE_SIMILARITY_TO_DISTANCE( NSTRUCTS , distmat ) ;
	//PRINT_DISTMAT( NSTRUCTS , distmat , 1 ) ;
	DE_ALLOC_1D( indv_maxs) ;
	
}


float SEQ_ALGN_SCORE( int a , int b ) 
{
	int blosum_gap_cost = -5 ;

	float arr[3] ;
	
	float **DP_MATRIX ;
	char **derivation;
	
	float max2( float [] , char * ) ;
	
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
		DP_MATRIX[y][x] = DP_MATRIX[y][x-1] - blosum_gap_cost ;
		derivation[y][x] = '1' ;
	}
	
	for( int x = 0 , y = 1 ; y < PROT_SIZES[a]+1 ; y++ )
	{
		DP_MATRIX[y][x] = DP_MATRIX[y-1][x] - blosum_gap_cost ;
		derivation[y][x] = '2' ;
	}

	//rest of the matrix
	float ew;
	for( int y = 1 ; y < PROT_SIZES[a]+1 ; y++ )
		for( int x = 1 ; x < PROT_SIZES[b]+1 ; x++ )
		{
			int i = 0 , j = 0 ;
			for( i = 0 ; i < 20 ; i++ ) if( PROT[a][y].res_name == aa_str[i] ) break ;
			for( j = 0 ; j < 20 ; j++ ) if( PROT[b][x].res_name == aa_str[j] ) break ;

			if( i == 20 || j == 20 ) ew = 0 ;
			else ew = SUBMAT[i][j] ;
				

			//diagonal derivation
			arr[0] = DP_MATRIX[y-1][x-1] + ew ;
			
			//horizontal derivation
			arr[1] = DP_MATRIX[y][x-1] +  blosum_gap_cost ;
			
			//vertical
			arr[2] = DP_MATRIX[y-1][x] + blosum_gap_cost ;
			
			DP_MATRIX[y][x] = max2( arr , &derivation[y][x] ) ;
		}
	//cout << "out of rest-DP\n" ;
	float cost = DP_MATRIX[ PROT_SIZES[a] ][ PROT_SIZES[b] ] ;

	DE_ALLOC_2D( DP_MATRIX  , (PROT_SIZES[a]+1) ) ;
	DE_ALLOC_2D( derivation , (PROT_SIZES[a]+1) ) ;
	//exit(0);
	return (cost ) ;
}

float SEQ_ALGN_SCORE_2( int a , int b ) 
{
	int ind = (NSTRUCTS-1)*a + b - ( a*(a+1)/2 + 1 ) ;
	struct merge_lib *merge_trav  ;
	merge_trav = Merged_Library[ind] ;
	float score = 0 ;
	for( int i = 0 ; i < merge_trav->sizes[0] ; i++ )
	{
		int A = i ;
		int B = merge_trav->mates[0][i] ;
		if( B == -99 ) continue;

		int x = 0 , y = 0 ; float ew = 0 ;
		for( x = 0 ; x < 20 ; x++ ) if( PROT[a][A].res_name == aa_str[x] ) break ;
		for( y = 0 ; y < 20 ; y++ ) if( PROT[b][B].res_name == aa_str[y] ) break ;

		if( x == 20 || y == 20 ) ew = 0 ;
		else ew = SUBMAT[x][y] ;

		score +=ew ;
	}
	return(score) ;
}

float STR_ALGN_SCORE( int a , int b , float *im , float T ) 
{
	int ind = (NSTRUCTS-1)*a + b - ( a*(a+1)/2 + 1 ) ;
	struct merge_lib *merge_trav  ;
	merge_trav = Merged_Library[ind] ;
	float score = 0 ;
	for( int i = 0 ; i < merge_trav->sizes[0] ; i++ )
	{
		int A = i ;
		int B = merge_trav->mates[0][i] ;
		if( B == -99 ) continue;
		score += ( (Extended_Edge_Weights[ind][A][B]*T)/im[ind] ) ;
	}
	return(score) ;
}


void STR_ALGN_INIT(  float *im , float *t ) 
{
	int cntr  = 0 ;
	*t = -999999 ;
	for( int i =0  ; i <  NSTRUCTS ; i++ )
		for( int j = i+1 ; j <  NSTRUCTS ; j++ )
		{
			im[cntr] = -999999 ;
			for( int  k = 0 ; k < PROT_SIZES[i] ; k++ )
				for( int l = 0 ; l <  PROT_SIZES[j] ; l++ )
				{
					if( Extended_Edge_Weights[cntr][k][l] > im[cntr] )
						im[cntr] = Extended_Edge_Weights[cntr][k][l] ;
				}
			if( im[cntr] > *t )
				*t = im[cntr] ;
			cntr++ ;
		}
}

float max2( float arr[] , char *ind )
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



void PRINT_DISTMAT( int NSTRUCTS , float **distmat , int FLAG) 
{
	cout << "\n============================================\n" ;
	if( FLAG == 0 )
		cout << "Similarity Matrix based on Global pairwise weighted bipartite matching scores:\n" ;
	else
		cout << "Tranformation into Distance Matrix(from Similarity) for guide tree construction:\n" ;
		
	cout << std::setprecision(3) << ios::fixed ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		for( int j = 0 ; j < NSTRUCTS ; j++ )
		{
			if( i == j ) 
				cout << setw(3) << ' ' << "----" << setw(3) << ' '  ;
			else
				cout << setw(10)<< std::setprecision(1) << std::fixed << distmat[i][j] ;
		}
		cout << "\n" ;
	}
}
void CHANGE_SIMILARITY_TO_DISTANCE( int NSTRUCTS , float **distmat ) 
{
	float max = -9999.99 ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		for( int j = i+1 ; j < NSTRUCTS ; j++ )
		{
			if( distmat[i][j] > max )
				max =distmat[i][j] ;
		}

	
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		for( int j = i+1 ; j < NSTRUCTS ; j++ )
		{
			distmat[i][j] = -log(distmat[i][j]/max) ;
			distmat[j][i] = distmat[i][j] ;
		}

	/*
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		for( int j = i+1 ; j < NSTRUCTS ; j++ )
		{
			float SimII = 0.0, SimJJ = 0.0 ;
			for( int k = 0 ; k < PROT_SIZES[i] ; k++ )
			{
				for( int l = 0 ; l < 20 ; l++ ) 
					if ( PROT[i][k].res_name == aa_str[l])
					{
						SimII += SUBMAT[l][l] ;
					       break;	
					}
			}
			
			for( int k = 0 ; k < PROT_SIZES[j] ; k++ )
			{
				for( int l = 0 ; l < 20 ; l++ ) 
					if ( PROT[j][k].res_name == aa_str[l])
					{
						SimJJ += SUBMAT[l][l] ;
					       break;	
					}
			}


			
			//cout << endl <<"****\n" << i << "-" << j << "|\t" << SimII << "," << SimJJ ;
			//cout << "," << (SimII+SimJJ)/2 ;
			//cout << "," << distmat[i][j] << "\n******\n" ;
			distmat[i][j] = distmat[i][j]/((SimII+SimJJ)/2) ;
			distmat[i][j] = -log(distmat[i][j]) ;
			distmat[j][i] = distmat[i][j] ;
		}
		distmat[i][i] = 0.0 ;
	}
	*/
}

void GENERATE_GUIDE_TREE( float **distmat , struct tree_info *Tree )
{
	NEIGHBOUR_JOINING_METHOD( NSTRUCTS , distmat  , Tree ) ;
	//UPGMA( NSTRUCTS , distmat , Tree ) ;
}
void PRINT_NODE_INFO( struct tree_info *Tree ) 
{
	int Tree_Size =  NSTRUCTS + (NSTRUCTS - 1) ;
	cout << "\nTree_status:\n" ;
	for( int i = 0 ; i< Tree_Size  ; i++ )
	{
		cout << setw(4) << Tree[i].node_num  ;
		cout << setw(3) << Tree[i].leaf_flag ;
		cout << setw(5) << Tree[i].child_node_nums[0] << " " << setw(3) << Tree[i].child_node_nums[1] ; 
		cout << setprecision(1) << setw(11) << Tree[i].branch_lengths[0] ; 
		cout << setprecision(1) << setw(11) << Tree[i].branch_lengths[1] << "\n" ; 
	}
		
}
void PROG_ALGN_along_GuideTree( struct tree_info *Tree ) 
{
	void COLLECT_NODE_COMPOSITION_FROM_GUIDE_TREE( struct  node_composition_from_Tree* , struct tree_info* , int ) ;
	void FINAL_ALIGNMENT_PJS_DP( struct node_composition_from_Tree [] , struct tree_info*  );
	struct node_composition_from_Tree *Nodes = new struct node_composition_from_Tree [2*NSTRUCTS-1];
	COLLECT_NODE_COMPOSITION_FROM_GUIDE_TREE( Nodes , Tree , NSTRUCTS ) ;
	if( gibberish)
	//if( 1)
	{
		for( int i = 0 ; i <  2*NSTRUCTS -1 ; i++ )
		{
			cout << "NODE_COMPOSITION INFO:\n" ;
			cout << setw(3) << i <<setw(4) << Nodes[i].num ;
			for( int j = 0 ; j <  Nodes[i].num ; j++ )
				cout << setw(4) << Nodes[i].list[j] ;
			cout << endl ;
		}
	}
	FINAL_ALIGNMENT_PJS_DP( Nodes , Tree ) ;
	//cout << endl << endl << endl ;
	void PRINT_ALIGNMENT_DP() ;
	//PRINT_ALIGNMENT_DP() ;
	//releasing memory
	for( int i = 0 ; i <  2*NSTRUCTS -1 ; i++ )
	{
		if( i < NSTRUCTS )
			delete Nodes[i].list ;
		else
			delete[] Nodes[i].list ;
	}
	delete[] Nodes ;
	
	
}
void COLLECT_NODE_COMPOSITION_FROM_GUIDE_TREE( struct  node_composition_from_Tree *Nodes , struct tree_info *Tree , int NSTRUCTS ) 
{
	for( int i = 0 ; i <  2*NSTRUCTS-1 ; i++ )
	{
		if(i< NSTRUCTS)
		{
			Nodes[i].num = 1 ;
			Nodes[i].list = new int ;
			Nodes[i].list[0] = i ;
		}
		else
		{
			Nodes[i].num = Nodes[Tree[i].child_node_nums[0]].num + Nodes[Tree[i].child_node_nums[1]].num ; 
			Nodes[i].list = new int[Nodes[i].num] ;
			
			int cntr = 0 ;
			for( int j = 0 ; j < Nodes[Tree[i].child_node_nums[0]].num  ; j++ , cntr++) 
				Nodes[i].list[cntr] = Nodes[Tree[i].child_node_nums[0]].list[j] ;
			for( int j = 0  ; j < Nodes[Tree[i].child_node_nums[1]].num  ; j++ , cntr++)
				Nodes[i].list[cntr] = Nodes[Tree[i].child_node_nums[1]].list[j] ;
		}
	}
	
}

void FINAL_ALIGNMENT_PJS_DP( struct node_composition_from_Tree *Nodes , struct tree_info *Tree  )
{
	int **tree_ALGN; 
	int *tree_ALGN_SIZES; 
	//Alloc
	tree_ALGN = new int* [NSTRUCTS] ;
	tree_ALGN_SIZES = new int [NSTRUCTS] ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		tree_ALGN[i] = new int [ PROT_SIZES[i] + 10 ] ;// plus some padding just in case!
	//init to resp. indices for the start
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		for( int j = 0 ; j < PROT_SIZES[i] ; j++ )
			tree_ALGN[i][j] = j ;
		tree_ALGN_SIZES[i] = PROT_SIZES[i] ;
	}
		
	
	int depth = NSTRUCTS ;
	while( depth < 2*NSTRUCTS -1 )
	{
		int nodeA ,  nodeB ;
		nodeA =  Tree[depth].child_node_nums[0] ;
		nodeB =  Tree[depth].child_node_nums[1] ;
		int algn_size1 = tree_ALGN_SIZES[ Nodes[nodeA].list[0] ] ;
	       	int algn_size2 = tree_ALGN_SIZES[ Nodes[nodeB].list[0] ] ; //size of alignment in each cluster
		float **edge_weights;
		void CALC_EDGE_WEIGHTS_for_FinAlgn( int , int , int , int  , struct node_composition_from_Tree [] , int** , float *** ) ;
		void CALC_EDGE_WEIGHTS_for_FinAlgn_nonSP( int , int , int , int  , struct node_composition_from_Tree [] , int** , float *** ) ;
		CALC_EDGE_WEIGHTS_for_FinAlgn( nodeA , nodeB , algn_size1 , algn_size2  , Nodes , tree_ALGN , &edge_weights ) ;
		//CALC_EDGE_WEIGHTS_for_FinAlgn_nonSP( nodeA , nodeB , algn_size1 , algn_size2  , Nodes , tree_ALGN , &edge_weights ) ;
		void DYNAMIC_PROGRAMMING_METHOD_for_FinAlgn( int  , int  , float ** , int ***, int *) ;
		int **result_pairalgn, result_size;
		//cout << "FinAlgn:" << nodeA << " " << nodeB << endl ;
		DYNAMIC_PROGRAMMING_METHOD_for_FinAlgn( algn_size1 , algn_size2 , edge_weights , &result_pairalgn , &result_size) ;
		DE_ALLOC_2D( edge_weights , algn_size1 ) ;
		//updating ALGN for representatives from each cluster. note that though reps are used
		// the edge_weights of the matching comes from taking into account all pairwise edge_weights
		// between every pair of points in each cluster.
		// ...and Absorb the changes into other members of the 2 clusters

		for( int j = 0 ; j < Nodes[nodeA].num ; j++ )
		{
			int *temp_algn = new int [result_size] ;
			for( int i = 0 ; i < result_size ; i++ )
			{
				if( result_pairalgn[0][i] == -99 )
					temp_algn[i] = -99;
				else
					temp_algn[i] =  tree_ALGN[ Nodes[nodeA].list[j] ][ result_pairalgn[0][i]] ;
			}
			//release memory for this array
			delete[] tree_ALGN[ Nodes[nodeA].list[j] ] ;
			//alloc new size
			tree_ALGN[ Nodes[nodeA].list[j] ] = new int [result_size] ;
			for( int i = 0 ; i < result_size ; i++ )
				tree_ALGN[ Nodes[nodeA].list[j] ][i] = temp_algn[i] ;
			tree_ALGN_SIZES[ Nodes[nodeA].list[j] ] = result_size ;
			delete[] temp_algn ;
		}
		
		for( int j = 0 ; j < Nodes[nodeB].num ; j++ )
		{
			int *temp_algn = new int [result_size] ;
			for( int i = 0 ; i < result_size ; i++ )
			{
				if( result_pairalgn[1][i] == -99 )
					temp_algn[i] = -99;
				else
					temp_algn[i] =  tree_ALGN[ Nodes[nodeB].list[j] ][ result_pairalgn[1][i]] ;
			}
			//release memory for this array
			delete[] tree_ALGN[ Nodes[nodeB].list[j] ] ;
			//alloc new size
			tree_ALGN[ Nodes[nodeB].list[j] ] = new int [result_size] ;
			for( int i = 0 ; i < result_size ; i++ )
				tree_ALGN[ Nodes[nodeB].list[j] ][i] = temp_algn[i] ;
			tree_ALGN_SIZES[ Nodes[nodeB].list[j] ] = result_size ;
			delete[] temp_algn ;
		}
		DE_ALLOC_2D( result_pairalgn , 2 ) ;
		
		depth++;
			
	}
	//copy tree_ALGN -> I_ALGN
	//alloc I_ALGN
	I_ALGN = new int* [NSTRUCTS] ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		I_ALGN[i] = new int [ tree_ALGN_SIZES[0] ] ;
	//copy
	for( int i = 0 ; i < NSTRUCTS ;i ++ )
		for( int j = 0 ; j < tree_ALGN_SIZES[0] ; j++ )
			I_ALGN[i][j] = tree_ALGN[i][j];
	ALGN_LEN = tree_ALGN_SIZES[0] ;

	//calculate ALGN (of symbols)
	//alloc ALGN
	ALGN = new char* [NSTRUCTS] ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		ALGN[i] = new char [ ALGN_LEN + 1] ;
	//calc
	for( int i = 0 ; i < NSTRUCTS ;i ++ )
	{
		for( int j = 0 ; j < ALGN_LEN ; j++ )
		{
			if( I_ALGN[i][j] != -99 ) ALGN[i][j] = PROT[i][I_ALGN[i][j]].res_name;
			else ALGN[i][j] = '-' ;
				
		}
		ALGN[i][ALGN_LEN] = '\0' ;
	}

	DE_ALLOC_2D( tree_ALGN , NSTRUCTS ) ;
	DE_ALLOC_1D( tree_ALGN_SIZES) ;

}
void CALC_EDGE_WEIGHTS_for_FinAlgn( int nodeA , int nodeB , int algn_size1 , int algn_size2 , struct node_composition_from_Tree *Nodes , int **tree_ALGN , float ***EW )
{
	//alloc memory for EW
	(*EW) = new float* [ (algn_size1) ] ;
	for( int i = 0 ; i < (algn_size1) ; i++ )
		(*EW)[i] = new float [ (algn_size2) ] ;
	//initialize EW
	for( int i = 0 ; i < (algn_size1) ; i++ )
		for( int  j = 0 ; j < (algn_size2) ; j++ )
			(*EW)[i][j] = 0.0 ;

	for( int i = 0 ; i < algn_size1 ; i++ )
		for( int  j = 0 ; j < algn_size2 ; j++ )
		{
			int val1, val2 , t_index ;
			for(int k = 0 ; k < Nodes[nodeA].num ; k++ )
				for( int l = 0 ; l < Nodes[nodeB].num ; l++ )
				{
					int temp_ind1 = Nodes[nodeA].list[k] ;
					int temp_ind2 = Nodes[nodeB].list[l] ;
					
					if( Nodes[nodeA].list[k] < Nodes[nodeB].list[l] )
					{
						t_index = (NSTRUCTS-1)*temp_ind1 + temp_ind2 - ( (temp_ind1*(temp_ind1+1)/2) + 1 ) ;
						val1 = tree_ALGN[ Nodes[nodeA].list[k] ][i] ;
						val2 = tree_ALGN[ Nodes[nodeB].list[l] ][j] ;
						if( val1 != -99 &&  val2 != -99 )
							(*EW)[i][j] += Extended_Edge_Weights[t_index][val1][val2] ; 
					}
					else
					{
						t_index = (NSTRUCTS-1)*temp_ind2 + temp_ind1 - ( (temp_ind2*(temp_ind2+1)/2) + 1 ) ;
						val1 = tree_ALGN[ Nodes[nodeA].list[k] ][i] ;
						val2 = tree_ALGN[ Nodes[nodeB].list[l] ][j] ;
						if( val1 != -99 &&  val2 != -99 )
							(*EW)[i][j] += Extended_Edge_Weights[t_index][val2][val1] ; 
					}
				}
		}
}
void CALC_EDGE_WEIGHTS_for_FinAlgn_nonSP( int nodeA , int nodeB , int algn_size1 , int algn_size2 , struct node_composition_from_Tree *Nodes , int **tree_ALGN , float ***EW )
{
	//alloc memory for EW
	(*EW) = new float* [ (algn_size1) ] ;
	for( int i = 0 ; i < (algn_size1) ; i++ )
		(*EW)[i] = new float [ (algn_size2) ] ;
	//initialize EW
	for( int i = 0 ; i < (algn_size1) ; i++ )
		for( int  j = 0 ; j < (algn_size2) ; j++ )
			(*EW)[i][j] = 0.0 ;
	
	//find most closest elements one from each cluster to act as representatives for alignment along guidetree
	int ind1 = 0 , ind2 = 0 ;
	float MIN = 999999999 ;
	for(int k = 0 ; k < Nodes[nodeA].num ; k++ )
		for( int l = 0 ; l < Nodes[nodeB].num ; l++ )
		{
			if ( ProgAlgn_distmat[ Nodes[nodeA].list[k] ][ Nodes[nodeB].list[l] ] < MIN )
			{
				MIN = ProgAlgn_distmat[ Nodes[nodeA].list[k] ][ Nodes[nodeB].list[l] ] ;
				ind1 = Nodes[nodeA].list[k] ;
				ind2 = Nodes[nodeB].list[l] ;
			}
			
		}


	for( int i = 0 ; i < algn_size1 ; i++ )
		for( int  j = 0 ; j < algn_size2 ; j++ )
		{
			int val1, val2 , t_index ;
					
					if( ind1 < ind2 )
					{
						t_index = (NSTRUCTS-1)*ind1 + ind2 - ( (ind1*(ind1+1)/2) + 1 ) ;
						val1 = tree_ALGN[ ind1 ][i] ;
						val2 = tree_ALGN[ ind2 ][j] ;
						if( val1 != -99 &&  val2 != -99 )
							(*EW)[i][j] += Extended_Edge_Weights[t_index][val1][val2] ; 
					}
					else
					{
						t_index = (NSTRUCTS-1)*ind2 + ind1 - ( (ind2*(ind2+1)/2) + 1 ) ;
						val1 = tree_ALGN[ ind1 ][i] ;
						val2 = tree_ALGN[ ind2 ][j] ;
						if( val1 != -99 &&  val2 != -99 )
							(*EW)[i][j] += Extended_Edge_Weights[t_index][val2][val1] ; 
					}
		}
}


void DYNAMIC_PROGRAMMING_METHOD_for_FinAlgn( int a  , int b  , float **edge_weights , int ***algn , int *algn_size ) 
{
	float arr[3] ;
	
	float **DP_MATRIX ;
	char **derivation;
	
	float max2( float [] , char * ) ;
	
	// Allocate sufficient space for the DP matrix
	DP_MATRIX = new float* [ a + 1 ] ;
	derivation = new char* [ a + 1 ] ;

	for( int i = 0 ; i <  a + 1 ; i ++ )
	{
		DP_MATRIX[i]  = new float [ b + 1  ] ;
		derivation[i] = new char  [ b + 1  ] ;
	}
	
	//boundary conditions
	DP_MATRIX[0][0]  = 0 ; 
	derivation[0][0] = '0' ;
	
	for( int x = 1 , y = 0 ; x < b + 1 ; x++ )
	{
		DP_MATRIX[y][x] = 0 ;
		derivation[y][x] = '1' ;
	}
	
	for( int x = 0 , y = 1 ; y < a + 1 ; y++ )
	{
		DP_MATRIX[y][x] = 0 ;
		derivation[y][x] = '2' ;
	}

	//rest of the matrix
	for( int y = 1 ; y < a + 1 ; y++ )
		for( int x = 1 ; x < b + 1 ; x++ )
		{
			//diagonal derivation
			arr[0] = DP_MATRIX[y-1][x-1] + edge_weights[y-1][x-1] ;
			
			//horizontal derivation
			arr[1] = DP_MATRIX[y][x-1]  ;
			
			//vertical
			arr[2] = DP_MATRIX[y-1][x]  ;
			
			DP_MATRIX[y][x] = max2( arr , &derivation[y][x] ) ;
		}
	//cout << "out of rest-DP\n" ;
	
	//find algn_len
	int LEN = 0 ;
	int y = a ;
	int x = b ;
	while( y >= 0 && x >=0 )
	{
		if( y == 0 && x == 0 )
			break;
		LEN++ ;
		if( derivation[y][x] == '0' )
		{
			y--; x-- ;
		}
		else if( derivation[y][x] == '1' )
			x-- ;
		else if( derivation[y][x] == '2' )
			y-- ;
	}
	//alloc space for algn
	(*algn) = new int* [2] ;
	for( int i = 0 ; i < 2 ; i++ )
		(*algn)[i] = new int [LEN+10] ; //plus some padding!

	//backtrack to find the optimal path
	y = a ;
	x = b ;
	int ind = 0 ;
	//cout << y << " " << x << "\n" ;
	while( y >= 0 && x >=0 )
	{
		if( y == 0 && x == 0 )
			break;
		if( derivation[y][x] == '0' )
		{
			(*algn)[0][ind] = y-1 ;
			(*algn)[1][ind] = x-1 ;
			ind++;y--; x-- ;
		}
		else if( derivation[y][x] == '1' )
		{
			(*algn)[0][ind] = -99 ;
			(*algn)[1][ind] = x-1 ;
			ind++; x-- ;
		}
		else if( derivation[y][x] == '2' )
		{
			(*algn)[0][ind] = y-1 ;
			(*algn)[1][ind] = -99 ;
			ind++; y-- ;
		}
	}

	int temp;
	//reverse A
	for( int i = 0 , j = ind - 1 ; i < j ; i++, j-- )
	{
		temp = (*algn)[0][i] ;
		(*algn)[0][i] = (*algn)[0][j] ;
		(*algn)[0][j] = temp;
	}
	//reverse B
	for( int i = 0 , j = ind - 1 ; i < j ; i++, j-- )
	{
		temp = (*algn)[1][i] ;
		(*algn)[1][i] = (*algn)[1][j] ;
		(*algn)[1][j] = temp;
	}

	//int ALIGNMENT_SIZE = ind ;
	*algn_size = ind ;


	if( gibberish )
	{
		cout << "********" << a << " " << b << "***************\n" ;
		for( int i = 0  ; i < PROT_SIZES[a] ; i++ )
			cout << setw(5) << i << "(" << (*algn)[0][i] << ")" ;
		cout << endl ;
		for( int i = 0  ; i < PROT_SIZES[b] ; i++ )
			cout << setw(5) << i << "(" << (*algn)[1][i] << ")" ;
		cout << endl ;
		cout << "***********************\n" ;
	}
	DE_ALLOC_2D(DP_MATRIX,  a+1 ) ;
	DE_ALLOC_2D(derivation, a+1 ) ;
}

void PRINT_ALIGNMENT_DP()
{
	//PRINT_ALIGNMENT
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		cout << i <<":\t" ;
		for(int j = 0 ; j < ALGN_LEN ; j++ )
			cout << setw(4) << ALGN[i][j] ;
		cout << endl ;
	}
	//exit(0);
}


