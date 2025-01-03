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
#include "upgma.h"
#include "de_alloc_routines.h"
#include "globals.h" 
void UPGMA( int NumSeqs , float ** distmat , struct tree_info Tree[]) 
{
	int Tree_Size ;
	int *current_nodes = new int [NumSeqs];
	int curr_NodeList_Size = NumSeqs ;
	float **raw_distmat ;
	raw_distmat = new float* [NumSeqs] ;
	for( int i = 0 ; i < NumSeqs ; i++ )
		raw_distmat[i] = new float [NumSeqs] ;
	
	float **changed_distmat ;
	changed_distmat = new float* [NumSeqs] ;
	for( int i = 0 ; i < NumSeqs ; i++ )
		changed_distmat[i] = new float [NumSeqs] ;
	//float R[NumSeqs] ; 

	
	void INITIALIZE_TREE_TO_LEAF_NODES_UPGMA( int , int* ,  struct tree_info [] , int * ) ;
	void DUPLICATE_DISTMAT_UPGMA( int  , float ** , float ** ) ;
	void CHOOSE_MINIMUM_PAIR_UPGMA( int , float ** , int * , int * ) ;
	void UPDATE_TREE_INFO_AND_DISTMAT_UPGMA( float ** , float**  , int* , int , int [] , struct tree_info [] , int , int , int *) ;
	void PRINT_NODE_INFO_UPGMA( int , int , int [] , struct tree_info[] , float ** , float ** , int ) ;
	
	
	INITIALIZE_TREE_TO_LEAF_NODES_UPGMA( NumSeqs , current_nodes , Tree  , &Tree_Size) ;
	DUPLICATE_DISTMAT_UPGMA( NumSeqs , raw_distmat , distmat ) ;
	//PRINT_NODE_INFO_UPGMA( Tree_Size , curr_NodeList_Size , current_nodes , Tree ) ;
	while( curr_NodeList_Size > 2 )
	{
		int ind1=-1,ind2=-1 ;
		CHOOSE_MINIMUM_PAIR_UPGMA( curr_NodeList_Size , raw_distmat , &ind1 , &ind2 ) ;
		UPDATE_TREE_INFO_AND_DISTMAT_UPGMA( distmat , raw_distmat, &Tree_Size , curr_NodeList_Size , current_nodes , Tree , ind1 , ind2 , &curr_NodeList_Size ) ;
		
	}
	//updating the last two remaining nodes into the Tree
	UPDATE_TREE_INFO_AND_DISTMAT_UPGMA( distmat , raw_distmat , &Tree_Size , curr_NodeList_Size , current_nodes , Tree , 0 , 1 , &curr_NodeList_Size) ;
	//PRINT_NODE_INFO_UPGMA( Tree_Size , curr_NodeList_Size , current_nodes , Tree , raw_distmat , changed_distmat , 1 ) ;
	DE_ALLOC_1D( current_nodes ) ;
	DE_ALLOC_2D( raw_distmat , NumSeqs) ;
	DE_ALLOC_2D( changed_distmat , NumSeqs) ;
}
void INITIALIZE_TREE_TO_LEAF_NODES_UPGMA( int NumSeqs , int curr_nodes[] , struct tree_info Tree[] , int *ts ) 
{
	for( int i = 0 ; i < NumSeqs ; i++ )
	{
		Tree[i].node_num = i ;
		Tree[i].leaf_flag = '1';
		Tree[i].child_node_nums[0] = -99 ; 
		Tree[i].child_node_nums[1] = -99 ; 
		Tree[i].branch_lengths[0] = -999.99 ; 
		Tree[i].branch_lengths[1] = -999.99 ; 
		Tree[i].node_size = 1 ; 
		Tree[i].node_list = new int ; Tree[i].node_list[0] = i ;
		curr_nodes[i] = i ;
	}
	*ts = NumSeqs ;
	
}
void DUPLICATE_DISTMAT_UPGMA( int NumSeqs , float **raw_distmat , float **distmat ) 
{
	for( int i = 0 ; i <  NumSeqs ; i++ )
		for( int j = 0 ; j < NumSeqs ; j++ )
			raw_distmat[i][j] = distmat[i][j];
}

void CHOOSE_MINIMUM_PAIR_UPGMA( int size , float **r_distmat , int *ind1 , int *ind2 ) 
{
	float MIN = r_distmat[0][1] ; *ind1 = 0 ; *ind2 = 1 ;
	
	for( int i = 0 ; i <  size ; i++ )
		for( int j = i+1 ;  j < size ; j++ )
			if(r_distmat[i][j] < MIN)	  
			{
				MIN = r_distmat[i][j] ;
				*ind1 = i ; *ind2=j ;
			}
}
void UPDATE_TREE_INFO_AND_DISTMAT_UPGMA( float **distmat , float **raw_distmat , int *t_index , int size , int current_nodes[] , struct tree_info Tree[] , int ind1 , int ind2 , int *update_size ) 
{
	//updating struct tree_info
	
		Tree[*t_index].node_num = *t_index ;
		Tree[*t_index].leaf_flag = '0';
		Tree[*t_index].child_node_nums[0] = current_nodes[ind1] ; 
		Tree[*t_index].child_node_nums[1] = current_nodes[ind2] ; 
		Tree[*t_index].branch_lengths[0] = raw_distmat[ind1][ind2]/2 ; 
		Tree[*t_index].branch_lengths[1] = raw_distmat[ind1][ind2]/2 ; 
		Tree[*t_index].node_size = Tree[Tree[*t_index].child_node_nums[0]].node_size + 
					   Tree[Tree[*t_index].child_node_nums[1]].node_size ; 
		Tree[*t_index].node_list = new int [ Tree[*t_index].node_size ] ;
		int cntr = 0 ;
		for( int j = 0 ; j < Tree[ Tree[*t_index].child_node_nums[0]].node_size  ; j++ , cntr++) 
			Tree[*t_index].node_list[cntr] = Tree[Tree[*t_index].child_node_nums[0]].node_list[j] ;
		for( int j = 0 ; j < Tree[ Tree[*t_index].child_node_nums[1]].node_size  ; j++ , cntr++) 
			Tree[*t_index].node_list[cntr] = Tree[Tree[*t_index].child_node_nums[1]].node_list[j] ;

		
	       	
	// updating current_nodes 
		int *temp_list = new int [size] ;
		for( int i = 0 , j = 0 ; i < size ; i++ )
		{
			if( i != ind1 && i !=ind2 )
				temp_list[j++] = current_nodes[i] ;
		}
		temp_list[size-2] = *t_index ;
	//copy_back the temp_list to current_list
		for( int  i = 0 ; i < size-1 ; i ++ )
			current_nodes[i] = temp_list[i] ;
	//updating raw_distmat
		float **temp_distmat = new float* [size];
		for( int  i = 0 ; i < size ; i++ )
			temp_distmat[i] = new float [size] ;

		for( int i = 0 , k = 0; i <  size ; i++ )
		{
			if( i == ind1 || i ==ind2 )
				continue ;
			for(int j = 0  , l = 0 ; j < size ; j++ )
			{
				if( j == ind1 || j ==ind2 )
					continue ;
				temp_distmat[k][l++] = raw_distmat[i][j] ;
			}
			k++ ;
		}
		// adding distances wrt to the new node
		for( int i = 0 ; i < size - 2 ; i++ )
		{
			temp_distmat[size-2][i] = 0;
			for( int k = 0 ; k < Tree[*t_index].node_size  ; k++ ) 
				for( int l = 0 ; l < Tree[ current_nodes[i] ].node_size  ; l++ ) 
					temp_distmat[size-2][i] += distmat[ Tree[*t_index].node_list[k] ][ Tree[ current_nodes[i] ].node_list[l] ];
			
			temp_distmat[size-2][i] /= ( Tree[*t_index].node_size + Tree[ current_nodes[i] ].node_size );
			temp_distmat[i][size-2] = temp_distmat[size-2][i] ;
			
		}
		
	//copy_back the temp_distmat to raw_distmat
		for( int i = 0 ; i<size-1 ; i++ )
			for(int j = 0 ; j<size-1 ; j++ )
				raw_distmat[i][j] = temp_distmat[i][j] ;
		
		*t_index = *t_index + 1 ;
		*update_size = *update_size - 1 ;
		DE_ALLOC_1D( temp_list ) ;
		DE_ALLOC_2D( temp_distmat , size ) ;
}
void PRINT_NODE_INFO_UPGMA( int Tree_Size , int size , int current_nodes[] , struct tree_info Tree[]  , float **raw_distmat , float **changed_distmat , int flag ) 
{
	cout << "\n**************" ;
	if(!flag) cout << "BEFORE" ;
	else 	 cout << " AFTER" ;
	cout << "*****************\n" ;
	cout << "\nRaw_Dist_Mat\n" << setprecision(1) ;
	for( int i = 0 ; i < size ; i++ )
	{
		for( int j = 0 ; j < size ; j++ )
		{
			if(i == j )
				cout << setw(4) << ' ' << "---" << setw(4) << ' '  ;
			else
				cout << setw(11) << raw_distmat[i][j] ;
		}
		cout << "\n" ;
	}
	if(!flag)
	{
	cout << "\nChanged_Dist_Mat\n" << setprecision(1) ;
	for( int i = 0 ; i < size ; i++ )
	{
		for( int j = 0 ; j < size ; j++ )
		{
			if(i == j )
				cout << setw(4) << ' ' << "---" << setw(4) << ' '  ;
			else
				cout << setw(11) << changed_distmat[i][j] ;
		}
		cout << "\n" ;
	}
	}
	
	cout << "\nCurrent Node LIST :" ;
	for( int i = 0 ; i < size ; i++ )
		cout << setw(5) << current_nodes[i] ;
	cout << endl ;
	//if(flag) exit(0) ;

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
