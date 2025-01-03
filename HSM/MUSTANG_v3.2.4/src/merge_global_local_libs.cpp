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
using std::endl;
#include <iomanip>
using std::setw ;
using std::setprecision ;
#include "merge_global_local_libs.h"
#include "globals.h"
#include "macros.h"
void MERGE_GLOBAL_AND_LOCAL_LIBRARIES() 
{
	const int loop_limit = (NSTRUCTS*(NSTRUCTS-1))/2 ;
	//Alloc space to Merged_Library
	Merged_Library = new struct merge_lib* [loop_limit] ;
	
	// first dump all the entries in global-lib into merged-lib.
	//MERGE_SIZE = 0 ;
	for( int I = 0 ; I < NSTRUCTS ; I ++ )
		for( int J = I+1 ; J< NSTRUCTS ;J++ )
	    	{
			int i =  (NSTRUCTS-1)*I + J - ( I*(I+1)/2 + 1 ) ;	
			Merged_Library[i] = new (struct merge_lib) ;
			Merged_Library[i]->mates[0] = new int [ Global_Library[i].sizes[0] ] ;
			Merged_Library[i]->mates[1] = new int [ Global_Library[i].sizes[1] ] ;
			Merged_Library[i]->res_indices[0] = new int [ Global_Library[i].sizes[0] ] ;
			Merged_Library[i]->res_indices[1] = new int [ Global_Library[i].sizes[1] ] ;
			Merged_Library[i]->sizes[0] = Global_Library[i].sizes[0]  ;
			Merged_Library[i]->sizes[1] = Global_Library[i].sizes[1]  ;
	
			Merged_Library[i]->seq_index[0] = I  ;
			Merged_Library[i]->seq_index[1] = J  ;
			Merged_Library[i]->link = NULL ;
 
			for( int j = 0 ; j < 2 ; j++ )
				for( int k = 0 ; k < Global_Library[i].sizes[j] ; k++ )
				{
					Merged_Library[i]->mates[j][k] = Global_Library[i].mates[j][k] ;
					Merged_Library[i]->res_indices[j][k] = k ;
				}
			//MERGE_SIZE++ ;
		}
	    /*
	// add new records into merged-lib from local-lib iff they are not subsumes in their global-counterparts
	struct merge_lib *merge_trav ;
	for( int I = 0 ; I < NSTRUCTS ; I ++ )
	  for( int J = I+1 ; J< NSTRUCTS ;J++ )
	    {
		int i =  (NSTRUCTS-1)*I + J - ( I*(I+1)/2 + 1 ) ;	
		merge_trav = Merged_Library[i] ;
	//for( int i = 0 ; i < loop_limit ; i++ )
	//{
		struct loc_lib *tmp_hdr, *tmp_trav ;
		int l = MERGE_SIZE ;
		tmp_hdr = NULL ;
		tmp_trav = NULL ;
		tmp_hdr =  Local_Library[i] ;
		tmp_trav = tmp_hdr ;
		while( tmp_trav!=NULL ) 
		{
			int mismatch_flag = FALSE ;
			for( int j = tmp_trav->frag_start[0] , k = 0 ; j < tmp_trav->frag_start[0] + tmp_trav->frag_size  ; j++ , k++  )
			{
				if( Merged_Library[i]->mates[0][j] != tmp_trav->frag_start[1]+k )
				{
					//cout << "HELLO"  << "\n" ;
					mismatch_flag = TRUE;
					break;
				}
			}
			if( mismatch_flag == TRUE )
			{
				merge_trav ->link = new (struct merge_lib) ;
				merge_trav = merge_trav->link;
				merge_trav->link = NULL ;
				merge_trav->mates[0] = new (int)[ tmp_trav->frag_size ] ;
				merge_trav->mates[1] = new (int)[ tmp_trav->frag_size ] ;
				merge_trav->res_indices[0] = new (int)[ tmp_trav->frag_size ] ;
				merge_trav->res_indices[1] = new (int)[ tmp_trav->frag_size ] ;
				merge_trav->sizes[0] = tmp_trav->frag_size ;
				merge_trav->sizes[1] = tmp_trav->frag_size ;
				merge_trav->seq_index[0] = I  ;
				merge_trav->seq_index[1] = J  ;
				for( int j = tmp_trav->frag_start[0] , k = 0 ; j < tmp_trav->frag_start[0] + tmp_trav->frag_size  ; j++ , k++  )
				{
					merge_trav->mates[0][k] = tmp_trav->frag_start[1] + k ;
					merge_trav->res_indices[0][k] = tmp_trav->frag_start[0] + k ;
					
					merge_trav->mates[1][k] = tmp_trav->frag_start[0] + k ;
					merge_trav->res_indices[1][k] = tmp_trav->frag_start[1] + k ;
				}
				MERGE_SIZE++ ;
				l++;
				//merge_trav= merge_trav->link;
			}
			//cout << Merged_Library[l].sizes[0] << endl ;
			tmp_trav =tmp_trav -> link ;
		}
	}
	*/
	
}

