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
#include<iostream>
using std::cin; 	using std::cout; 	using std::cerr; 
using std::flush; 	using std::endl; 	using std::ios;
using std::fixed;


#include<iomanip>
using std::setprecision ; using std::setw ; 

#include<fstream>
using std::ifstream ;
using std::ofstream ;

#include<cstdio> 
#include<cstdlib> 
#include<cctype> 
#include<cstring> 
#include<ctime> 

#include "macros.h"
#include "globals.h"
#include "output_algn.h"
#include "CmdLineParser.h"
#include "read_structures.h"
#include "distmat.h"
#include "sse_RK.h"
#include "primary_lib_gen.h"
#include "extended_lib_gen.h"
#include "progress_align.h"
#include "superpose_on_core.h"
#include "de_alloc_routines.h"


time_t rundate;
//char **struct_paths ;
/* ****************************************************
 * LOCAL FUNCTIONS
 * ****************************************************/
void PRINT_HEADER()
{
	if(!meditate)
	{
		cout << "\n" ;
		cout << "\033[0;37;1;41m";
		cout << "      MUSTANG (" << VERSION << "): A MUltiple STuructural AligNment alGorithm.        \033[0;0m\n" ;
		//cout << "\033[0;49m" ;
		cout << "\033[1;33;4;1;46m";
		cout << " Authors:     " ; 
		cout << "\033[0;1;37;4;1;46m" ;
		cout << "    A S Konagurthu,  J Whisstock, and P J Stuckey, A M Lesk. \033[0m\n";
		//cout << "\033[0m" ;
		cout << "\033[0;33;4;1;46m";
		cout << " Reference:   " ; cout << "\033[0;37;4;1;46m" ;
		cout << "    A S Konagurthu et al. Proteins 64(3):559-574, 2006       \033[0m\n";
		cout << "\033[0m" ;
		cout << "\n\n" ;
	}
}
void PARSE_COMMAND_LINE( int argc , char  **argv )
{
	if(!meditate)
        	cout << "Parsing the Command Line..." << flush ;
        if( argc != 2 )
        {
                cerr << "Wrong Command Line syntax! Use Syntax Given Below:\n <executable> <Info-file>\n" << endl ;
                exit(0) ;
        }
       // check whether the files
       // listed in the Command-Line exists.
       // ABORT IF WRONG
	ifstream infile1( argv[1] , ios::in ) ;
	if( !infile1 )
	{
                cerr << "CMD-LINE ARG ERROR:\n\tFile- " << argv[1] << " does not exist!\n" ; 
		infile1.close();
                exit(2) ;
        }
	if(!meditate)
	{
		cout << setw(56) << " ";
		cout << "[ \033[32;4mOK\033[0m ]\n" ;
		
	}

}

void DE_ALLOC_DISTMAT( )
{
	 for( int I = 0 ; I < NSTRUCTS ; I++ )
	 {
		for( int i = 0 ; i < PROT_SIZES[I] ; i++ )
			delete[] distance_matrices[I][i] ;
		delete[] distance_matrices[I] ;
	 }
	 delete[] distance_matrices ;

}
void DE_ALLOC_GLOBAL_EDGE_WEIGHTS( )
{
	 for( int I = 0 ; I < NSTRUCTS ; I++ )
		 for( int J = I + 1 ; J < NSTRUCTS ; J++ )
		 {
			int a= I,b=J ;
			int g_index = (NSTRUCTS-1)*a + b - ( a*(a+1)/2 + 1 ) ;
			for( int i = 0 ; i < PROT_SIZES[a] ; i++ )
				delete[] Edge_Weights[g_index][i] ;
			delete[] Edge_Weights[g_index] ;
		 }
	 delete[] Edge_Weights ;
}
void DE_ALLOC_EXTENDED_EDGE_WEIGHTS()
{
	 for( int I = 0 ; I < NSTRUCTS ; I++ )
		 for( int J = I + 1 ; J < NSTRUCTS ; J++ )
		 {
			int a= I,b=J ;
			int g_index = (NSTRUCTS-1)*a + b - ( a*(a+1)/2 + 1 ) ;
			for( int i = 0 ; i < PROT_SIZES[a] ; i++ )
				delete[] Extended_Edge_Weights[g_index][i] ;
			delete[] Extended_Edge_Weights[g_index] ;
		 }
	 delete[] Extended_Edge_Weights ;
}
void DE_ALLOC_GLOBAL_LIB()
{
	int Size =  ( (NSTRUCTS * (NSTRUCTS-1))/2 ) ;
	for( int i = 0 ; i < Size ; i++ )
	{
		for( int j = 0 ; j < 2; j++ )
			delete[] Global_Library[i].mates[j] ;
		
		delete[]  Global_Library[i].mates ;
	}
	delete[] Global_Library ;
}
void DE_ALLOC_LOCAL_LIB()
{
	int Size =  ( (NSTRUCTS * (NSTRUCTS-1))/2 ) ;
	for( int i = 0 ; i < Size ; i++ )
	{
		//cout << "\n" << i << " " <<std::flush;
		struct loc_lib *ptr , *nxt ;
		ptr = Local_Library[i] ;
		if(ptr!=NULL)
			nxt = ptr -> link ; 
		while( ptr!=NULL )
		{
		//	cout << "*"<<std::flush ;
			delete ptr ;
			ptr = nxt ;
			if(ptr!=NULL)
				nxt = ptr -> link ; 
		}
	}
	delete[] Local_Library ;
}
void DE_ALLOC_MERGE_LIB()
{
	 int Size =  ( (NSTRUCTS * (NSTRUCTS-1))/2 ) ;
	for( int i = 0 ; i < Size ; i++ )
	{
		struct merge_lib *ptr , *nxt ;
		ptr = Merged_Library[i] ;
		nxt = ptr -> link ; 
		while( ptr!=NULL )
		{
			delete[] ptr->mates[0] ; delete[] ptr->mates[1] ;
			delete[] ptr->res_indices[0] ; delete[] ptr->res_indices[1] ; 
			delete ptr ;
			ptr = nxt ;
			if(ptr!=NULL)
				nxt = ptr -> link ; 
		}
	}
	delete[] Merged_Library ;
}

void DE_ALLOC_MAXIMAL_FRAGMENT_LIB()
{
	delete[] Max_Frags_Library ;
}



//void ASSIGN_SECONDARY_STR_IDENTIFIERS() 
void SSE_RAMA() 
{
	//alloc secondary_stuct_identifier ;
	secondary_stuct_identifier = new int* [NSTRUCTS];
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		secondary_stuct_identifier[i] = new int [PROT_SIZES[i]] ;
	//initialize it!
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		for( int j = 0 ; j < PROT_SIZES[i] ; j++ ) secondary_stuct_identifier[i][j] = -1 ;
	
	int region_found = NO ;
	int Eregion_indx ;
	//int AA_indx ;
	//int submat_indx ;

	//find region in ramachandran plot.
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		for( int j = 0 ; j < PROT_SIZES[i] ; j++ )
		{
			Eregion_indx = -1 ;
			//compute Eregion_indx
			region_found = NO ;
		       	for( int l = 0 ; l < 7 ; l++ )
			{
				if(
					efimov_regions[l].Phi_Ulimit >= PROT[i][j].Phi &&
					efimov_regions[l].Phi_Llimit <= PROT[i][j].Phi &&
					efimov_regions[l].Psi_Ulimit >= PROT[i][j].Psi &&
					efimov_regions[l].Psi_Llimit <= PROT[i][j].Psi
				  )
				{
					region_found = YES ;
					Eregion_indx = efimov_regions[l].Region_Code ;
					break;
				}
			}	
			if( region_found == NO )
				Eregion_indx = 5 ;

			secondary_stuct_identifier[i][j] = Eregion_indx ;
		}
	

	}

	//Alloc SS_identifier (init taken care of by its constructor)
	SS_identifier = new struct SecondaryStruct_identifier* [NSTRUCTS] ;
	for( int i = 0 ; i < NSTRUCTS ; i++ ) 
		SS_identifier[i] = new struct SecondaryStruct_identifier [PROT_SIZES[i]] ;
	
	for( int i = 0 ; i < NSTRUCTS ; i++ ) 
		for( int j = 1 ; j < PROT_SIZES[i]-1 ; j++ ) 
		{
			//for the moment just keep track of helices only
			if( secondary_stuct_identifier[i][j] == 0 || secondary_stuct_identifier[i][j] == 6 )
			{
				//check if it is a continuation
				if( SS_identifier[i][j-1].SS_code == secondary_stuct_identifier[i][j] )
				{
					SS_identifier[i][j].SS_code = SS_identifier[i][j-1].SS_code ;
					SS_identifier[i][j].start = SS_identifier[i][j-1].start ;
					SS_identifier[i][j].end = SS_identifier[i][j-1].end ;
				}
				else
				{
					int k ;
					//else check if there exists a stretch of residues of the same type
					for( k = j+1 ; k < PROT_SIZES[i]-1 ; k++ )
						if( secondary_stuct_identifier[i][k] != secondary_stuct_identifier[i][j] )
							break ;
					if ( (k-j) >= 4 ) 
					{
						SS_identifier[i][j].SS_code = secondary_stuct_identifier[i][j] ;
						SS_identifier[i][j].start = j ;
						SS_identifier[i][j].end = k-1 ;
					}
						
				}
			}
			//add Psheets as well
			else if( secondary_stuct_identifier[i][j] == 1 || secondary_stuct_identifier[i][j] == 2 )
			{
				//check if it is a continuation
				if( SS_identifier[i][j-1].SS_code == 1 || SS_identifier[i][j-1].SS_code == 2 )
				{
					
					//SS_identifier[i][j].SS_code = SS_identifier[i][j-1].SS_code ;
					SS_identifier[i][j].SS_code = secondary_stuct_identifier[i][j] ;
					SS_identifier[i][j].start = SS_identifier[i][j-1].start ;
					SS_identifier[i][j].end = SS_identifier[i][j-1].end ;
				}
				else
				{
					int k ;
					//else check if there exists a stretch of residues of the same type
					for( k = j+1 ; k < PROT_SIZES[i]-1 ; k++ )
						if( secondary_stuct_identifier[i][k] != 1 && 
						    secondary_stuct_identifier[i][k] != 2 
						  ) break ;
					if ( (k-j) >= 3 ) 
					{
						SS_identifier[i][j].SS_code = secondary_stuct_identifier[i][j] ;
						SS_identifier[i][j].start = j ;
						SS_identifier[i][j].end = k-1 ;
					}
						
				}
			}
		}
}
/* *******************************************
 * GLOBAL FUNCTION MAIN                      *
 * ******************************************* */
int main( int argc , char *argv[] )
{
	clock_t start = clock(); 
	PRINT_HEADER();
	time(&rundate);
	MUSTANG_CMDLINE_PARSER( argc , argv ) ;
	READ_STRUCTURES( struct_paths ) ;
	CALCULATE_DISTANCE_MATRICES() ;
	/*SSE_RK() ;*/   /* (DEFUALT) Identify secondary structural motifs using Richards and Kundrot 
		        distnace matrix approach                                       */
	//exit(0);
	SSE_RAMA() ; /* Identify secondary structural motifs using ramachandran angles */
	
	PRIMARY_LIBRARY_GENERATION() ; // creates both global and local libs and merges them.
	EXTENDED_LIBRARY_GENERATION() ;
	PROGRESSIVE_ALIGNMENT_USING_EXTENDED_EDGE_WEIGHTS() ;
	SUPERPOSE_ON_CORE( struct_paths );
	if( param_Fhtml_flag) OUTPUT_ALGN_HTML( rundate );
	if( param_Fps_flag) OUTPUT_ALGN_POSTSCRIPT( rundate );
	if( param_Fpir_flag) OUTPUT_ALGN_PIR();
	if( param_Ffasta_flag) OUTPUT_ALGN_FASTA();
	if( param_Fmsf_flag) OUTPUT_ALGN_MSF( rundate );


	// free memory
	DE_ALLOC_DISTMAT();
	DE_ALLOC_GLOBAL_EDGE_WEIGHTS();
	DE_ALLOC_EXTENDED_EDGE_WEIGHTS();
	DE_ALLOC_GLOBAL_LIB();
	//DE_ALLOC_LOCAL_LIB();
	DE_ALLOC_MAXIMAL_FRAGMENT_LIB();
	DE_ALLOC_MERGE_LIB();
	DE_ALLOC_2D( ALGN , NSTRUCTS ) ;
	DE_ALLOC_2D( I_ALGN , NSTRUCTS ) ;
	DE_ALLOC_2D( ALGN_QUALITY , NSTRUCTS ) ;
	DE_ALLOC_1D( STRUCT_PATH ) ;
	DE_ALLOC_2D( struct_names, NSTRUCTS ) ;
	DE_ALLOC_2D( struct_paths, NSTRUCTS ) ;
	DE_ALLOC_2D( PROT , NSTRUCTS ) ;
	DE_ALLOC_1D( PROT_SIZES ) ;
	DE_ALLOC_2D( secondary_stuct_identifier , NSTRUCTS ) ;
	DE_ALLOC_2D( SS_identifier , NSTRUCTS ) ;
	cout <<  "\nAll Done! Thank you for using " << "\033[1m" << "MUSTANG"  << "\033[0m.\n" ;
	clock_t ends = clock();
	cout.setf(ios::fixed) ;
   	cout << "Running Time : " << setprecision(2) << (double) (ends - start) / CLOCKS_PER_SEC << "  s." <<  endl;
    cout << "\033[36;4m.                                -=o0o=-                                  .\033[0m\n" ;
   cout << "\033[36;3m|  Visit \033[0m\033[33;3mhttp://lcb.infotech.monash.edu.au\033[0m\033[36;3m for other supporting programs  |\033[0m\n" ;
   cout << "\033[34;3m|            Some suggested utilities for you to try out:                 |\n";
   cout << "\033[33;3m| SST\033[0m\033[36;3m: A Bayesian approach to assigning secondary structures.             \033[0m\033[33;3m|\033[0m\n";
   cout << "\033[33;3m| MMLigner\033[0m\033[36;3m: Pairwise structural alignment using information & compression.\033[0m\033[33;3m|\033[0m\n";
   cout << "\033[33;3m| Super\033[0m\033[36;3m: Rapidly screen entire PDB for well-superposable fragments.       \033[0m\033[33;3m|\033[0m\n";
   cout << "\033[36;4m.                                -=o0o=-                                  .\033[0m\n" ;


     	return 0;
}

