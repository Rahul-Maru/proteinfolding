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
/*
 * DATE:        4  Sept 2006
 * REVISION 1:  19 Oct  2006   (bug-fix when NCORE = 0)
*/
#include <iostream>
using std::cin;         using std::cout;        using std::cerr;
using std::flush;       using std::endl;        using std::ios;
using std::fixed;


#include <iomanip>
using std::setprecision ; using std::setw ;

#include <fstream>
#include <cstring>
using std::ifstream ;
using std::ofstream ;

#include "macros.h"
#include "globals.h"
#include "alloc_routines.h"
#include "de_alloc_routines.h"
#include "init_routines.h"
#include "superpose_on_core.h"
#include "superpose_2.h"
#include "read_structures.h"
#include "multiple_superposition.h"
#include "3D_manip_functions.h"

int **core_columns ;
int *core_columns_2 , NCORE = 0 ;
int **algn_to_seq_hash ;
/*
   This function defines core between each 1st and ith structures (i = 2,NSTRUCS).
   */
void FIND_CORE_COLUMNS_IN_ALGN()  
{
	ALLOC_2D( &core_columns , NSTRUCTS ,  ALGN_LEN ) ;
	INIT_2D( core_columns , NSTRUCTS , ALGN_LEN , 0 ) ;
	
	for( int j = 0 ; j < ALGN_LEN ; j++ )
		if( ALGN[0][j] != '_' )
			core_columns[0][j] = YES ;
			
				
	for( int i = 1 ; i < NSTRUCTS ; i ++)
		for( int j = 0 ; j < ALGN_LEN ; j++ )
		{
			if( ALGN[0][j] != '-' && ALGN[i][j] != '-' )
				core_columns[i][j] = YES ;
		}
	
}
/*
   This function alternatively defines core as all ungapped columns in the final alignment
   */
void FIND_CORE_COLUMNS_IN_ALGN_2()  
{
	ALLOC_1D( &core_columns_2 ,  ALGN_LEN ) ;
	INIT_1D( core_columns_2 , ALGN_LEN , 0 ) ;
	NCORE = 0 ;
	
	for( int j = 0 ; j < ALGN_LEN ; j++ )
	{
		int gap_flag = OFF ;
		for( int i = 0 ; i < NSTRUCTS ; i++ )
		{
			if( ALGN[i][j] == '-' ) 
			{ 
				gap_flag = ON ; 
				break ; 
			}
		}
		if( gap_flag == OFF ) 
		{ 
			core_columns_2[j] = YES ; 
			NCORE++ ; 
		} 
	}
}
void CALC_ALGN_TO_SEQ_HASH() 
{
	ALLOC_2D( &algn_to_seq_hash , NSTRUCTS, ALGN_LEN ) ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		int cntr =0 ;
		for( int j = 0 ; j < ALGN_LEN ; j++ )
			if( ALGN[i][j] != '-' ) algn_to_seq_hash[i][j] = cntr++ ;
	}
}

void PRINT_ROTMATS_RMSDS( int NSTRUCTS, float **CMs, float ***ROTMATS, float **RMSDS) {
	cout << "Writing rotation matrices and RMSDs:" ;
	char temp_name[200] ;
	//create path/name for outputfile
	char output_file_name[200] ;
	strcpy( output_file_name , OUTPUT_FILENAME_PREFIX ) ;
	strcpy(temp_name , output_file_name );
	strcat(temp_name , ".rms_rot" ) ;
	if(!meditate)
	{
		     cout<< "(...writing " 
			 << "\033[31;1m"   
	       	         << temp_name  
		         << "\033[0m" << ')'
			 << setw(20-strlen(temp_name))<< " "  ;
	}
	ofstream out( temp_name , ios::out ) ;
	out << "RMSD matrix (based on multiple superpostion):\n" ;
	out << "     " ;
	for( int j = 0 ; j < NSTRUCTS ; j++ ) {
		out << setw(3) << j+1 << "  " ;
	}
	out << endl ;
	out << "     " ;
	for( int j = 0 ; j < NSTRUCTS ; j++ ) {
		out << "-----" ;
	}
	out << endl ;
	for( int i = 0 ; i < NSTRUCTS ; i++ ) {
			out << setw(3) << i+1 << "|" ;
		for( int j = 0 ; j < NSTRUCTS ; j++ ) {
			if(i==j) out << "  ---" ;
			else { out << setw(5) 
				<< setprecision(1) 
				<< fixed 
				<< RMSDS[i][j];
			}
		}
		out <<  endl ;
	}
	out << endl ; 

	for( int i = 0 ; i < NSTRUCTS ; i++ ) {
		out << "Traslation vector for struct " 
			<< i+1 << " (" << struct_names[i] << ")" 
			<< ":\n" ;  
		for( int j = 0 ; j < 3 ; j++  ) {
			out << setw(8) 
				<< setprecision(3) 
				<< fixed 
				<< CMs[0][j] - CMs[i][j];
		}
		out << endl ;
		out << "Rotation matrix for struct " 
			<< i+1 << " (" << struct_names[i] << ")" 
			<< ":\n" ;  
		for( int j = 0 ; j < 3 ; j++  ) {
			for( int k = 0 ; k < 3 ; k++ ) {
				out << setw(8) 
					<< setprecision(3) 
					<< fixed 
					<< ROTMATS[i][j][k] ;
			}
			out << endl;
		}
		out << endl ;
	}
	out.close() ;
	cout << "[ \033[32;4mOK\033[0m ]\n" ;
}

void PRODUCE_SINGLE_SUPERPOSED_PDB( )
{
	char temp_name[200] ;
	//create path/name for outputfile
	//char output_file_name[200] = "alignment" ;
	char output_file_name[200] ;
	strcpy( output_file_name , OUTPUT_FILENAME_PREFIX ) ;
	strcpy(temp_name , output_file_name );
	strcat(temp_name , ".pdb" ) ;
	if(!meditate)
	{
		if( strlen(temp_name) < 18 ) {
		     cout<< "(...writing " 
			 << "\033[31;1m"   
	       	         << temp_name  
		         << "\033[0m" << ')'
			 << setw(20-strlen(temp_name))<< " "  ;
		}
		else {
			cout << "\n" ;
			cout << " (...writing "  ;
		        cout << "\033[31;1m"   
	       	             << temp_name 
		             << "\033[0m" 
			     << ")"
			     << setw(55-strlen(temp_name))<< " "  ;
		}
			
	}	

	//ofstream transf( "./results/new_sup.pdb" , ios::out ) ;
	ofstream transf( temp_name , ios::out ) ;
	struct Complete_pdb_info *pdb_trav ;
	int cntr = 0 ;
	char chain_id = 'A' ;

	char temp_num[6] ;
	char temp_ins=' ';

	// Fill Headers
	transf.unsetf(ios::left);
	transf.unsetf(ios::right);
	transf << "REMARK  Produced by MUSTANG " << VERSION << "\n" ;
	transf << "REMARK  Authors: A S Konagurthu, J C Whisstock, P J Stuckey, and A M Lesk\n" ;
	transf << "REMARK\n" ;
	transf << "REMARK  Structures: ";
	int temp_cntr = 0 ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		if(++temp_cntr % 5 == 0 )
			transf << "\nREMARK              ";
		//else
		transf << struct_names[i] << " " ;
	}
	transf << "\nREMARK\n" ;
	transf << "REMARK All structures were superposed on the coordinate frame of: " << struct_names[0] ;
	transf << "\nREMARK\n";

	// first the Rigid Coords
	for( int j = 0 ; j < PDB_SIZES[0] ; j++ )
	{
		pdb_trav = COMPLETE_PDBs[0][j] ;
		while( pdb_trav != NULL )
		{
			temp_ins=' ';
			cntr++;
			transf.setf(ios::left) ;
			transf << setw(6) <<"ATOM" ;//1-6
			transf.unsetf(ios::left) ; 
			transf.setf(ios::right) ; 
			transf << setw(5) << cntr ;//7-11
			transf.unsetf(ios::right); 
			transf << " "  ; //12
			transf <<  pdb_trav->atom_name ; //13-17
				
			transf << setw(3) << pdb_trav->residue ; //18-20
			transf << " "  ; //21
			transf << chain_id  ; //22

			for( int k = 0 , l = 0 ; k < (signed)strlen(pdb_trav->residue_num) ; k++ )
			{
				if( pdb_trav->residue_num[k] >= 48 && pdb_trav->residue_num[k] <=57 )
					temp_num[ l++ ] = pdb_trav->residue_num[k] ;
				else
					temp_ins=pdb_trav->residue_num[k];
				
				if( k == (signed)strlen(pdb_trav->residue_num) - 1 )
					temp_num[l]='\0' ;
			}
			
			transf.setf(ios::right) ;
			transf << setw(4) << temp_num ; //23-26
			transf << temp_ins  ; //27
			transf << "   "  ; //28-30
			transf << setprecision(3) ;
			transf.setf(ios::fixed) ;
			transf << setw(8) << pdb_trav->coords[0] ; //31-38
			transf << setw(8) << pdb_trav->coords[1] ; //39-46
			transf << setw(8) << pdb_trav->coords[2] ; //47-54
			transf << setprecision(2) << setw(6) << pdb_trav->occupancy ; //55-60
			transf << setprecision(2) << setw(6) << pdb_trav->B_factor  ; //61-66
			transf << endl;
			transf.unsetf(ios::right) ;
			
			
			pdb_trav = pdb_trav->link ;
		}
		
	}
	transf << "TER                  "  << chain_id << endl ;
	//filling up with rest of the structures
	for( int i = 1 ; i < NSTRUCTS ; i++ )
	{
		cntr = 0 ;
		chain_id++ ;
		for( int j = 0 ; j < PDB_SIZES[i] ; j++ )
		{
			pdb_trav = COMPLETE_PDBs[i][j] ;
			while( pdb_trav != NULL )
			{
				temp_ins=' ';
				cntr++;
				transf.setf(ios::left) ;
			        transf<< setw(6) <<"ATOM" ;//1-6
			       	transf.unsetf(ios::left); 
				transf.setf(ios::right) ; 
				transf << setw(5) << cntr ;//7-11
				transf.unsetf(ios::right); 
				transf << " "  ; //12
				transf <<  pdb_trav->atom_name ; //13-17
				//transf << " "  ; //17
				transf << setw(3) << pdb_trav->residue ; //18-20
				transf << " "  ; //21
				transf << chain_id  ; //22

				for( int k = 0 , l = 0 ; k < (signed)strlen(pdb_trav->residue_num) ; k++ )
				{
					if( pdb_trav->residue_num[k] >= 48 && pdb_trav->residue_num[k] <=57 )
						temp_num[ l++ ] = pdb_trav->residue_num[k] ;
					else
						temp_ins=pdb_trav->residue_num[k];
					
					if( k == (signed)strlen(pdb_trav->residue_num) - 1 )
						temp_num[l]='\0' ;
				}
				
				transf.setf(ios::right) ;
			       	transf << setw(4) << temp_num ;//23-26
				transf.unsetf(ios::right); 
				transf << temp_ins  ; //27
				transf << "   "  ; //28-30
				transf << setprecision(3) ;
				transf.setf(ios::fixed) ;
				transf << setw(8) << pdb_trav->sup_coords[0] ; //31-38
				transf << setw(8) << pdb_trav->sup_coords[1] ; //39-46
				transf << setw(8) << pdb_trav->sup_coords[2] ; //47-54
				transf << setprecision(2) << setw(6) << pdb_trav->occupancy ; //55-60
				transf << setprecision(2) << setw(6) << pdb_trav->B_factor  ; //61-66
				transf << endl;
				transf.unsetf(ios::right) ;

				pdb_trav = pdb_trav->link ;
			}
			
		}
		transf << "TER                  "  << chain_id<< endl ;

	}
	transf << "END\n"<< endl ;
	transf.close();
}
void DE_ALLOC_COMPLETE_PDBs()
{
	struct Complete_pdb_info *pdb_trav , *pdb_trav_nxt ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		for( int j = 0 ; j < PDB_SIZES[i] ; j++  )
		{
			pdb_trav = COMPLETE_PDBs[i][j] ;
			while( pdb_trav != NULL )
			{
				pdb_trav_nxt = pdb_trav->link ;
				delete pdb_trav ;
				pdb_trav = pdb_trav_nxt ;
			}
		}
		delete[] COMPLETE_PDBs[i] ;
	}
	delete[] COMPLETE_PDBs ;
}


void CREATE_ALGN_QUALITY_MASK () 
{
	//alloc memory for alignment quality
	ALLOC_2D( &ALGN_QUALITY , NSTRUCTS , ALGN_LEN ) ;

	//NOTE: ALGN_QUALITY is a boolean mask over the ALGN. 
	// 1 -signifies good quality
	// 0 -signifies bad  quality

	//init all to 1 ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		for( int j = 0 ; j < ALGN_LEN ; j++ )
			ALGN_QUALITY[i][j] = 1 ;

	int nres = 0 ; // # of residues
	int ncpairs ;  // # of correct pairs
	struct Complete_pdb_info *pdb_trav1 , *pdb_trav2 ;
	float diameter ;
	for( int J = 0 ; J < ALGN_LEN ; J++ )
	{
		nres = 0 ;
		
		for( int i = 0 ; i < NSTRUCTS ; i++ ) if ( I_ALGN[i][J] != -99 ) nres++ ;
		if( nres < 2 ) continue ;

		for( int i1 = 0 ; i1 < NSTRUCTS ; i1++ )
		{
			ncpairs = 0 ;
			if( I_ALGN[i1][J] == -99 )  continue ;
			pdb_trav1 = COMPLETE_PDBs[i1][ I_ALGN[i1][J] ] ;
			while( pdb_trav1 != NULL )
			{
				if( !strcmp( pdb_trav1->atom_name , " CA  ") ) break;
				else pdb_trav1 = pdb_trav1->link ;
			}
			if( pdb_trav1 == NULL ) continue ;

			for( int i2 = 0 ; i2 < NSTRUCTS ; i2++ )
			{
				if( i1 == i2 ) continue ;
				if( I_ALGN[i2][J] == -99 ) continue ;
				
				pdb_trav2 = COMPLETE_PDBs[i2][ I_ALGN[i2][J] ] ;
				while( pdb_trav2 != NULL )
				{
					if( !strcmp( pdb_trav2->atom_name , " CA  ") ) break;
					else pdb_trav2 = pdb_trav2->link ;
				}
				if( pdb_trav2 == NULL ) continue ;

				diameter = normAminusB( pdb_trav1->sup_coords , pdb_trav2->sup_coords ) ;

				if( diameter < CA_CA_DIAMETER ) ncpairs++ ;
			}
			float temp;
			if( (temp = (float)ncpairs/(float)(nres-1)) < 0.5 )
				ALGN_QUALITY[i1][J] = 0 ;
		}

	}

}

void SUPERPOSE_ON_CORE( char **paths )
{
	FIND_CORE_COLUMNS_IN_ALGN_2() ;
	
	if( NCORE < 6 ){
		if( param_s_flag && NCORE > 0) {
			cout << "Superposing structures on the core:                                              "<< "[ \033[31mFAILED\033[0m ] " ;
			cout << "\033[41;30mToo few matched columns to superpose on!\033[0m\n";
		}
		else if( param_s_flag && NCORE == 0) {
			cout << "Superposing structures on the core:                                              "<< "[ \033[31mFAILED\033[0m ] " ;
			cout << "\033[41;30mNo matched columns to superpose on!\033[0m\n";
		}
		param_D_flag = 0;
		return ;
	}
	
	CALC_ALGN_TO_SEQ_HASH() ;
	READ_ENTIRE_PDBS( paths ) ;
	if(!meditate || !param_s_flag ) cout << "Superposing structures on the core: " << std::flush ;

	///*
	//optimal multiple simultaneous superposition
	float ***Coords , **CMs , ***ROTMATS, **RMSDS;
	//alloc the above
	Coords = new float** [NSTRUCTS] ;
	CMs = new float* [NSTRUCTS] ;
	ROTMATS = new float** [NSTRUCTS] ;
	RMSDS = new float* [NSTRUCTS] ;
	
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		Coords[i] = new float* [NCORE] ;
		   CMs[i] = new float [3] ;
		   ROTMATS[i] = new float* [3] ;
		 RMSDS[i] = new float [NSTRUCTS] ;
	}
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		for( int j = 0 ; j < NCORE ; j++ )
			Coords[i][j] = new float [3] ;

		for( int j = 0 ; j < 3 ; j++ )
		   	ROTMATS[i][j] = new float [3] ;

	}
	

	//gather core coords from all structures
	int Size = 0 ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		Size = 0 ;
		for( int j = 0 ; j < ALGN_LEN ; j++ )
		{
			if( core_columns_2[j] == ON )
			{
				for( int m = 0 ; m < 3 ; m++ )
					Coords[i][Size][m] = PROT[i][algn_to_seq_hash[i][j]].CA_coords[m] ;
				Size++;
			}
		}
	}

	// Call to the Multiple simultaneous superposition routine
	M_SUPERPOSE( NSTRUCTS , NCORE , Coords , CMs , ROTMATS, RMSDS ); 

	//superpose all structures using ROTMATS
	//void rotate_vector_2(float [] , const float [][3] , const int ) ;
	void rotate_vector_2(float [] , float** , const int ) ;
	for( int j = 1 ; j < NSTRUCTS ; j++ ) 
	{
		float temp_coords[3] ;
		struct Complete_pdb_info *pdb_trav ;
		for( int k = 0 ; k < PDB_SIZES[j] ; k ++ )
		{
			pdb_trav = COMPLETE_PDBs[j][k] ;
			while( pdb_trav != NULL )
			{
				for( int l = 0 ; l < 3 ; l++ )
					temp_coords[l] = pdb_trav->coords[l] ;
				
				for( int l = 0 ; l < 3 ; l++ )
					temp_coords[l] = temp_coords[l] - CMs[j][l] ;

				rotate_vector_2( temp_coords , ROTMATS[j] , 3 ) ;
		
				for( int l = 0 ; l < 3 ; l++ )
					temp_coords[l] = temp_coords[l] + CMs[0][l] ;
				for( int l = 0 ; l < 3 ; l++ )
					pdb_trav->sup_coords[l] = temp_coords[l] ;

					//for( int l = 0 ; l <  3 ; l++ )
					//	cout << setw(9) << pdb_trav->coords[l] ;
					//cout << "|" ;
					//for( int l = 0 ; l <  3 ; l++ )
					//	cout << setw(9) << pdb_trav->sup_coords[l] ;
					//cout << "\n" ;
				pdb_trav = pdb_trav->link;
			}
		}		
	}
	
	//since first structure is rigid, simply copy all pdb_trav->coords into pdb_trav->sup_coords for this structure
	{ // this is just a plain block used to limit the scope of local variables defined in this.
		//float temp_coords[3] ;
		struct Complete_pdb_info *pdb_trav ;
		for( int k = 0 ; k < PDB_SIZES[0] ; k ++ )
		{
			pdb_trav = COMPLETE_PDBs[0][k] ;
			while( pdb_trav != NULL )
			{
				for( int l = 0 ; l < 3 ; l++ )
					pdb_trav->sup_coords[l] = pdb_trav->coords[l] ;

				pdb_trav = pdb_trav->link;
			}
		}		
	}

	CREATE_ALGN_QUALITY_MASK () ; 

	if( param_s_flag) PRODUCE_SINGLE_SUPERPOSED_PDB() ;

	//DE_ALLOC_2D(core_columns , NSTRUCTS) ;
	DE_ALLOC_3D(Coords , NSTRUCTS , NCORE ) ;
	DE_ALLOC_1D(core_columns_2 ) ;
	DE_ALLOC_COMPLETE_PDBs() ;
	DE_ALLOC_1D(PDB_SIZES) ;
	DE_ALLOC_2D(algn_to_seq_hash,NSTRUCTS) ;
	if(!meditate || param_s_flag )
		cout << "[ \033[32;4mOK\033[0m ]\n" ;
	if( param_r_flag ) PRINT_ROTMATS_RMSDS( NSTRUCTS,CMs, ROTMATS,RMSDS ) ;
	DE_ALLOC_3D(ROTMATS , NSTRUCTS , 3 ) ;
	DE_ALLOC_2D(RMSDS , NSTRUCTS ) ;
	DE_ALLOC_2D(CMs , NSTRUCTS ) ;
}

void rotate_vector_2( float v[] , float **R , const int S )
{
	float *s = new float [S] ;
	for( int i = 0 ; i < S ; i++ )
	{
		s[i]=v[i] ;
		v[i]=0.0 ;
	}
	for( int i = 0 ; i < S ; i++ )
		for( int j = 0 ; j < S ; j++ )
			v[i] = v[i] + s[j] * R[i][j] ;

	DE_ALLOC_1D(s) ;
}

