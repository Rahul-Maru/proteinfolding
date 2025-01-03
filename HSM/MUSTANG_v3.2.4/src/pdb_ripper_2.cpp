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


#include <iomanip>
using std::setprecision ; using std::setw ; 

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using std::ifstream ;

#include "macros.h"
#include "globals.h"
#include "pdb_ripper.h"
#include "alloc_routines.h"
#include "init_routines.h"
int CHECK_NCHAINS_IN_PDB( char *fname )
{
	int NCHAINS = 0 ;
	int chain_flag = OFF ;
	char buffer[90];
	ifstream pdb( fname , ios::in ) ;
	if( !pdb )
	{
       		cerr << "Pdb File Error: Structure " << fname << " does not exist in the path specified.\n" ; 
		pdb.close();
       		exit(1) ;
       	}
	while ( !pdb.eof() )
	{
		pdb.getline( buffer , 90 ) ;
		if( pdb.eof() || 
		    (buffer[0] == 'E' && buffer[1] == 'N' && buffer[2] == 'D' && buffer[3] == ' ') 
		  ) break ;

		if( buffer[0] == 'A' && buffer[1] == 'T' && buffer[2] == 'O' && buffer[3] == 'M' && chain_flag == OFF )
		{
			chain_flag = ON ;
		}
			
		if( chain_flag == ON && buffer[0] == 'T' && buffer[1] == 'E' && 
		    buffer[2] == 'R' && buffer[3] == ' ' /*&& buffer[ 21 ] == chain_id */
		  )
		{
			chain_flag = OFF ;
			NCHAINS++ ;
		}
	}
	return NCHAINS ;
}
int CHECK_NRESIDUES_IN_PDB( char *fname )
{
	int NRESIDUES = 0 ;
	char buffer[90];
	ifstream pdb( fname , ios::in ) ;
	if( !pdb )
	{
       		cerr << "Pdb File Error: Structure " << fname << " does not exist in the path specified.\n" ; 
		pdb.close();
       		exit(1) ;
       	}
	char prev_res_num[8] = "" ;
	char curr_res_num[8] = "" ;
	while ( !pdb.eof() )
	{
		pdb.getline( buffer , 90 ) ;
		
		
		if( pdb.eof() || 
		    (buffer[0] == 'E' && buffer[1] == 'N' && buffer[2] == 'D' && buffer[3] == ' ' )
		  ) break ;

		if( buffer[0] == 'T' && buffer[1] == 'E' && buffer[2] == 'R' && buffer[3] == ' ' )
			break ;
		
		if( buffer[0] == 'A' && buffer[ 1 ] == 'T' && buffer[ 2 ] == 'O' && buffer[ 3 ] == 'M' )
			if( ( buffer[13] == 'N' && buffer[14] == ' ' ) ||
			    ( buffer[13] == 'C' && buffer[14] == 'A' ) ||
			    ( buffer[13] == 'C' && buffer[14] == ' ' ) ||
			    ( buffer[13] == 'O' && buffer[14] == ' ' ) )
			{
				for( int i = 22 , j = 0 ; i < 28 ; i++ , j++ ) curr_res_num[j] = buffer[i] ;
				curr_res_num[6] = '\0' ;
				if( strcmp(curr_res_num, prev_res_num) ) NRESIDUES++ ;
			}
		strcpy( prev_res_num , curr_res_num ) ;
		
	}
	return NRESIDUES ;
}

void PARSE_PDB_STRUCTURE( char *fname , int Sindx ) 
{
	if(!meditate) 
		cout << "Parsing the PDB file: \033[35m" << setw(15) << struct_names[Sindx] << "\033[0m";
	
	PROT_SIZES[Sindx] = 0 ;
	//char strbuffer[4][90];
	char buffer[90] ;
	char tempstr[20] = "" ;
	char temp_name[4] ;
	char temp_num[10] ;
	float temp_coord[3] ;
	float temp_occupancy ;
	int index = 0 ;
	//char prev_num[9] = "-1";
	int atom_handle = 0 ;

	ifstream pdb( fname , ios::in ) ;
	if( !pdb )
	{
		cerr << "Pdb File Error: Structure " << fname << " does not exist in the path specified.\n" ; 
		pdb.close();
		exit(1) ;
	}

	while ( !pdb.eof() )
	{
		pdb.getline( buffer , 90 ) ;
		if ( (buffer[0] == 'T' && buffer[1] == 'E' && buffer[2] == 'R')  ) break ;
		if ( (buffer[0] != 'A' || buffer[1] != 'T' || buffer[2] != 'O' || buffer[3] != 'M')  ) continue ;
		else if( 
				( buffer[13] == 'N' && buffer[14] == ' ' ) || 
				( buffer[13] == 'C' && buffer[14] == 'A' ) || 
				( buffer[13] == 'C' && buffer[14] == ' ' ) || 
				( buffer[13] == 'O' && buffer[14] == ' ' )  
		       )
		{
			     if( buffer[13] == 'N' && buffer[14] == ' ' ) atom_handle = 0 ;
			else if( buffer[13] == 'C' && buffer[14] == 'A' ) atom_handle = 1 ;
			else if( buffer[13] == 'C' && buffer[14] == ' ' ) atom_handle = 2 ;
			else if( buffer[13] == 'O' && buffer[14] == ' ' ) atom_handle = 3 ;


			
			for ( int i = 17 ,  j = 0 ; i < 20 ; i++ )
			{
				tempstr[ j++ ] = buffer[ i ] ;
				if ( i == (20 - 1) )
					tempstr[ j ] = '\0' ;
			}
			sscanf( tempstr , "%s" , temp_name ) ;
			
			for ( int i = 22 ,  j = 0 ; i < 28 ; i++ )
			{
				tempstr[ j++ ] = buffer[ i ] ;
				if( i == (28 - 1) )
					tempstr[ j ] = '\0' ;
			}
			//sscanf( tempstr , "%d" , &temp_num ) ;
			strcpy( temp_num , tempstr) ;
			
			for ( int i = 54 ,  j = 0 ; i < 60 ; i++ )
			{
				tempstr[ j++ ] = buffer[ i ] ;
				if( i == (60 - 1) )
					tempstr[ j ] = '\0' ;
			}
			//sscanf( tempstr , "%d" , &temp_num ) ;
			sscanf( tempstr , "%f" , &temp_occupancy ) ;


			//X coord of  CA 
			for ( int i = 30 ,  j = 0 ; i < 38 ; i++ )
			{
				tempstr[ j++ ] = buffer[ i ] ;
				if( i == (38 - 1) )
					tempstr[ j ] = '\0' ;
			}
			sscanf( tempstr , "%f" , &temp_coord[0] ) ;

			//Y coord of  CA
			for ( int i = 38 ,  j = 0 ; i < 46 ; i++ )
			{
				tempstr[ j++ ] = buffer[ i ] ;
				if( i == (46 - 1) )
					tempstr[ j ] = '\0' ;
			}
			sscanf( tempstr , "%f" , &temp_coord[1] ) ;

			//Z coord of  CA
			for ( int i = 46 ,  j = 0 ; i < 54 ; i++ )
			{
				tempstr[ j++ ] = buffer[ i ] ;
				if( i == (54 - 1) )
					tempstr[ j ] = '\0' ;
			}
			sscanf( tempstr , "%f" , &temp_coord[2] ) ;
		
			if( index != 0 && (!strcmp( PROT[Sindx][index-1].res_num , temp_num )) )
			{
				if( temp_occupancy > PROT[Sindx][index-1].backbone_occupancies[atom_handle]+float_epsilon )
				{
					if( atom_handle == 0 )
					{
						for( int i = 0 ; i < 3 ; i++ )
							PROT[Sindx][index-1].N_coords[i]  = temp_coord[i] ;
					}
					else if( atom_handle == 1 )
					{
						for( int i = 0 ; i < 3 ; i++ )
							PROT[Sindx][index-1].CA_coords[i]  = temp_coord[i] ;
					}	
					else if( atom_handle == 2 )
					{
						for( int i = 0 ; i < 3 ; i++ )
							PROT[Sindx][index-1].C_coords[i]  = temp_coord[i] ;
					}
					else if( atom_handle == 3 )
					{
						for( int i = 0 ; i < 3 ; i++ )
							PROT[Sindx][index-1].O_coords[i]  = temp_coord[i] ;
					}
					PROT[Sindx][index-1].backbone_occupancies[atom_handle] = temp_occupancy ;
				}
			}
			else
			{

				//fillin in the residue number(with insertion code) as it appears in PDB
				strcpy( PROT[Sindx][index].res_num , temp_num ) ;
				PROT[Sindx][ index ].res_indx = index ;
			
				// associating single letter code for this residue
				for( int i = 0 ; i < 23 ; i++ )
				{
					if( !strcmp(res_names[i].three_letter_code , temp_name) )
					{
						PROT[Sindx][index].res_name = res_names[i].single_letter_code ;
						break;
					}
				}
				if( atom_handle == 0 )
				{
					for( int i = 0 ; i < 3 ; i++ )
						PROT[Sindx][index].N_coords[i]  = temp_coord[i] ;
				}
			        else if( atom_handle == 1 )
				{
					for( int i = 0 ; i < 3 ; i++ )
						PROT[Sindx][index].CA_coords[i]  = temp_coord[i] ;
				}	
				else if( atom_handle == 2 )
				{
					for( int i = 0 ; i < 3 ; i++ )
						PROT[Sindx][index].C_coords[i]  = temp_coord[i] ;
				}
				else if( atom_handle == 3 )
				{
					for( int i = 0 ; i < 3 ; i++ )
						PROT[Sindx][index].O_coords[i]  = temp_coord[i] ;
				}
				index++ ;
				PROT_SIZES[Sindx]++ ;
				//cout << index << endl ;
			}
		}
	}

		
	if(!meditate)
	{
		cout << "  (...contains \033[1m" << setw(4) << PROT_SIZES[Sindx] << "\033[0m residues.)";
		cout << "  [ \033[32;4mOK\033[0m ]\n" ;

	}
	
}

/*..........................................................................*/
void PARSE_ENTIRE_PDB_STRUCTURE( char *fname , int Sindx ) 
{
	char strbuffer[90];
	char tempstr[20] = "" ;
	char temp_atom_name[6] ;
	char temp_res_name[4] ;
	char temp_res_num[6] ;
	float temp_coord[3] ;
	float temp_occupancy;
	float temp_bfactor;
	int index = 0 ;
	char prev_res_name[6]="";
	struct Complete_pdb_info *pdb_trav = NULL ;
	ifstream pdb( fname , ios::in ) ;
	if ( !pdb )
	{
		cerr << "\nError in opening the PDB:" <<  fname << endl;
		exit( 1 ) ;
	}
	while ( !pdb.eof() )
	{
		pdb.getline( strbuffer , 90 ) ;
		//cout << (int)strbuffer[0] << endl ;
		if ( strbuffer[ 0 ] == 'T' && strbuffer[ 1 ] == 'E' && strbuffer[ 2 ] == 'R' )
			break ;
		if ( strbuffer[ 0 ] != 'A' || strbuffer[ 1 ] != 'T' || strbuffer[ 2 ] != 'O' || strbuffer[ 3 ] != 'M' )
			continue;
		
		if( index == 0 && COMPLETE_PDBs[Sindx][index] == NULL )
		{
			COMPLETE_PDBs[Sindx][index] =  new struct Complete_pdb_info ;
			pdb_trav = COMPLETE_PDBs[Sindx][index] ;
			pdb_trav -> link =  NULL ;
			for ( int i = 22 ,  j = 0 ; i < 27 ; i++ )
			{
				prev_res_name[ j++ ] = strbuffer[ i ] ;
				if( i == (27 - 1) )
					prev_res_name[ j ] = '\0' ;
			}
			PDB_SIZES[Sindx]++ ;
		}
		else
		{
			for ( int i = 22 ,  j = 0 ; i < 27 ; i++ )
			{
				tempstr[ j++ ] = strbuffer[ i ] ;
				if( i == (27 - 1) )
					tempstr[ j ] = '\0' ;
			}
			if( strcmp( prev_res_name , tempstr ) == 0 ) // same residue
			{
				pdb_trav -> link = new struct Complete_pdb_info ;
				pdb_trav = pdb_trav -> link ;
				pdb_trav -> link =  NULL ;
			}
			else // not same
			{
				index++ ;
				PDB_SIZES[Sindx]++ ;
				//if( index > PROT_SIZES[Sindx] ) // HERE!!!use if statement on PROT_SIZE to execute a breakpoint
				//	cout << index << " " << strbuffer << endl ;
				//cout << index << endl ;
				COMPLETE_PDBs[Sindx][index] =  new struct Complete_pdb_info ;
				pdb_trav = COMPLETE_PDBs[Sindx][index] ;
				pdb_trav -> link =  NULL ;
				for ( int i = 22 ,  j = 0 ; i < 27 ; i++ )
				{
					prev_res_name[ j++ ] = strbuffer[ i ] ;
					if( i == (27 - 1) )
						prev_res_name[ j ] = '\0' ;
				}

			}

		}
		// fill in pdb_info into pdb_trav

		//atom_name
		for ( int i = 12 ,  j = 0 ; i < 17 ; i++ )
		{
			tempstr[ j++ ] = strbuffer[ i ] ;
			if ( i == (17 - 1) )
				tempstr[ j ] = '\0' ;
		}
		strcpy(temp_atom_name, tempstr ) ;
			
		//residue	
		for ( int i = 17 ,  j = 0 ; i < 20 ; i++ )
		{
			tempstr[ j++ ] = strbuffer[ i ] ;
			if ( i == (20 - 1) )
				tempstr[ j ] = '\0' ;
		}
		sscanf( tempstr , "%s" , temp_res_name ) ;
		//residue_num
		
		for ( int i = 22 ,  j = 0 ; i < 27 ; i++ )
		{
			tempstr[ j++ ] = strbuffer[ i ] ;
			if( i == (27 - 1) )
				tempstr[ j ] = '\0' ;
		}
		sscanf( tempstr , "%s" , temp_res_num ) ;
		//X coord of CA
		for ( int i = 30 ,  j = 0 ; i < 38 ; i++ )
		{
			tempstr[ j++ ] = strbuffer[ i ] ;
			if( i == (38 - 1) )
				tempstr[ j ] = '\0' ;
		}
		sscanf( tempstr , "%f" , &temp_coord[0] ) ;

		//Y coord of CA
		for ( int i = 38 ,  j = 0 ; i < 46 ; i++ )
		{
			tempstr[ j++ ] = strbuffer[ i ] ;
			if( i == (46 - 1) )
				tempstr[ j ] = '\0' ;
		}
		sscanf( tempstr , "%f" , &temp_coord[1] ) ;

		//Z coord of CA
		for ( int i = 46 ,  j = 0 ; i < 54 ; i++ )
		{
			tempstr[ j++ ] = strbuffer[ i ] ;
			if( i == (54 - 1) )
				tempstr[ j ] = '\0' ;
		}
		sscanf( tempstr , "%f" , &temp_coord[2] ) ;
		//occupancy
		for ( int i = 54 ,  j = 0 ; i < 60 ; i++ )
		{
			tempstr[ j++ ] = strbuffer[ i ] ;
			if( i == (60 - 1) )
				tempstr[ j ] = '\0' ;
		}
		sscanf( tempstr , "%f" , &temp_occupancy ) ;
		//B-factor
		for ( int i = 60 ,  j = 0 ; i < 66 ; i++ )
		{
			tempstr[ j++ ] = strbuffer[ i ] ;
			if( i == (66 - 1) )
				tempstr[ j ] = '\0' ;
		}
		sscanf( tempstr , "%f" , &temp_bfactor ) ;


		//finally fillin the pdb_info datastructs
		strcpy( pdb_trav -> atom_name , temp_atom_name );
		strcpy( pdb_trav -> residue , temp_res_name );
		strcpy( pdb_trav -> residue_num , temp_res_num );

		for( int i = 0 ; i < 3 ; i++ )
			pdb_trav -> coords[i] = temp_coord[i] ;
		
		pdb_trav -> occupancy = temp_occupancy ;
		pdb_trav -> B_factor = temp_bfactor ;
	}
	pdb.close() ;
}
