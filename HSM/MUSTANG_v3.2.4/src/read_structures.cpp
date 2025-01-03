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
using std::cin;         using std::cout;        using std::cerr;
using std::flush;       using std::endl;        using std::ios;
using std::fixed;


#include <iomanip>
using std::setprecision ; using std::setw ;

#include <fstream>
#include <cstdlib>
#include <cstring>
using std::ifstream ;
using std::ofstream ;


#include "macros.h"
#include "globals.h"
#include "read_structures.h"
#include "pdb_ripper.h"
#include "alloc_routines.h"
#include "de_alloc_routines.h"
#include "init_routines.h"
#include "3D_manip_functions.h"

void CALCULATE_SEQUENCES();
void CALCULATE_BACKBONE_DIHEDRALS(); 

void READ_STRUCTURES( char **paths )
{
	//ALLOC
	PROT = new struct pdb_info* [NSTRUCTS] ;
	PROT_SIZES = new int [NSTRUCTS] ;
	for( int Sindx = 0 ; Sindx < NSTRUCTS ; Sindx++ )
	{
		int nchains = CHECK_NCHAINS_IN_PDB( paths[Sindx] ) ;
		if( nchains > 1 )
		{
			cerr << "Can't Proceed. The structure: " << paths[Sindx] << " contains more than one chain.\n" ;
			cerr << "At this moment this program requires the input pdb structures to " ;
			cerr << "contain just one chain.\n\n Hint! Edit the structure file by retaining " ;
		        cerr << "only the chain of interest\n\n" ;
			cerr << "If you are interested in more that one chain ";
			cerr << "edit the input file such that you append all the interested chains together " ;
			cerr << "by eliminating the TER lines between them.\n\n";
			cerr << "Warning: these appended chains will be treated as contiguous.\n" ;

			for( int i = 0 ; i < Sindx ; i++ ) delete[] PROT[i] ;
			delete[] PROT ;
			exit(1) ;
		}
		int nresidues = CHECK_NRESIDUES_IN_PDB( paths[Sindx] ) ;
		//cerr << endl <<  nresidues ; 
		PROT[Sindx] = new struct pdb_info [nresidues+1] ;
		PARSE_PDB_STRUCTURE( paths[Sindx] , Sindx ) ;
		//cerr << " " << PROT_SIZES[Sindx] << endl ;
		if( PROT_SIZES[Sindx] == 0 ) 
		{
			cerr << "Can't Proceed. One or more input structures is unreabadle\n" ;
			exit(1) ;
		}
		if( PROT_SIZES[Sindx] != nresidues )
		{
			cerr << nresidues << " " << Sindx << " " <<  PROT_SIZES[Sindx] << endl ;
			exit(0) ;
		}
	}
	
	
	//CALCULATE_SEQUENCES();
	CALCULATE_BACKBONE_DIHEDRALS() ;
	//if( meditate == FALSE ) cout << "[ \033[32;1m" << "OK \033[0m]\n" ; 
}

/*
void CALCULATE_SEQUENCES()
{
	SEQS = new char* [NSTRUCTS] ;
	SEQ_LENS = new int [NSTRUCTS] ;
	
	for( int Sindx = 0 ; Sindx <  NSTRUCTS ; Sindx++ )
	{
		SEQS[Sindx] = new char [PROT_SIZES[Sindx]+1] ;
		for( int i = 0 ; i < PROT_SIZES[Sindx] ; i++ )
			SEQS[Sindx][i] = PROT[Sindx][i].res_name ;
		SEQS[Sindx][ PROT_SIZES[Sindx] ] = '\0' ;
		SEQ_LENS[Sindx] = PROT_SIZES[Sindx] ;
	}
}
*/

void CALCULATE_BACKBONE_DIHEDRALS() 
{
	for( int Sindx = 0 ; Sindx <  NSTRUCTS ; Sindx++ )
	{
		for ( int res = 0 ; res < PROT_SIZES[Sindx] ; res++ )
		{
		     if( res != 0 && res != PROT_SIZES[Sindx] - 1 )
		     {
			//Omega
			compute_dihedral_angle( PROT[Sindx][res-1].CA_coords , PROT[Sindx][res-1].C_coords   ,
					        PROT[Sindx][res].N_coords, PROT[Sindx][res].CA_coords, 
						PROT[Sindx][res].Omega
					      );
			//Phi
			compute_dihedral_angle( PROT[Sindx][res-1].C_coords, PROT[Sindx][res].N_coords   , 
					        PROT[Sindx][res].CA_coords , PROT[Sindx][res].C_coords   , 
						PROT[Sindx][res].Phi
					      );
			//Psi
			compute_dihedral_angle( PROT[Sindx][res].N_coords  , PROT[Sindx][res].CA_coords  ,  
					        PROT[Sindx][res].C_coords  , PROT[Sindx][res+1].N_coords ,
						PROT[Sindx][res].Psi
					      );
		     }
		     else
		     {
			     if(res==0)
				//Psi
				compute_dihedral_angle( PROT[Sindx][res].N_coords  , PROT[Sindx][res].CA_coords  ,  
					                PROT[Sindx][res].C_coords  , PROT[Sindx][res+1].N_coords , 
							PROT[Sindx][res].Psi
						      );
			     else
			     {
				//Phi
				compute_dihedral_angle( PROT[Sindx][res-1].C_coords, PROT[Sindx][res].N_coords   , 
					                PROT[Sindx][res].CA_coords , PROT[Sindx][res].C_coords   , 
							PROT[Sindx][res].Phi
						      );
				//Omega
				compute_dihedral_angle( PROT[Sindx][res-1].CA_coords , PROT[Sindx][res-1].C_coords   ,
					       		PROT[Sindx][res].N_coords, PROT[Sindx][res].CA_coords, 
							PROT[Sindx][res].Omega
						      );
			     }

		     }
		}
	}

}
//*/

/*.......................................................................*/
///*
void READ_ENTIRE_PDBS( char **paths ) 
{
	const int ADD_OFFSET  = 100 ; // a crazy dumb hack around the incorrigible pdb format!
	
	//Alloc space for PDB_SIZES
	PDB_SIZES = new int [NSTRUCTS] ;
	INIT_1D( PDB_SIZES , NSTRUCTS,  0 ) ;
	
	//Alloc space for COMPLETE_PDB datastructure
	COMPLETE_PDBs = new struct Complete_pdb_info** [NSTRUCTS] ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		COMPLETE_PDBs[i] = new struct Complete_pdb_info* [PROT_SIZES[i]+ADD_OFFSET] ;
	
	//initialize the linked list ptrs to NULL
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		for( int j = 0 ; j < PROT_SIZES[i] ; j++ )
			COMPLETE_PDBs[i][j] = NULL ;
	
	//parse entire pdb
	for( int Sindx = 0 ; Sindx < NSTRUCTS ; Sindx++ )
	{
		PARSE_ENTIRE_PDB_STRUCTURE( paths[Sindx] , Sindx ) ;
	}
}
//*/
