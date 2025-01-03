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
#include<cmath> 

#include "macros.h"
#include "globals.h"
#include "sse_RK.h"
#include "de_alloc_routines.h"
#include "3D_manip_functions.h"

#define HELIX_MISMATCH_LIMIT      1.0
#define H310EXTN_MISMATCH_LIMIT   1.0
#define SHEET_MISMATCH_LIMIT      1.0

/*                    REFERENCE STRUCTURE DISTANCES AND COORDINATES TAKEN FROM:
   F. M. RICHARDS and C. E. KUNDROT,  PROTEINS:Structure, Function, and Genetics  3:71-84 (1988).  */

/* DATALPHA stores the distance masks for a 50-residue ideal alpha-helix*/
float DATALPHA[] = {  0.00,  3.75,  5.36,  5.02,  6.11,  8.53,  9.75, 10.43, 12.18, 14.09, 
		     15.15, 16.32, 18.19, 19.77, 20.85, 22.33, 24.13, 25.49, 26.70, 28.36, 
		     30.01, 31.27, 32.65, 34.35, 35.84, 37.11, 38.64, 40.30, 41.67, 43.06, 
		     44.64, 46.20, 47.52, 48.97, 50.61, 52.07, 53.41, 54.95, 56.55, 57.93, 
		     59.33, 60.93, 62.45, 63.81, 65.29, 66.89, 68.33, 69.72, 71.27, 72.821 
                   };

/* DATBETA stores the distance masks for a 20-residue ideal strand from a parallel beta-sheet */
float DATBETA[]  = {  0.00,  3.75,  6.47,  9.89, 12.94, 
		     16.28, 19.40, 22.72, 25.87, 29.17, 
		     32.34, 35.62, 38.81, 42.09, 45.28, 
		     48.55, 51.74, 55.01, 58.21, 61.481 
		   };

/* DATTURN stores the triangular mask (3,5) for a sharp beta-turn */
float DATTURN[]  = { 0.00, 0.00, 0.00, 3.70, 0.00,
		     0.00, 5.60, 3.70, 0.00, 5.00,
		     5.10, 3.70, 5.40, 6.90, 6.10  
		   };

/* DATBB stores the distances between one CA atom and the atom in a neighbor strand in a beta-sheet */
float DATBB[]    = {  4.95,  6.20,  8.45, 11.45, 14.55, 
		     17.65, 20.90, 24.15, 27.45, 30.65, 
		     34.00, 37.25, 40.55
		   };


/* DATA310 stores the X,Y,Z coordinates for an alpha helix (residues 1-30) followed by a 3-10 helix 
   (residues 31-50). The axes are colinear with each other and with the X coordinate axis. On execution
   the program computes a rederence distance matrix from this coordinate set.*/
float DATA310[50][3] = {
	{ -0.042, -1.388,  1.834}, {  1.472, -1.616, -1.673}, {  3.013,  1.872, -1.334}, {  4.464,  1.004,  2.097},
	{  5.939, -2.228,  0.679}, {  7.427, -0.325, -2.287}, {  8.917,  2.295,  0.069}, { 10.402, -0.459,  2.268}, 
	{ 11.884, -2.174, -0.812}, { 13.360,  1.146, -2.011}, { 14.844,  2.143, -0.899}, { 16.334, -1.730,  1.525},
	{ 17.815, -1.221, -1.965}, { 19.289,  0.238, -2.309}, { 20.777,  0.548,  2.244}, { 22.266, -2.289,  0.155},
	{ 23.743, -1.898, -1.281}, { 25.222,  2.260,  0.581}, { 26.712, -0.916,  2.104}, { 28.197,  1.792,  1.456},
	{ 29.673,  1.605, -1.704}, { 31.158,  1.456,  1.818}, { 32.647, -1.992,  1.096}, { 34.126, -0.717, -2.193},
	{ 35.607,  2.325, -0.409}, { 37.096,  0.079,  2.304}, { 38.572, -2.247, -0.349}, { 40.053,  0.737, -2.226},
	{ 41.558,  2.124,  1.009}, { 43.051, -1.299,  1.839}, { 44.986, -0.755, -1.696}, { 46.921,  1.843,  0.220},
	{ 48.857, -1.132,  1.471}, { 50.792, -0.684, -1.726}, { 52.728,  1.833,  0.295}, { 54.663, -1.192,  1.424},
	{ 56.598, -0.613, -1.752}, { 58.534,  1.819,  0.370}, { 60.469, -1.249,  1.373}, { 62.404, -0.541, -1.776},
	{ 64.340,  1.803,  0.444}, { 66.275, -1.304,  1.321}, { 68.210, -0.467, -1.797}, { 70.146,  1.783,  0.518},
	{ 72.081, -1.357,  1.266}, { 74.017, -0.393, -1.814}, { 75.952,  1.760,  0.591}, { 77.887, -1.408,  1.210},
	{ 79.823, -0.318, -1.829}, { 81.758,  1.734,  0.662}
} ;

float UTRNGLMAT[50][50] ;
void CALC_UTRNGLMAT_310EXTN() 
{
	for( int i = 0 ; i < 50 ; i++ )
	{
		for( int j = i+1 ; j < 50 ; j++ )
		{
			UTRNGLMAT[i][j] = normAminusB( DATA310[i] , DATA310[j] ) ;
			UTRNGLMAT[j][i] = UTRNGLMAT[i][j]; 
		}
		UTRNGLMAT[i][i] = 0 ;
	}
}

void CHECK_310_EXTN( int ISTRUC , int IRES , int ldcntr ) 
{
	int I1 = 30 - ldcntr  ; int J1 = 31  ;
	int I2 = IRES         ; int J2 = IRES +ldcntr  ;
	
	int cntr = 0 ;
	int terminal_flag = OFF;
	for( int k = 0 ; k+J2 < PROT_SIZES[ISTRUC] ; k++ )
	{
		cntr++ ;
		for( int i = 0 ; i <  ldcntr ; i++ )
		{
			for( int j = i ; j < ldcntr ; j++ )
			{
				if( fabs( distance_matrices[ISTRUC][I2+j][J2+k] - 
							  UTRNGLMAT[I1+j][J1+k] ) > H310EXTN_MISMATCH_LIMIT )
				{
					//cout << "HERE\n" ;
					terminal_flag = ON ; 
					break ;
				}
			}
			if( terminal_flag == ON ) break;
		}
		if( terminal_flag == ON ) break;
	}

	if( cntr > 1 )
	{
		for( int k = 0 ; k < ldcntr+cntr-1 ; k++ )
		{
			SS_identifier[ISTRUC][IRES+k].SS_code =  2 ;
			SS_identifier[ISTRUC][IRES+k].start   = IRES ;
			SS_identifier[ISTRUC][IRES+k].end     = IRES+ldcntr+cntr-2 ;
		}
		IRES += (ldcntr+cntr-2) ; 
	}
}

void Assign_CisPeptides( int ISTRUC )
{
	for( int JRES = 0 ; JRES < PROT_SIZES[ISTRUC] ; JRES++ )
	{
		if( distance_matrices[ISTRUC][JRES][JRES+1] <= 3.1 && 
		    distance_matrices[ISTRUC][JRES][JRES+1] >= 2.7   
		  )
		{
			SS_identifier[ISTRUC][JRES].SS_code =  1 ;
			SS_identifier[ISTRUC][JRES].start   = -9 ;
			SS_identifier[ISTRUC][JRES].end     = -9 ;
		}
	}


}

void Assign_Alpha_310_helices( int ISTRUC )
{
	float *ldmask = new float [ PROT_SIZES[ISTRUC] ] ;
	int ldcntr = 0 ;
	int Hterminal_flag = OFF;
	for( int IRES = 0 ; IRES < PROT_SIZES[ISTRUC] ; IRES++ )
	{
		ldcntr = 0 ;
		Hterminal_flag = OFF;
		for( int JRES = IRES+1 ; JRES < PROT_SIZES[ISTRUC] ; JRES++ ) /* helix atleast 5 residues length */
		{
			ldcntr++ ;
			
			/* get distances in the right hand subcolumn. i.e. the subcolumn along
			   the vetical line (ires,jres)--(jres,jres) in the (upper) triangle of 
			   distances given by (ires,ires), (ires,jres), and (jres,jres).
			*/  
			for( int k = JRES, l = 0 ; k >= IRES ; k-- , l++ )
				ldmask[l] =  distance_matrices[ISTRUC][k][JRES] ;

			//check this distance mask against that of an ideal helix
			for( int k = 0 ; k <= JRES-IRES ; k++ )
			{
				//cout << "\n" << IRES << " "<< JRES<< " | " << fabs( ldmask[k] - DATALPHA[k] ) << " " << HELIX_MISMATCH_LIMIT << "\n" ;
				//sleep(1);
				if( fabs( ldmask[k] - DATALPHA[k] ) > HELIX_MISMATCH_LIMIT )
				{
					//cout << "HERE\n" ;
					Hterminal_flag = ON ; 
					break ;
				}
			}
			if(Hterminal_flag == ON) {
					//cout << "HERE1\n" ;
				break;
			}
		}
		if( (Hterminal_flag == ON || IRES + 1 == PROT_SIZES[IRES])  && ldcntr > 5 )
		{
			CHECK_310_EXTN( ISTRUC , IRES , ldcntr-1 ) ;
			for( int k = 0 ; k < ldcntr-1 ; k++ )
			{
				SS_identifier[ISTRUC][IRES+k].SS_code =  2 ;
				SS_identifier[ISTRUC][IRES+k].start   = IRES ;
				SS_identifier[ISTRUC][IRES+k].end     = IRES + ldcntr -2 ;
			}
			IRES += ldcntr - 2 ; 
		}
	}
	DE_ALLOC_1D( ldmask ) ;

}
/*
void Assign_SharpTurns( int ISTRUC )
{
	for( int JRES = 0 ; JRES < PROT_SIZES[ISTRUC] ; JRES++ )
	{
		int assign_flag = ON ;
		for( int k = 0 ; k < 5 ; k++ )
			if( SS_identifier[ISTRUC][JRES+k].SSE_RK == 2 ) 
			{
				assign_flag = OFF ;
				break;
			}
		if(assign_flag == ON)
		{
			cout << "THINK!!!!"	;
		}
		
	}

}
*/
void Assign_Strand( int ISTRUC )
{
	float *ldmask = new float [ PROT_SIZES[ISTRUC] ] ;
	int ldcntr = 0 ;
	int Bterminal_flag = OFF;
	for( int IRES = 0 ; IRES < PROT_SIZES[ISTRUC] ; IRES++ )
	{
		if( SS_identifier[ISTRUC][IRES].SS_code != -1 ) continue; 
		ldcntr = 0 ;
		Bterminal_flag = OFF;
		for( int JRES = IRES+1 ; JRES < PROT_SIZES[ISTRUC] ; JRES++ ) /* sheet atleast 4 residues length */
		{
			if( SS_identifier[ISTRUC][JRES].SS_code != -1 ) {
					Bterminal_flag = ON ; 
				break; 
			}
			ldcntr++ ;
			
			/* get distances in the right hand subcolumn. i.e. the subcolumn along
			   the vetical line (ires,jres)--(jres,jres) in the (upper) triangle of 
			   distances given by (ires,ires), (ires,jres), and (jres,jres).
			*/  
			for( int k = JRES, l = 0 ; k >= IRES ; k-- , l++ )
				ldmask[l] =  distance_matrices[ISTRUC][k][JRES] ;

			//check this distance mask against that of an ideal helix
			for( int k = 0 ; k <= JRES-IRES ; k++ )
			{
				//cout << "\n" << IRES << " "<< JRES<< " | " << fabs( ldmask[k] - DATALPHA[k] ) << " " << HELIX_MISMATCH_LIMIT << "\n" ;
				//sleep(1);
				if( fabs( ldmask[k] - DATBETA[k] ) > SHEET_MISMATCH_LIMIT )
				{
					//cout << "HERE\n" ;
					Bterminal_flag = ON ; 
					break ;
				}
			}
			if(Bterminal_flag == ON) {
					//cout << "HERE1\n" ;
				break;
			}
		}
		if( Bterminal_flag == ON && ldcntr > 3 )
		{
			for( int k = 0 ; k < ldcntr-1 ; k++ )
			{
				SS_identifier[ISTRUC][IRES+k].SS_code =  4 ;
				SS_identifier[ISTRUC][IRES+k].start   = IRES ;
				SS_identifier[ISTRUC][IRES+k].end     = IRES + ldcntr -2 ;
			}
			IRES += ldcntr - 2 ; 
		}
	}
	DE_ALLOC_1D( ldmask ) ;
}

void SSE_RK()
{
	CALC_UTRNGLMAT_310EXTN() ;
	//Alloc SS_identifier (init taken care of by its constructor)
	SS_identifier = new struct SecondaryStruct_identifier* [NSTRUCTS] ;

	for( int ISTRUC = 0 ; ISTRUC < NSTRUCTS ; ISTRUC++ )
	{
		SS_identifier[ISTRUC] = new struct SecondaryStruct_identifier [PROT_SIZES[ISTRUC]] ;
		Assign_CisPeptides( ISTRUC );
		Assign_Alpha_310_helices( ISTRUC );
		//Assign_SharpTurns( ISTRUC );
		Assign_Strand( ISTRUC );
	}
	///*
	//print 
	for( int i = 0 ; i < NSTRUCTS ; i++ ) 
	{
		//cout << struct_paths[i] << endl ;
		for( int j = 1 ; j < PROT_SIZES[i]-1 ; j++ ) 
		{
			//print debug info
			cout << ios::fixed ;
			cout << setw(4) << j << setw(9) << PROT[i][j].res_num << "(" << PROT[i][j].res_name << ")" ;
			cout << " " << setw(8) << setprecision(1) << PROT[i][j].Phi ;
			cout << " " << setw(8) << setprecision(1) << PROT[i][j].Psi ;
			cout << " " << " ||< " << SS_identifier[i][j].SS_code ;
			cout << "(" << setw(4) << SS_identifier[i][j].start ;
			cout << "--" << setw(4) << SS_identifier[i][j].end << ") >||" ;
			
			//if( secondary_stuct_identifier[i][j] == 0 ) cout << " - HELIX\n" ;
			//else if(secondary_stuct_identifier[i][j] == 1 ) cout << " - SHEET\n" ;
			//else cout << " - NONE\n" ;
			
			switch( SS_identifier[i][j].SS_code )
			{
				case 1: cout << " -  CIS-PEPTIDE\n";
					break ;
				case 2: cout << " -  ALPHA-HELIX\n";
					break ;
				case 3: cout << " -    ESHEET\n";
					break ;
				case 4: cout << " -    STRAND\n";
					break ;
				case 5: cout << " -     NONE\n";
					break ;
				case 6: cout << " -      LHELIX\n";
					break ;
				default: cout << " -    NONE\n";
					break ;
			}
		}
		cout << "\n-----\n\n\n" ;
	}
	//*/


}

