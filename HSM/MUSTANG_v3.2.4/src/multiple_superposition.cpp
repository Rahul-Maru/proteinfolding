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
/* Routine for Multiple Simultaneous Superposition of structures minimizing 
 * sum-of-pairs residuals. The original paper explaining this algorithm is:
 *      R. Diamond
 *	ON THE MULTIPLE SIMULTANEOUS SUPERPOSITION OF MOLECULAR STRUCTURES 
 *	BY RIDY BODY TRANSFORMATIONS.
 *	PROTEIN SCIENCE
 *	1992 1: 1279--1287.
 *   
 * Needs a multiple alignment (...and of course the coordinates of 
 * the involved structures).
 *	
 * Created: 24 May 2005.
 *	
 * Author:  Arthur M Lesk (Original fortran code)  
 *	    Arun S Konagurthu(C/C++ translated version).
 */   
#include<iostream>
using std::cin; 	using std::cout; 	using std::cerr; 
using std::flush; 	using std::endl; 	using std::ios;
using std::fixed;


#include<iomanip>
using std::setprecision ; using std::setw ; 

#include <cstdlib>

#include "jacobi.h"
#include "alloc_routines.h"
#include "de_alloc_routines.h"
#include <cmath>
#define THRESH 0.0001

void MOVE_TO_CENTER_OF_GRAVITY( int NSTRUC , int NRES , float ***R , float **CMs )
{
	float RECNRS = 1.0/(float)(NRES) ;
	float XBAR , YBAR , ZBAR ;
		
	for( int ISTRUC = 0 ; ISTRUC < NSTRUC ; ISTRUC++ )
	{
		XBAR = 0.0 ;
		YBAR = 0.0 ;
		ZBAR = 0.0 ;
		
		for( int I = 0 ; I < NRES ; I++ )
		{
			XBAR = XBAR + R[ISTRUC][I][0] ;
			YBAR = YBAR + R[ISTRUC][I][1] ;
			ZBAR = ZBAR + R[ISTRUC][I][2] ;
		}
		
		XBAR = RECNRS*XBAR ; 
		YBAR = RECNRS*YBAR ;
		ZBAR = RECNRS*ZBAR ;

		CMs[ISTRUC][0] = XBAR ;	
		CMs[ISTRUC][1] = YBAR ;	
		CMs[ISTRUC][2] = ZBAR ;	
		
		for( int I = 0 ; I < NRES ; I++ )
		{
			R[ISTRUC][I][0] = R[ISTRUC][I][0] - XBAR ;
			R[ISTRUC][I][1] = R[ISTRUC][I][1] - YBAR ;
			R[ISTRUC][I][2] = R[ISTRUC][I][2] - ZBAR ;
		}
	}
}

void EVV4X4( double PP[][4] , double **RHO , int ISTRUC ) 
{
	float EVECS[4][4], EVALUS[4] ;
	JACOBI_ROTATIONS( PP , EVALUS , EVECS ) ;
	
   	RHO[ISTRUC][0] = EVECS[0][0] ;
   	RHO[ISTRUC][1] = EVECS[1][0] ;
   	RHO[ISTRUC][2] = EVECS[2][0] ;
   	RHO[ISTRUC][3] = EVECS[3][0] ;

	//print eigen vectors
			if(0)
			{
				cout << std::setprecision(4) ;
				cout << "\t Eigenvalues of the Quaternion matrix:" ;
				for( int i = 0 ; i < 4 ; i++ )
					cout << std::setw(17) << EVALUS[i] ;
				cout << endl ;
				cout << "\t Eigenvectors of the Quaternion Matrix:\n" ;
				for ( int i = 0 ; i < 4 ; i++ )
				{
					cout<<"\t" ;
					for( int j = 0 ; j < 4 ; j++ )
						cout << std::setw(12) << EVECS[i][j] ;
					cout << "\n" ;
				}
						
				//for( int i = 0 ; i < JACOBI_SIZE ; i++ )
				//	cout << std::setw(7) << eigen_values[i] ;
				//cout << endl ;
			}

}

void MTRXYZ(float CURVUE[3][3] , float RXV , float RYV, float RZV )
{

/*
C     DECOMPOSES CURVUE INTO EQUIVALENT ROTATION ANGLES
C
C   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
C
C
*/
#define RTOD    57.295776367 
	double THRESH_2 =  1.0E-7  ;
	int tempnum1 , tempnum2 ;
	float ALTRX, ALTRY, ALTRZ , A31;
	double RXMRZ , RXPRZ , RX , RY, RZ , ELTMAX , X ;
	RXV = 0.0 ;
	RYV = 0.0 ;
	RZV = 0.0 ;

/*     FIRST ORTHONORMALIZE */

	A31 = CURVUE[0][2] ;
	if( fabs(fabs(A31) - 1.0) < THRESH_2 ) 
		goto label_10 ;

	if(A31 - 1.0) {
		goto label_1 ;
		goto label_10 ;
		goto label_2 ;
	}

label_1:
    	RXMRZ = atan2( CURVUE[2][1]-CURVUE[1][0],CURVUE[1][1]+CURVUE[2][0] ) ;
	goto label_3 ;

label_2:
    	RXMRZ = atan2(CURVUE[1][0]-CURVUE[2][1],-CURVUE[1][1]-CURVUE[2][0]) ;

label_3:
    	if(A31 + 1.0) 
    	{
		goto label_4  ;
	    	goto label_10 ;
	    	goto label_5  ;
    	}

label_4:
    	RXPRZ = atan2(-CURVUE[2][1]-CURVUE[1][0],-CURVUE[1][1]+CURVUE[2][0]) ;
        goto label_6  ;

label_5:	
    	RXPRZ = atan2(CURVUE[2][1]+CURVUE[1][0],CURVUE[1][1]-CURVUE[2][0]) ;

label_6:
    	RX = 0.5*(RXPRZ + RXMRZ) ;
      	RZ = 0.5*(RXPRZ - RXMRZ) ;

/*     CHOOSE ONE OF FOUR POSSIBLE WAYS TO CALCULATE RY FOR */
/*     NUMERICAL STABILITY                                  */

      	ELTMAX = CURVUE[0][0] ;
      	X = cos(RZ) ;
      	if(fabs(CURVUE[0][1]) <= fabs(ELTMAX)) goto label_7 ;
      	ELTMAX = CURVUE[0][1] ;
      	X = -sin(RZ) ;

label_7:
    	if(fabs(CURVUE[1][2]) <= fabs(ELTMAX)) goto label_8 ;
      	ELTMAX = CURVUE[1][2] ;
      	X = -sin(RX) ;

label_8:
    	if(fabs(CURVUE[2][2]) <= fabs(ELTMAX)) goto label_9 ;
      	ELTMAX = CURVUE[2][2] ;
      	X = cos(RX) ;

label_9:
    	RX = RTOD*RX ;
      	RY = RTOD*atan2((double)A31,(double)ELTMAX/X) ;
      	RZ = RTOD*RZ ;
      	goto label_20 ;

/*
C   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
C
C    DEGENERATE CASE: RY = 90 OR -90
C
*/
label_10:
   	RX = 0.0 ;
      	RY = (CURVUE[0][2] >= 0.0 ? 90.0 : -90.0 ) ;
      	RZ = RTOD*atan2(CURVUE[1][0],CURVUE[1][1]) ;


/*   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .*/

label_20:
	tempnum1 = (int)((RX + 540.0)*100000) ;
	tempnum2 = (int)((360.0)*100000)      ;
	ALTRX = (float)(tempnum1 % tempnum2)  ;
	tempnum1 = (int)((540.0 -RY)*100000)  ;
	ALTRY = (float)(tempnum1 % tempnum2)  ;
	tempnum1 = (int)((RZ + 540.0)*100000) ;
	ALTRZ = (float)(tempnum1 % tempnum2)  ;
/*	
C
C     DOES THIS SET GIVE THE SMALLER NUMBERS
C
*/
      	if(fabs(RX) + fabs(RY) + fabs(RZ) <= fabs(ALTRX) + fabs(ALTRZ) + fabs(ALTRZ)) 
		goto label_30 ;
/*      
C
C     USE ALTERNATE SET
C
*/
      	RX = ALTRX ;
      	RY = ALTRY ;
      	RZ = ALTRZ ;

label_30:
   	RXV = RX ;
      	RYV = RY ;
      	RZV = RZ ;

}


void  M_SUPERPOSE( int NSTRUC , int NRES , float ***R , float **CMs , float ***ROTMATS , float **RMSDS )
{
	double LAMBDA , MU , NU , SIGMA , M2TRM ;
	//float  ROWT[3][3] , RXV , RYV , RZV ;
	double **E0 , *E0A , *EA , **E ;
	//Alloc the above
	ALLOC_2D( &E0 , NSTRUC , NSTRUC ) ;
	ALLOC_1D( &E0A , NSTRUC ) ;
	ALLOC_1D( &EA , NSTRUC ) ;
	ALLOC_2D( &E , NSTRUC , NSTRUC ) ;
	
	double ****P , PP[4][4] , M[3][3] , **RHO, RHOMAT[4][4] ;
	//Alloc the pointers above
	P = new double*** [NSTRUC] ;
	for( int i = 0 ; i < NSTRUC ; i++ ) P[i] = new double** [NSTRUC] ;
	for( int i = 0 ; i < NSTRUC ; i++ ) 
		for( int j = 0 ; j < NSTRUC ; j++ )
			P[i][j] = new double* [4] ;
	for( int i = 0 ; i < NSTRUC ; i++ ) 
		for( int j = 0 ; j < NSTRUC ; j++ )
			for( int k = 0 ; k < 4 ; k++ )
				P[i][j][k] = new double [4] ;

	ALLOC_2D( &RHO , NSTRUC , 4 ) ;
	//  .  .  .  .  .  .  .  .
	double ROT[3][3] ;
	
/* PHASE I. INITIALIZATION */
      	if( NSTRUC < 2 ) 
	{
		cerr << "Too few structures to superpose: NSTRUCTS = " << NSTRUC << "\n" ;
		exit(0) ;
	}

	MOVE_TO_CENTER_OF_GRAVITY( NSTRUC , NRES , R , CMs ) ;


      	double RNXNM1 = ((double)1.0/(double)(NRES*(NSTRUC-1))) ;

	/* CALC E0IJ AND P */
	
      	for( int ISTRUC = 0 ; ISTRUC < NSTRUC ; ISTRUC++ )
      	{
      		E0A[ISTRUC] = 0.0 ;
      		for( int JSTRUC = 0 ; JSTRUC < NSTRUC ; JSTRUC++ )
		{
      			E0[ISTRUC][JSTRUC] = 0.0 ;

      	 		for( int I = 0; I < 3 ; I++ )
	 			for( int J = 0; J < 3 ; J++ )
    					 M[I][J] = 0.0 ;

			for( int K = 0 ; K < NRES ; K++ )
			{
				for( int I = 0 ; I < 3 ; I++ )
				{
					E0[ISTRUC][JSTRUC] += (( R[JSTRUC][K][I] - R[ISTRUC][K][I] ) * 
						               ( R[JSTRUC][K][I] - R[ISTRUC][K][I] ) ) ;
					for( int J = 0; J < 3 ; J++ )
						M[I][J] += ( R[JSTRUC][K][J] *  R[ISTRUC][K][I]) ;
				}
			}

			if( ISTRUC != JSTRUC) 
				E0A[ISTRUC] += E0[ISTRUC][JSTRUC] ;

			M2TRM = ( M[0][0] + M[1][1] + M[2][2] ) + ( M[0][0] + M[1][1] + M[2][2] ) ;

			for( int I = 0; I < 3 ; I++ )
			 	for( int J = 0 ; J < 3 ; J++ ) 
					P[ISTRUC][JSTRUC][I][J] = M[I][J] + M[J][I] ;

			P[ISTRUC][JSTRUC][0][0] -=  M2TRM ;
			P[ISTRUC][JSTRUC][1][1] -=  M2TRM ;
			P[ISTRUC][JSTRUC][2][2] -=  M2TRM ;
			P[ISTRUC][JSTRUC][3][0]  = M[2][1] - M[1][2] ;
			P[ISTRUC][JSTRUC][0][3]  = M[2][1] - M[1][2] ;
			P[ISTRUC][JSTRUC][3][1]  = M[0][2] - M[2][0] ;
			P[ISTRUC][JSTRUC][1][3]  = M[0][2] - M[2][0] ;
			P[ISTRUC][JSTRUC][3][2]  = M[1][0] - M[0][1] ;
			P[ISTRUC][JSTRUC][2][3]  = M[1][0] - M[0][1] ;
			P[ISTRUC][JSTRUC][3][3]  = 0.0 ;
		}
	}



/* PHASE II.  SUPERPOSE EACH ONE ON STRUCTURE 1 */
      	RHO[0][0] = 0.0 ;
	RHO[0][1] = 0.0 ;
      	RHO[0][2] = 0.0 ;
     	RHO[0][3] = 1.0 ;
	
      	for( int ISTRUC = 1 ; ISTRUC < NSTRUC ; ISTRUC++ )
	{
      		for( int I = 0 ; I < 4 ; I++ )
			for( int J = 0 ; J < 4 ; J++ )
   				PP[I][J] = P[0][ISTRUC][I][J] ;
		
		/*     DIAGONALIZE PP AND PUT HIGHEST EIGENVECTOR INTO RHO[ISTRUC][*]  */
		if(0) 	cout << "********" << NSTRUC << "*********\n" ;
      		EVV4X4(PP,RHO,ISTRUC) ;
	}
	if(0) cout << endl << endl ;


/* PHASE III.  ITERATIVE LOOP OVER SUPERPOSITIONS */
      	double SUMEA = 0.0 , OSUMEA = 0.0 ;
      	for( int ICYC = 0 ; ICYC < 6 ; ICYC++ ) 
      	{
      		for( int ISTRUC = 0 ; ISTRUC < NSTRUC ; ISTRUC++ ) 
		{
      			for( int I = 0 ; I < 4 ; I++ ) 
				for( int J = 0 ; J < 4 ; J++ )
					PP[I][J] = 0.0 ;
      			EA[ISTRUC] = E0A[ISTRUC] ;
 
      			for( int JSTRUC = 0 ; JSTRUC < NSTRUC ; JSTRUC++ )
			{
				if(JSTRUC == ISTRUC) break ; 
 
				LAMBDA = RHO[JSTRUC][0] ;
				MU = RHO[JSTRUC][1] ;
				NU = RHO[JSTRUC][2] ;
				SIGMA = RHO[JSTRUC][3] ;
				RHOMAT[0][0] = SIGMA ;
				RHOMAT[1][0] = -NU ;
				RHOMAT[2][0] = MU ;
				RHOMAT[3][0] = LAMBDA ;
				RHOMAT[0][1] = NU ;
				RHOMAT[1][1] = SIGMA ;
				RHOMAT[2][1] = -LAMBDA ;
				RHOMAT[3][1] = MU ;
				RHOMAT[0][2] = -MU ;
				RHOMAT[1][2] = LAMBDA ;
				RHOMAT[2][2] = SIGMA ;
				RHOMAT[3][2] = NU ;
				RHOMAT[0][3] = -LAMBDA ;
				RHOMAT[1][3] = -MU ;
				RHOMAT[2][3] = -NU ;
				RHOMAT[3][3] = SIGMA ;
				
				for( int I = 0 ; I < 4 ; I++ )
				  for( int J = 0 ; J < 4 ; J++ )
				    for( int K = 0 ; K < 4 ; K++ )
				      for( int L = 0 ; L < 4 ; L++ )
					 PP[I][J] += (RHOMAT[K][J] * P[JSTRUC][ISTRUC][L][K] * RHOMAT[L][I]) ;
			}
      			EVV4X4(PP,RHO,ISTRUC) ;
			
			/*CALCULATE EA = E0A - 2*RHO*PII*RHO */
	
			for( int K = 0 ; K < 4 ; K++ )
				for( int L = 0 ; L < 4 ; L++ )
   					EA[ISTRUC] = EA[ISTRUC]  - 2.0*RHO[ISTRUC][L]*PP[K][L]*RHO[ISTRUC][K] ;
			
      			EA[ISTRUC] = sqrt( EA[ISTRUC]*RNXNM1 ) ;
		}


		
/* PHASE IV.  CALCULATE RESULTS */
		
		SUMEA = 0.0 ;
      		for( int ISTRUC = 0 ; ISTRUC < NSTRUC ; ISTRUC++ )
			SUMEA = SUMEA + EA[ISTRUC] ;
		
		if(0)
		{
			cout << setw(5) << "ICYC" << endl ;
			cout << "ICYC = " << ICYC  ;
		}
	       	
		
/* PHASE V. CONVERGENCE TEST */
      		if( ICYC > 0 ) 
		{
			if( fabs(OSUMEA - SUMEA) <= THRESH ) break; /* CONVERGED OR MAXIMUM NUMBER OF CYCLES */
		}
		else	OSUMEA = SUMEA ;

		if(0) cout << endl << endl << "^^^^^ " << OSUMEA << "   " << SUMEA << "^^^^^\n\n" ;
	}

	
/* PHASE VI. TRANSFORM COORDINATES */	

//      WRITE (*,203)
//  203 FORMAT(1X/' OPTIMAL ROTATIONS ...'/1X)
 
	double SS ;
	for( int I = 0 ; I < NSTRUC ; I++ )
	{
		/* GET ROTATION MATRIX */
 
		LAMBDA = RHO[I][0] ;
		MU     = RHO[I][1] ;
		NU     = RHO[I][2] ;
		SIGMA  = RHO[I][3] ;
		SS = SIGMA * SIGMA ;
 
		ROT[0][0] = LAMBDA*LAMBDA - MU*MU - NU*NU + SS ;
		ROT[1][0] = 2*(LAMBDA*MU - SIGMA*NU) ;
		ROT[2][0] = 2*(LAMBDA*NU + SIGMA*MU) ;
		ROT[0][1] = 2*(LAMBDA*MU + SIGMA*NU) ;
		ROT[1][1] = -LAMBDA*LAMBDA + MU*MU - NU*NU + SS ;
		ROT[2][1] = 2*(MU*NU - SIGMA*LAMBDA) ;
		ROT[0][2] = 2*(LAMBDA*NU - SIGMA*MU) ;
		ROT[1][2] = 2*(MU*NU + SIGMA*LAMBDA) ;
		ROT[2][2] = -LAMBDA*LAMBDA - MU*MU + NU*NU + SS ;
		
		//--------------
		//store rotmats for use outside the routine
		for( int i = 0 ; i < 3; i++ )
			for( int j = 0 ; j < 3 ; j++ )
				ROTMATS[I][j][i] = ROT[i][j] ;

		//print 'em
		if(0)
		{
			cout << endl ;	
			cout << endl ;	
			cout << setprecision(3) << ios::fixed ;
			for( int i = 0 ; i < 3; i++ )
			{
				for( int j = 0 ; j < 3 ; j++ )
					cout << setw(10) << ROTMATS[I][i][j] ;
				cout << endl ;	
			}
		}
		
		//--------------
 
		/* GET EQUIVALENT ANGLES */
 
/*  added
		for( int K = 0 ; K < 3 ; K++ )
			for( int L = 0 ; L < 3 ; L++ )
				ROWT[K][L] = ROT[K][L] ;
 
		MTRXYZ(ROWT,RXV,RYV,RZV) ;
 added */
 
//WRITE (*,204) I,RXV,RYV,RZV
//204 FORMAT(' STRUCTURE',I4,'  ROTATION(DEGREES) = ',3F10.5)

		/* TRANSFORM COORDINATES */
		float XX, YY, ZZ ;
		for( int J = 0 ; J < NRES ; J++ )
		{
			XX = R[I][J][0] ;
			YY = R[I][J][1] ;
			ZZ = R[I][J][2] ;
			R[I][J][0] = ROT[0][0]*XX + ROT[1][0]*YY + ROT[2][0]*ZZ ;
			R[I][J][1] = ROT[0][1]*XX + ROT[1][1]*YY + ROT[2][1]*ZZ ;
			R[I][J][2] = ROT[0][2]*XX + ROT[1][2]*YY + ROT[2][2]*ZZ ;
		}
	}
/*   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   . */
	/* CALCULATE RMS DEVIATIONS */
	float RECNRS = (float)1.0/float(NRES) ;

	for( int I = 0 ; I < NSTRUC ; I++ )
		for( int J = 0 ; J < NSTRUC ; J++ )
		{
			E[J][I] = 0.0 ;
			if(I == J) break ;
			for( int K = 0 ; K < NRES ; K++ )
			{
				E[J][I] = E[J][I] + (R[I][K][0] - R[J][K][0])*(R[I][K][0] - R[J][K][0])
					          + (R[I][K][1] - R[J][K][1])*(R[I][K][1] - R[J][K][1])
					          + (R[I][K][2] - R[J][K][2])*(R[I][K][2] - R[J][K][2]) ;
			}
			E[J][I] = sqrt(RECNRS*E[J][I]) ;
			E[I][J] = E[J][I] ;
		}

	//COPY RMSDS TO REFLECT BACK IN THE CALLER
	for( int I = 0 ; I < NSTRUC ; I++ ) {
		for( int J = 0 ; J < NSTRUC ; J++ ) {
			RMSDS[I][J] = (float)E[I][J] ;
		}
	}
	/*
	cout << "\nRMSD Matrixs\n" ;
	for( int I = 0 ; I < NSTRUC ; I++ )
	{
		for( int J = 0 ; J < NSTRUC ; J++ )
		{
			cout << setw(9) << E[I][J] ;
		}
		cout << endl ;
	}
	*/
	//cout << setprecision(3) ;
	//cout << " (Net RMSD = " ;
	float RMSD = 0 ;
	for( int I = 0 ; I < NSTRUC ; I++ )
		for( int J = I+1 ; J < NSTRUC ; J++ )
			RMSD += (E[I][J] * E[I][J])  ;
	//cout << RMSD ;
	RMSD = sqrt( (RMSD*2)/(NSTRUC *(NSTRUC-1))) ;
	//cout << setw(9) << RMSD << ")";
	//cout << " (Avg RMSD = " ;
	RMSD = 0 ;
	for( int I = 0 ; I < NSTRUC ; I++ )
		for( int J = I+1 ; J < NSTRUC ; J++ )
			RMSD += E[I][J]  ;
	//cout << RMSD ;
	RMSD = (RMSD*2)/(NSTRUC *(NSTRUC-1)) ;
	//cout << setw(9) << RMSD << ") -- " << NRES ;

			
//WRITE (*,205) (I,I=1,NSTRUC)
//205 FORMAT(1X/' MATRIX OF R.M.S. DEVIATIONS ...'/1X/6X,15I8/1X)
//DO 720 I = 1,NSTRUC
//WRITE (*,206) FILENM(I)(1:4),(E(I,J),J=1,NSTRUC)
//206 FORMAT(1X,A4,1X,10F7.3/(6X,10F7.3)) 

	DE_ALLOC_2D( E0 , NSTRUC ) ;
	DE_ALLOC_2D( E , NSTRUC ) ;
	DE_ALLOC_2D( RHO , NSTRUC ) ;
	DE_ALLOC_1D( E0A ) ;
	DE_ALLOC_1D( EA ) ;

	//dealloc 4D  P
	for( int i = 0 ; i < NSTRUC ; i++ )
	{
		for( int j = 0 ; j < NSTRUC ; j++ )
		{
			for( int k = 0 ; k < 4; k++ )
				delete[] P[i][j][k] ;
			delete[] P[i][j] ;
		}
		delete[] P[i] ;
	}
	delete[] P ;
	if(0) cout << endl << endl << "DONE!!!!!\n" ;
}
