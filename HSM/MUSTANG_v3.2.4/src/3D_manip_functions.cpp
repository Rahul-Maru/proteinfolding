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
/* This file contains some useful functions for 3D Geometry manipulation */ 
#include <iostream>
using std::cin; 	using std::cout; 	using std::cerr; 
using std::flush; 	using std::endl; 	using std::ios;
using std::fixed;


#include <iomanip>
using std::setprecision ; using std::setw ; 

#include <fstream>
#include <cstring>
using std::ifstream ;
using std::ofstream ;

#include <cmath>
#include <cstdlib>
#include <string>
#include <cstdio>
#include "3D_manip_functions.h"


/* This function computes direction cosines of a vector given its direction ratios. */ 
void compute_d_cosines( float dratios[3], float dcosines[3] ) 
{
   float sqsum = 1 / ( sqrt( pow(dratios[0],2) + pow(dratios[1],2) + pow(dratios[2],2) ) ) ; 
   for( int i = 0 ; i < 3 ; i++ ) dcosines[i] = dratios[i] * sqsum ; 
}

/* This functions computes the unit vector along the normal to a plane formed by two lines A and B
*/
void compute_normal( float dratiosA[3], float dratiosB[3], float dcosnormal[3] ) 
{
   float dratnormal[3] ; 
   void compute_cross_product( float[], float[], float[] ) ; 
   compute_cross_product(dratiosA,dratiosB,dratnormal) ; 
   compute_d_cosines( dratnormal, dcosnormal ) ; 
}

/* This function computes the cross product of two vectors A and B */
void compute_cross_product( float dratiosA[3], float dratiosB[3], float dratnormal[3] ) 
{ 
   dratnormal[0] = ( dratiosA[1] * dratiosB[2] ) - ( dratiosA[2] * dratiosB[1] ) ; 
   dratnormal[1] = ( dratiosA[2] * dratiosB[0] ) - ( dratiosA[0] * dratiosB[2] ) ; 
   dratnormal[2] = ( dratiosA[0] * dratiosB[1] ) - ( dratiosA[1] * dratiosB[0] ) ; 
}
  
/* This function computes the dot product of two vectors A and B */
void compute_dot_product( float dratiosA[3], float dratiosB[3], float& dotproduct ) 
{
  dotproduct = 0.0 ; 
  for( int i = 0 ; i < 3 ; i++ ) 
  dotproduct += dratiosA[i] * dratiosB[i] ; 
}

/* This function computes the box product ( scalar triple product ) of three vectors A,B,C */
void compute_box_product( float dratiosA[3], float dratiosB[3], float dratiosC[3], float& boxproduct )
{
  float cross[3] ; 
  compute_cross_product(dratiosB,dratiosC,cross); 
  compute_dot_product(dratiosA,cross,boxproduct); 
}

/* This function computes the rotation matrix for rotation through an angle theta about a line 
   whose direction cosines are given by dcosines 
*/
void compute_rotation_matrix( float dcosines[3], float theta, float rotation_matrix[3][3] ) 
{
   float cost = cos(theta) ; 
   float sint = sin(theta) ; 
   float a1 = dcosines[0] * sint ; 
   float a2 = dcosines[1] * sint ; 
   float a3 = dcosines[2] * sint ; 
   float b1 = dcosines[1] * dcosines[2] * ( 1 - cost ) ; 
   float b2 = dcosines[2] * dcosines[0] * ( 1 - cost ) ; 
   float b3 = dcosines[0] * dcosines[1] * ( 1 - cost ) ; 
   for( int i = 0 ; i < 3 ; i++ ) 
    rotation_matrix[i][i] = cost + ( dcosines[i] * dcosines[i] * ( 1 - cost ) ) ; 
   rotation_matrix[0][1] = b3 - a3 ; 
   rotation_matrix[1][0] = b3 + a3 ; 
   rotation_matrix[2][0] = b2 - a2 ; 
   rotation_matrix[0][2] = b2 + a2 ; 
   rotation_matrix[1][2] = b1 - a1 ; 
   rotation_matrix[2][1] = b1 + a1 ; 
}

/* This function rotates initial_vector to final_vector using the rotation matrix 
*/
void rotate_vector( float rotation_matrix[3][3], float initial_vector[3], float final_vector[3] ) 
{
  for( int i = 0 ; i < 3 ; i++ ) 
  { final_vector[i] = 0.0 ; 
     for( int j = 0 ; j < 3 ; j++ ) 
       final_vector[i] += ( initial_vector[j] * rotation_matrix[i][j] ) ; 
  }
}

/* This function rotates a point about the origin using the rotation matrix */
void rotate_point( float rotation_matrix[3][3], float xq, float yq, float zq, float& newx, float& newy, float& newz ) 
{
 float initial_vector[3], final_vector[3] ; 
 initial_vector[0] = xq ; 
 initial_vector[1] = yq ; 
 initial_vector[2] = zq ; 
 rotate_vector( rotation_matrix, initial_vector, final_vector) ;
 newx = final_vector[0] ; 
 newy = final_vector[1] ; 
 newz = final_vector[2] ; 
} 

/* This function performs tetrahedral fixing, i.e given coordinates x1 of the central atom and 
   coordinates x2 and x3 of two attached atoms, and the tetrahedral angle theta, computes the
   coordinates x4 and x5 of the other two atoms attached to the central atom 
*/
void tetrahedral_fix( float x1[3], float x2[3], float x3[3], float theta, float bond_length1, float bond_length2, float x4[3], float x5[3] )
{
  theta = theta * (M_PI/180) ; 
  float angle1 = theta / 2.0 ;
  float v1[3], v2[3], v3[3], v4[3], wrat[3], wcos[3], u1[3], u2[3], wfinal[3] ;
  int i = 0 ; 
  for( i = 0 ; i < 3 ; i++ ) { v1[i] = x2[i] - x1[i] ; v2[i] = x3[i] - x1[i] ; } 
  compute_d_cosines( v1, u1 ) ; 
  compute_d_cosines( v2, u2 ) ; 
  for( i = 0 ; i < 3 ; i++ ) 
  wrat[i] = - 0.5 * ( u1[i] + u2[i] ) ; 
  
  compute_normal( v1, v2, v3 ) ; 
  compute_normal( v3, wrat, v4 ) ; 
  compute_d_cosines( wrat, wcos ) ; 
  float angle2 = -(angle1) ; 
  float rotation_matrix[3][3] ; 
  compute_rotation_matrix( v4, angle1, rotation_matrix ) ;          
  rotate_vector( rotation_matrix, wcos, wfinal ) ; 
  for( i = 0 ; i < 3 ; i++ )
  x4[i] = x1[i] + ( wfinal[i] * bond_length1 ) ; 
  compute_rotation_matrix( v4, angle2, rotation_matrix ) ; 
  rotate_vector( rotation_matrix, wcos, wfinal ) ; 
  for( i = 0 ; i < 3 ; i++ ) 
  x5[i] = x1[i] + ( wfinal[i] * bond_length2 ) ; 
}

/* This function computes the dihedral angle between the planes ABC and BCD defined by the points
   A, B, C and D 
*/
void compute_dihedral_angle( float A[3], float B[3], float C[3], float D[3], float& dihedral_angle)
{
  int i ; 
  float sign = 0.0 ; 
  float AB[3], BC[3], CD[3], normalABC[3], normalBCD[3] ; 
  for( i = 0 ; i < 3 ; i++ ) 
   { 
     AB[i] = B[i] - A[i] ; 
     BC[i] = C[i] - B[i] ; 
     CD[i] = D[i] - C[i] ; 
   }
  compute_normal( AB, BC, normalABC ) ; 
  compute_normal( BC, CD, normalBCD ) ;
  compute_dot_product( normalABC, normalBCD, dihedral_angle ) ; 
  if( dihedral_angle > 1 ) 
  {
	  //cout << setprecision(13) <<  (double)dihedral_angle << endl ; 
	  dihedral_angle = 1 ;
  	  //exit(0);
  }
  else if( dihedral_angle < -1 ) 
  {
	  //cout << "*" << setprecision(13) <<  (double)dihedral_angle << endl ; 
	  dihedral_angle = -1 ;
	  //exit(0);
  }
  dihedral_angle = acos( dihedral_angle ) ; 
  compute_box_product( BC, normalABC, normalBCD, sign ) ; 
  dihedral_angle = ( sign > 0 ) ? dihedral_angle : -dihedral_angle ; 
  dihedral_angle *= 180/M_PI ; 
}
  
/* This function transforms an array of points X to a new coordinate system EX where the point 
   NAT1 is origin, NAT2 lies on the x axis, and NAT3 lies on the xy plane 
*/
void coordinate_transformer( float X[4][500],float EX[4][500], int N, int NAT1, int NAT2, int NAT3){ 
  int i,j,k ; 
  float T[4][4] ; 
  float Tee[4][4] ; 
  float Tphi[4][4] ; 
  float en[4] ; 
  float sq[4] ; 
  for( k = 1 ; k <= 3 ; k++ ) 
  {
    sq[k] = X[k][NAT1] ; 
    for( j = 1 ; j <= N ; j++ ) 
     {  
        X[k][j] -= sq[k] ; 
     }
  }
  
 for( i = 1 ; i <= 3 ; i++ ) 
     sq[i] = X[i][NAT2] * X[i][NAT2] ; 
  
  sq[1] = sqrt( sq[1] + sq[2] + sq[3] ) ; 
  T[1][1] = X[1][NAT2] / sq[1] ; 
  sq[2] = sqrt( sq[2] + sq[3] ) ; 
  sq[3] = 1.0 / sq[2] ; 
  sq[1] = 1.0 - ( T[1][1] * T[1][1] ) ; 
  sq[1] = sqrt( sq[1] ) ; 
  float cssq = 0.5 * ( 1.0 + T[1][1] ) ; 
  float snsq = 1 - cssq ;  

  for( k = 2 ; k <= 3 ; k++ ) 
  { en[k] = X[k][NAT2] * sq[3] ;
    T[1][k] = en[k] * sq[1] ; 
    T[k][1] = - T[1][k] ; 
    sq[k] = en[k] * en[k] ; 
  }

  sq[1] = ( sq[3] - sq[2] ) * snsq ; 
  T[3][3] = cssq - sq[1] ; 
  T[2][2] = cssq + sq[1] ; 
  T[2][3] = snsq * en[2] * en[3] ; 
  T[2][3] *= -2.0 ; 
  T[3][2] = T[2][3] ; 

  for( i = 1 ; i <= 3 ; i++ ) 
  { 
   en[i] = 0 ; 
   for( k =1 ; k <= 3 ; k++ ) 
   en[i] += T[i][k] * X[k][NAT3] ; 
  }
  sq[2] = ( en[2] * en[2] ) + ( en[3] * en[3] ) ; 
  sq[2] = 1.0 / sqrt( sq[2] ) ; 
  Tphi[2][2] = en[2] * sq[2] ; 
  Tphi[3][3] = Tphi[2][2] ; 
  Tphi[2][3] = en[3] * sq[2] ; 
  Tphi[3][2] = -Tphi[2][3] ; 
  Tphi[1][1] = 1 ; 
  Tphi[1][2] = 0 ; 
  Tphi[1][3] = 0 ; 
  Tphi[2][1] = 0 ; 
  Tphi[3][1] = 0 ; 
  for( i = 1 ; i <= 3 ; i++ ) 
 {
   for( j = 1 ; j <= 3 ; j++ ) 
     {
      Tee[i][j] = 0 ; 
      for( k = 1 ; k <= 3 ; k++ ) 
      Tee[i][j] += Tphi[i][k] * T[k][j] ; 
     }
 }
 
 for( j = 1 ; j <= N ; j++ ) 
 {  
     for( i = 1 ; i <= 3 ; i++ ) 
     {
       EX[i][j] = 0 ; 
       for( k = 1 ; k <= 3 ; k++ ) 
       EX[i][j] += Tee[i][k] * X[k][j] ; 
     }
 }

}

/* This function is no longer used. This is a remanant from my graph theoretic based structure
   prediction method. The function has been rewritten in C++ using Kearsley's method of superposition
   which is lots lots elegant.

   This function performs superposition - it finds the best superposition of the N point set XB 
onto the set XA and outputs the corresponding rotation matrix R, the translation vector V, and the 
root - mean - square error in the fit RMSE 
*/ 
void superpose( int N, float XA[4][500], float XB[4][500], float R[4][4], float V[4], float& RMSE)
{ ofstream f("in00.dat"); 
 ofstream f1("in01.dat"); 
 f << N << "\n" ; 
 int i,j ;
 for( i = 1 ; i <= N ; i++ )
 {
    for( j = 1 ; j <= 3 ; j++ )
     { f << XA[j][i] << " " ; 
       f1 << XB[j][i] << " " ; 
     }
     f << "\n" ; 
     f1 << "\n" ; 
 }
f.close();
f1.close();
cout << "\n Invoking Fortran..... " << "\n" ; 
char command[20] = "./super" ;
int tmp = system(command); 
if(!tmp) {
	cerr << "something wrong!\n" ;
	exit(1) ;
}
cout << "\n Back to C++...... " << "\n" ; 
cout << "\n Writing output to file.... " ;  
ifstream f2("out00.dat"); 
 for( i = 1 ; i <= 3 ; i++ )
 {
    for( j = 1 ; j <= 3 ; j++ )
     { f2 >> R[i][j] ; 
     }
 }   
 for( j = 1 ; j <= 3 ; j++ ) 
 { f2 >> V[j] ; }
 f2 >> RMSE ; 
 cout << "\n Done. " ; 
 f2.close();
}

char* itoa( int input_num) 
{
  int trial_value = input_num ;
  char* output_string = new char [100] ; 
  char output[100] ;   
  int digit, i = 0  ; 
  while( trial_value )
  {
    digit = trial_value % 10 ; 
    trial_value /= 10 ; 
    output[i++] = (digit + 48) ;
  }
  for( int j = 0 ; j < i ; j++ ) 
  *( output_string + j ) = output[i-1-j] ; 
  *( output_string + i ) = '\0' ; 
  return output_string ; 
}
           
float normAminusB( float A[3] , float B[3] )
{
	float dratios[3] ;
	dratios[0] = A[0] - B[0] ;
	dratios[1] = A[1] - B[1] ;
	dratios[2] = A[2] - B[2] ;

	return sqrt( pow(dratios[0],2) + pow(dratios[1],2) + pow(dratios[2],2) )  ;
}
