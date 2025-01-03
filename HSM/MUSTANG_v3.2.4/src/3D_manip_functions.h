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

#ifndef THREE_D_MANIP_H
#define THREE_D_MANIP_H

/* This function computes direction cosines of a vector given its direction ratios. */ 
void compute_d_cosines( float dratios[3], float dcosines[3] ); 

/* This functions computes the unit vector along the normal to a plane formed by two lines A and B*/
void compute_normal( float dratiosA[3], float dratiosB[3], float dcosnormal[3] ) ;

/* This function computes the cross product of two vectors A and B */
void compute_cross_product( float dratiosA[3], float dratiosB[3], float dratnormal[3] ) ;
  
/* This function computes the dot product of two vectors A and B */
void compute_dot_product( float dratiosA[3], float dratiosB[3], float& dotproduct ); 

/* This function computes the box product ( scalar triple product ) of three vectors A,B,C */
void compute_box_product( float dratiosA[3], float dratiosB[3], float dratiosC[3], float& boxproduct );

/* This function computes the rotation matrix for rotation through an angle theta about a line 
   whose direction cosines are given by dcosines 
*/
void compute_rotation_matrix( float dcosines[3], float theta, float rotation_matrix[3][3] ) ;

/* This function rotates initial_vector to final_vector using the rotation matrix 
*/
void rotate_vector( float rotation_matrix[3][3], float initial_vector[3], float final_vector[3] ); 

/* This function rotates a point about the origin using the rotation matrix */
void rotate_point( float rotation_matrix[3][3], float xq, float yq, float zq, float& newx, float& newy, float& newz ); 

/* This function performs tetrahedral fixing, i.e given coordinates x1 of the central atom and 
   coordinates x2 and x3 of two attached atoms, and the tetrahedral angle theta, computes the
   coordinates x4 and x5 of the other two atoms attached to the central atom 
*/
void tetrahedral_fix( float x1[3], float x2[3], float x3[3], float theta, float bond_length1, float bond_length2, float x4[3], float x5[3] );

/* This function computes the dihedral angle between the planes ABC and BCD defined by the points
   A, B, C and D 
*/
void compute_dihedral_angle( float A[3], float B[3], float C[3], float D[3], float& dihedral_angle);
  
/* This function transforms an array of points X to a new coordinate system EX where the point 
   NAT1 is origin, NAT2 lies on the x axis, and NAT3 lies on the xy plane 
*/
void coordinate_transformer( float X[4][500],float EX[4][500], int N, int NAT1, int NAT2, int NAT3);

/* This function performs superposition - it finds the best superposition of the N point set XB 
onto the set XA and outputs the corresponding rotation matrix R, the translation vector V, and the 
root - mean - square error in the fit RMSE 
*/ 
void superpose( int N, float XA[4][500], float XB[4][500], float R[4][4], float V[4], float& RMSE);

/* This function finds the norm of difference between 2 vectors- Let A and B be 2 vectors, it
returns the norm of the vector A-B
*/
float normAminusB( float A[3] , float B[3] ) ;

char* itoa( int input_num) ;
#endif
           
