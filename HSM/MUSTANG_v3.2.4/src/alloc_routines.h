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
#ifndef MY_ALLOC_H
#define MY_ALLOC_H
// alloc 2D with fixed size on every 2nd dimension.
template <class Type>
void ALLOC_2D( Type ***arr , int dim1 , int dim2 )
{
	(*arr) = new Type* [dim1] ;
	for( int i = 0 ; i < dim1 ; i++ )
		(*arr)[i] = new Type [ dim2 ];
}


// alloc 2D with variable size of 2nd dimension.
template <class Type>
void ALLOC_2D( Type ***arr , int dim1 , int *dim2 )
{
	(*arr) = new Type* [dim1] ;
	for( int i = 0 ; i < dim1 ; i++ )
		(*arr)[i] = new Type [ dim2[i] ];
}
// alloc 2D with variable size of 2nd dimension plus and extra padding.
template <class Type>
void ALLOC_2D( Type ***arr , int dim1 , int *dim2 , int padding )
{
	(*arr) = new Type* [dim1] ;
	for( int i = 0 ; i < dim1 ; i++ )
		(*arr)[i] = new Type [ dim2[i] + padding ];
}

template <class Type>
void ALLOC_1D( Type **arr , int dim ) 
{
	(*arr) = new Type [dim] ;
}

#endif
