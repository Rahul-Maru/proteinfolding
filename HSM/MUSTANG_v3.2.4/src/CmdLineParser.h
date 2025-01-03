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
   CREATED BY: ARUN S KONAGURTHU
	 DATE: 06/2005
	 
   A simple class for parsing command lines which purely suits my needs.
*/  

#ifndef CMDLPARSE
#define CMDLPARSE

class CmdLine
{
	private:
		int  NOPTS ;
		int  allocSize ;
		char **option ;
		char ***args ;
		int  *nargs ;
	public:
		 CmdLine() ; /* constructor */
		~CmdLine() ; /* destructor  */
		
/*************   void CmdLine::setNOPTS( int num ) *****************************
   This member function assigns the number of options available in the cmdline.
*/   
		void setNOPTS( int num ) ;
		
/*************   void CmdLine::allocMemory( void ) **********************************
   This member function allocs memory( first-dimension only) to the private member 
   pointers.
*/   
		void allocMemory( void ) ;

/*************   void CmdLine::parseCmdLine( int argc , char **argv ) **********
   This member function counts and stores the various options supplied in the
   command line as well as associate with each of the options their respective
   arguments. This will be later processed according to the requirement of the
   program.
*/   
		void parseCmdLine( int argc , char **argv ) ;


/*************   int CmdLine::getNoptions( void ) ******************************
   This member function returns the number of options supplied in the command 
   line
*/ 
		int getNoptions();

/*************   char** CmdLine::getOptions( void ) ******************************
   This member function returns the complete array of options supplied in the command 
   line
*/ 
		char** getOptions();

		
/*************   int CmdLine::isOption( const char *param_opt ) ****************
   This member function checks if the supplied argument is a command line 
   option. returns 1 if it is; returns 0 if it is NOT.   
   
   Options are detected by the presence of '-' as a prefix followed by 
   alphabetical characters. Ex: -x or -xx ...
*/ 
		int isOption( const char *arg ) ;

/*************   int CmdLine::isOptionDuplicate( const char *arg , int N) *************
   This member function checks if the supplied option has already been stored. If it is
   then it returns the index of the Options array; returns -1 otherwise.
*/ 
		int isOptionDuplicate( const char *arg , int N ) ;

/*************   int CmdLine::isOptionAvailable( const char *param_opt ) *******
   This member function checks if the supplied argument is available in the 
   command line. returns 1 if it does; returns 0 if it does NOT.   
*/ 
		int isOptionAvailable( const char *param_opt ) ;


/*************   int CmdLine::getNArgs( const char *param_opt ) ****************
   This member function returns the number of arguments associated with a given 
   option.
*/ 
		int getNArgs( const char *param_opt ) ;

/*************   char* CmdLine::getArgs( const char *param_opt ) **************
   This member function returns the argument associated with the supplied
   option and argmunent index.
*/ 
		char* getArgs( const char *param_opt , int option_indx ) ;


/*************   void CmdLine::putArgs( const char *arg , const int indx ) ****
   This member function does the following:
   (1) checks if args already exists. If it does resize the arglist and store. 
   If not, create enough space for the new one
*/ 
		void putArgs( const char *arg , int indx ) ;

};



void MUSTANG_CMDLINE_PARSER( int argc, char **argv ) ;
#endif

