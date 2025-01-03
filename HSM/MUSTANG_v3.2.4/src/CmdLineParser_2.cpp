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

#include  "globals.h"
#include  "macros.h"
#include "CmdLineParser.h"
#include "de_alloc_routines.h"
/* constructor of CmdLine class*/
CmdLine::CmdLine()
{
	NOPTS = 0 ;
	option  = NULL ;
	args     = NULL ;
	nargs    = NULL ;
}
		
/* destructor of CmdLine class*/
CmdLine::~CmdLine()
{
	//dealloc 'options'
	for( int i = 0 ; i < allocSize ; i++ )
		delete[] option[i] ;
	delete[] option ;

	//dealloc 'args'
	for( int i = 0 ; i < allocSize ; i++ )
	{
		for( int j = 0 ; j < nargs[i] ;  j++ )
			delete[] args[i][j] ;
		delete[] args[i] ;
	}
	delete[] args ;

	//delete 'nargs'
	delete[] nargs ;
}

/*************   int CmdLine::setNOPTS( int num ) *****************************
   This member function assigns the number of options available in the cmdline.
*/   
void CmdLine::setNOPTS( int num ) 
{
	NOPTS = num ;
}


/*************   void CmdLine::allocMemory( void ) **********************************
   This member function allocs memory to the private members.
   (1) allocs 2D array for "option" ; initializes to NULLs
   (2) allocs only the first dimension for "args".
   (3) allocs 1D array for "nargs"; initializes to 0
*/   
void CmdLine::allocMemory( void )
{
	//2D
	option = new char* [NOPTS] ;
	for( int i = 0 ; i < NOPTS ; i++ ) option[i] = new char [50] ;

	//first dimension only
	args = new char** [NOPTS] ;

	//1D
	nargs = new int [NOPTS]  ;

	//initializations
	for( int i = 0 ; i < NOPTS ; i++ ) 
	{
		option[i][0] = '\0' ;
		nargs[i] = 0 ;
	}
}	


/*************   void CmdLine::parseCmdLine( int argc , char **argv ) **********
   This member function counts, and stores the various options supplied in the
   command line as well as associate with each of the options their respective
   arguments. This will be later processed according to the requirement of the
   program.
*/   
void CmdLine::parseCmdLine( int argc , char **argv )
{
	//first count the number of options available in the cmdline
	int cntr = 0 ;
	for( int i = 1 ; i < argc ; i++ ) if( isOption(argv[i]) ) cntr++ ;
	
	setNOPTS( cntr ) ;
	allocSize = cntr ;
	allocMemory() ;

	//now scan the command line and remember what you read; make associations
	//between options and their arguments
	int ocntr = 0 ;
	int indx = -9 ;
	cntr = 0 ;
	for( int i = 1 ; i < argc ; i++ )
	{
		if( isOption( argv[i]) )
		{
			indx = isOptionDuplicate( argv[i] ,  cntr++ ) ;
			if( indx == -1 ) 
			{
				indx = ocntr++ ;
				strcpy( option[indx] , argv[i] ) ;
			}
		}
		// since this is not what we call an option it must be the argument
		// part of it. 
		else	
		{
			if( indx == -9 )
			{
				cerr << "mustang: Unable to associate the argument \'" << argv[i] << "\' with any valid option!\n" ;
				cerr << "Try \'mustang --help\' for more information.\n" ;
				exit(0) ;
			}
			else	putArgs( argv[i] , indx ) ; 
		}
	}
}


/*************   int CmdLine::isOption( const char *arg ) *********************
   This member function checks if the supplied argument is a command line 
   option. returns 1 if it is; returns 0 if it is NOT.   
   
   Options are detected by the presence of '-' as a prefix followed by 
   alphabetical characters. Ex: -x or -xx ...
*/ 
int CmdLine::isOption( const char *arg )
{
	if( arg == NULL || strlen(arg) < 2 )
		return 0;

   	if( arg[0] == '-' )
   	{
		// this could also be a negative number
      		return( !isdigit(arg[1]) );
	}
   	else return 0;
}


/*************   int CmdLine::getNoptions( void ) *****************************
   This member function returns the number of options supplied in the command 
   line
*/ 
int CmdLine::getNoptions()
{
	return NOPTS ;
}

/*************   char** CmdLine::getOptions( void ) ******************************
   This member function returns the complete array of options supplied in the command 
   line
*/ 
char** CmdLine::getOptions()
{
	//make of copy of options
	char **coption = new char* [NOPTS] ;
	for( int i = 0 ; i < NOPTS ; i++ ) coption[i] = new char [50];

	for( int i = 0 ; i < NOPTS ; i++ ) 
		strcpy( coption[i] , option[i] ) ;
	return coption ;
}

		
/*************   int CmdLine::isOptionAvailable( const char *param_opt ) ******
   This member function checks if the supplied argument is available in the 
   command line. returns 1 if it does; returns 0 if it does NOT.   
*/ 
int CmdLine::isOptionAvailable( const char *param_opt ) 
{
	for( int i = 0 ; i < NOPTS ; i++ )
		if( !strcmp(param_opt , option[i] ) )	return 1 ;

	return 0;
}
		

/*************   int CmdLine::isOptionDuplicate( const char *param_opt , int N ) ***************
   This member function checks if the supplied option has already been stored. If it is
   then it returns the index of the Options array; returns -1 otherwise.
*/ 
int CmdLine::isOptionDuplicate( const char *arg , int N )
{
	for( int i = 0 ; i < N ; i++ )
		if( !strcmp( arg , option[i] ) ) {
			this->NOPTS-- ;
			return i ;
		}

	return -1 ;
}


/*************   int CmdLine::getNArgs( const char *param_opt ) ***************
   This member function returns the number of arguments associated with a given 
   option.
*/ 
int CmdLine::getNArgs( const char *param_opt )
{
	int indx = 0 ;
	for( int i = 0 ; i < NOPTS ; i++ ) 
		if( !strcmp( param_opt , option[i] ) ) 
		{   
			indx = i ;
			break ;
		}

	return nargs[indx] ;
}


/*************   char* CmdLine::getArgs( const char *param_opt ) **************
   This member function returns the argument associated with the supplied
   option and argmunent index.
*/ 
char* CmdLine::getArgs( const char *param_opt , int option_indx ) 
{
	int indx = 0 ;
	for( int i = 0 ; i < NOPTS ; i++ ) 
		if( !strcmp( param_opt , option[i] ) ) 
		{   
			indx = i ;
			break ;
		}
	
	//copy argument to return
	char *carg = new char [ strlen(args[indx][option_indx])+1 ] ;
	strcpy( carg , args[indx][option_indx] ) ;
	
	return carg ;
}


/*************   void CmdLine::putArgs( const char *argmt , const int indx ) ****
   This member function does the following:
   (1) checks if args already exists. If it does resize the arglist and store. 
       If not, create enough space for the new one.
*/ 
void CmdLine::putArgs( const char *argmt , int indx ) 
{
	if( nargs[indx] == 0 )
	{
		args[indx] = new char* [1] ;
		args[indx][0] = new char [ strlen(argmt)+1 ] ;
		strcpy( args[indx][0] , argmt ) ;
		nargs[indx]++ ;
	}
	else
	{
		//alloc_space to copy previous args
		char **cargs = new char* [ nargs[indx] ] ;
		for( int i = 0 ; i < nargs[indx] ; i++ ) cargs[i] = new char [ strlen(args[indx][i])+1 ] ;

		//copy args[indx]
		for( int i = 0 ; i < nargs[indx] ; i++ ) 
			strcpy( cargs[i] , args[indx][i] );


		//remove memory and resize to accommodate the new one
		for( int i = 0 ; i < nargs[indx] ; i++ ) delete[] args[indx][i] ;
		delete[] args[indx] ;


		args[indx] = new char* [ nargs[indx]+1 ] ;
		for( int i = 0 ; i < nargs[indx] ; i++ )  args[indx][i] = new char [ strlen(cargs[i])+1 ] ;

		
		//copy back the stuff
		for( int i = 0 ; i < nargs[indx] ; i++ ) 
			strcpy( args[indx][i] , cargs[i] ) ;
			
		args[indx][ nargs[indx] ] = new char [ strlen(argmt) + 1] ;
		strcpy ( args[indx][ nargs[indx] ] , argmt ) ;


		//dealloc the copy
		for( int i = 0 ; i < nargs[indx] ; i++ )
			delete[]  cargs[i] ;
		delete[]  cargs ;

		//increment number of arguments
		nargs[indx]++ ;
	}
}

struct program_options
{
	char opt[50] ;
	char args[100] ;
	char description[1000] ;
} PRGOPTS[10] = {
   { "-p", " <path>",                      
	"Path to the directory holding the (PDB) structures to be aligned." },
   { "-i", " <struct-1> <struct-2> ...",
	"Input structures to be aligned. Note: if -p option is used in the\n \
	 command line, supply only the file names of the structures; if not\n\
	 give the absolute/relative path of each of the input structures." },
   { "-f", " <description file>",    
	"This option is used to AVOID entering the path (-p)  and file name\n\
	 (-i) details in the command line. Instead, to keep the command line\n\
	 short, the user can enter the path and file name details in a \n\
	 \"description\" file and supply it in the command line. The format \n\
	 of the \"description file\" is given in  \'README.txt\'.\n\
	 Note: the options { -p , -i} and  {-f} are mutually exclusive." },
   { "-o", " <output identifier>", 
	 "A common identifier for various outputs of the program. \n\
	 Appropriate extentions (e.g. <identifier>.html, <identifier>.pdb,\n\
	 <identifier>.msf) will be added to this identifier depending on \n\
	 the options the user specifies in the command line.  DEFAULT output\n\
	 identifier: \'results\'" },
   { "-F", " [<html>/<pir>/<fasta>/<msf>]", 
	 "Alignment output format. The choices are: \'html\', \'fasta\', \n\
	 \'pir\', \'msf\'.  DEFAULT format: \'html\'" },
   { "-D", " [CA-CA diameter]", 
	 "Produce an HTML file where the the residues are reported in lower\n\
	 case with grey background when the aligned(superposed) CA-CA diameter\n\
	 of residues in a column of alignment is > the threshold."},
   { "-s", " [ON/OFF]",    
	 "Generate a PDB file containing optimal superposition of all the\n\
	 structures based on the alignment.  DEFAULT:  \'ON\'" },
   { "-r", " [ON/OFF]", 
         "Print RMSD table, Rotation matrices, translation vectors of structures. \n\
	 DEFAULT: \'OFF\'"},
   { "--help", "" ,  "display this help and exits." },
   { "--version", "" ,  "output version information and exits." }
};
const int NPRGOPTS = 10 ;
int CMDLINE_ERROR_FLAG = OFF ;

void CMDL_HELP()
{
	cerr << "MUSTANG Options: \n" ;
	for( int i = 0 ; i < NPRGOPTS ; i++ )
	{
		cerr.setf(ios::left);
		cerr << "    " <<  PRGOPTS[i].opt;
	       	cerr << " " << PRGOPTS[i].args  << "\n        ";
		cerr << PRGOPTS[i].description << endl ;
	}
	
}

void HOWL( char *o )
{
	static int cntr = 0 ;
	if( cntr == 0 ) 
	{
		cerr << "MUSTANG "
		     << VERSION 
		     << " Command Line error!\n" ;
	}

	cerr << "::" << setw(2) << cntr+1 <<"::" << " invalid option:   " << o << endl ;
	cntr++ ;
}
void CHECK_FOR_UNKNOWN_OPTS( CmdLine& CL ) 
{
	int NOPTS = CL.getNoptions() ;
	char **options = CL.getOptions() ;
	for( int i = 0 ; i < NOPTS ; i++ )
	{
		int flag = OFF ;
		for( int j = 0 ; j < NPRGOPTS ; j++ )
		{
			if( !strcmp( options[i] , PRGOPTS[j].opt ) ) 
			{
				flag = ON ;
				break ;
			}
		}
		if( flag == OFF ) 
		{
			if(!CMDLINE_ERROR_FLAG) CMDLINE_ERROR_FLAG = ON ;
			HOWL( options[i] ) ;
		}
			
	}
	
	DE_ALLOC_2D( options , NOPTS ) ;

	if(CMDLINE_ERROR_FLAG)
	{
		cerr << "Try \'mustang --help\' for more information.\n" ;
		exit(0);
	}
}

void CHECK_FOR_CONFLICTS( CmdLine& CL ) 
{
	// options {-i , -p} and {-f} are mutually exclusive
	int iflag = OFF ;
	int pflag = OFF ;
	int fflag  = OFF ;
	
	int NOPTS = CL.getNoptions() ;
	char **options = CL.getOptions() ;
	for( int i = 0 ; i < NOPTS ; i++ )
	{
		if( !strcmp( options[i] , "-p") ) pflag = ON ;
		if( !strcmp( options[i] , "-i") ) iflag = ON ;
		if( !strcmp( options[i] , "-f") ) fflag = ON ;
	}
	DE_ALLOC_2D( options , NOPTS ) ;

	if( (pflag||iflag) && fflag)
	{
		cerr << "mustang: Cannot use one or both of the -p or -i options in conjunction with -f option\n" ;
		cerr << "       : Use either one or both of the -p or -i options or just the -f option\n" ;
		exit(0) ;
	}

}

void PARSE_Info_file( char *fname  ) 
{
	char buffer[1000] ;
	ifstream infile(fname , ios::in ) ;
	//Read the first line; Path where pdb files can be found
	infile.getline(buffer , 1000 ) ;
	 //alloc STRUCT_PATH
	 STRUCT_PATH = new char [ strlen(buffer) +1 ] ;
	 //copy
	 strcpy( STRUCT_PATH , buffer ) ;
	//Read  NSTRUCTS
	 infile.getline(buffer , 1000 ) ;
	 sscanf(buffer,"%d",&NSTRUCTS) ;
	 
	
	//Read names of the structures
	struct_names = new char* [NSTRUCTS] ;
	struct_paths = new char* [NSTRUCTS] ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		infile >> buffer ;
		struct_names[i] = new char [ strlen(buffer) +1 ] ;
		strcpy( struct_names[i] , buffer ) ;

		strcpy( buffer, STRUCT_PATH ) ;
		strcat( buffer, struct_names[i] ) ;
		struct_paths[i] = new char [ strlen(buffer) +1 ] ;
		strcpy( struct_paths[i] , buffer ) ;
	}

	// Check if structures exist at this path
	char buff[500] ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		strcpy( buff , STRUCT_PATH );
		strcat( buff , struct_names[ i ] );
		ifstream pdbfile( buff , ios::in ) ;
		if( !pdbfile )
			cerr << "Error Opening File:" << buff << endl ;
		pdbfile.close() ;
	}
	infile.close() ;
}

void Parse_InfoFile( char *fname )
{
	char buffer[1000] = "" ;
	char cpybuffer[1000] = "" ;
	int pathFlag = NO ;
	int nstruc = 0 ;
	int maxstrlen = 0 ;
	int pathlen = 0 ;
	ifstream infile(fname , ios::in ) ;
	if( !infile )
	{
		cerr << "\nError opening infile:" << fname << endl ;
		exit(1);
	}

	while( !infile.eof() )
	{
		infile.getline(buffer , 1000 ) ;
		strcpy( cpybuffer , "" ) ;
		char *ptr = strtok( buffer, " \t") ;
		
		while( ptr != NULL )
		{
			strcat( cpybuffer , ptr ) ;
			ptr = strtok( NULL , " \t" ) ;
		}
		
		
		if( cpybuffer[0] == '>' ) 
		{
			pathFlag = YES ;
			if( (signed)strlen(cpybuffer) > pathlen ) 
				pathlen = strlen(cpybuffer) ;
			else
				cout << "GAD SAVE ME!!!\n" ;
			
		}
		if( cpybuffer[0] == '+' ) 
		{
			nstruc++ ;
			if( (signed)strlen(cpybuffer) > maxstrlen ) 
				maxstrlen = strlen(cpybuffer) ;
		}
	}
	infile.close() ;
	NSTRUCTS = nstruc ;
	if( NSTRUCTS < 2 )
	{
		cerr << "mustang: unable to find enough structures in the info file: " << fname << ".\n         Are you using the right format???\n" ;
		cerr << "Try \'mustang --help\' for more information.\n" ;
		exit(0);
	}

	struct_names = new char* [ nstruc ] ;
	for( int i = 0 ; i < nstruc ; i++ )
		struct_names[i] =  new char [maxstrlen + pathlen ] ;


	ifstream infile1(fname , ios::in ) ;
	nstruc = 0 ;
	while( !infile1.eof() )
	{
		infile1.getline(buffer , 1000 ) ;
		strcpy( cpybuffer , "" ) ;
		char *ptr = strtok( buffer, " \t") ;
		
		while( ptr != NULL )
		{
			strcat( cpybuffer , ptr ) ;
			ptr = strtok( NULL , " \t" ) ;
		}
		
		
		if( cpybuffer[0] == '>' ) 
		{
			STRUCT_PATH = new char [ strlen(cpybuffer) +1 ] ;
			for( int i = 1 ; i <= (signed)strlen(cpybuffer) ; i++ )
				STRUCT_PATH[i-1] = cpybuffer[i] ;
			strcat( STRUCT_PATH, "/" ) ;
		}
		if( cpybuffer[0] == '+' ) 
		{
			for( int i = 1 ; i <= (signed)strlen(cpybuffer) ; i++ )
				struct_names[nstruc][i-1] = cpybuffer[i] ;

			nstruc++ ;
		}
	}
	infile1.close() ;

	
	struct_paths = new char* [NSTRUCTS] ;
	if( pathFlag )
	{
		char buffer [1000] ;
		for( int i = 0 ; i < NSTRUCTS ; i++ )
		{
			strcpy( buffer, STRUCT_PATH ) ;
			strcat( buffer, struct_names[i] ) ;
			struct_paths[i] = new char [ strlen(buffer) +1 ] ;
			strcpy( struct_paths[i] , buffer ) ;
		}
	}
	else 
	{
		char buffer [1000] ;
		for( int i = 0 ; i < NSTRUCTS ; i++ )
		{
			strcpy( buffer, struct_names[i] ) ;
			struct_paths[i] = new char [ strlen(buffer) +1 ] ;
			strcpy( struct_paths[i], buffer ) ;
		}
	}

	// Check if structures exist at this path
	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		ifstream pdbfile( struct_paths[i] , ios::in ) ;
		if( !pdbfile )
			cerr << "Error Opening File:" << struct_paths[i] << endl ;
		pdbfile.close() ;
	}

	//exit(0) ;

}

void CHECK_OPTION_ARGS( CmdLine& CL ) 
{
	int NOPTS = CL.getNoptions() ;
	char **options = CL.getOptions() ;
	for( int i = 0 ; i < NOPTS ; i++ )
	{
		int indx  = 0 ;
		for( int j = 0 ; j < NPRGOPTS ; j++ )
		{
			if( !strcmp( options[i] , PRGOPTS[j].opt ) ) 
			{
				indx = j ;
				break ;
			}
		}

		char *buffer;
		int n ;
		switch( indx )
		{
			case 0: // option -p
				if( CL.getNArgs( "-p") != 1 ) 
				{
					cerr << "mustang: -p option takes exactly one argument.\n" ;
					cerr << "Try \'mustang --help\' for more information.\n" ;
					exit(0) ;
				}
				buffer = CL.getArgs( "-p" , 0 ) ;
				STRUCT_PATH = new char [ strlen(buffer) +1 ] ;
	 			strcpy( STRUCT_PATH , buffer ) ;
				DE_ALLOC_1D( buffer) ;
				param_p_flag =  YES ;
				break;
				
			case 1: // option -i
				if( (NSTRUCTS = CL.getNArgs( "-i")) < 2 ) 
				{
					cerr << "mustang: -i option takes atleast 2 arguments.\n" ;
					cerr << "Try \'mustang --help\' for more information.\n" ;
					exit(0) ;
				}
				struct_names = new char* [NSTRUCTS] ;
				for( int i = 0 ; i < NSTRUCTS ; i++ )
				{
					buffer = CL.getArgs( "-i" , i ) ;
					struct_names[i]= new char [ strlen(buffer) +1 ] ;
	 				strcpy( struct_names[i] , buffer ) ;
					DE_ALLOC_1D( buffer) ;
				}
				param_i_flag =  YES ;
				break;

			case 2: // option -f
				if( CL.getNArgs( "-f") != 1 ) 
				{
					cerr << "mustang: -f option takes exactly one argument.\n" ;
					cerr << "Try \'mustang --help\' for more information.\n" ;
					exit(0) ;
				}
				buffer = CL.getArgs( "-f" , 0 ) ;
				//PARSE_Info_file( buffer ) ;
				Parse_InfoFile( buffer ) ;
				DE_ALLOC_1D( buffer) ;
				param_f_flag = YES ;
				break;
					
			case 3: // option -o
				if( CL.getNArgs( "-o") != 1 ) 
				{
					cerr << "mustang: -o option takes exactly one argument.\n" ;
					cerr << "Try \'mustang --help\' for more information.\n" ;
					exit(0) ;
				}
				buffer = CL.getArgs( "-o" , 0 ) ;
				//OUTPUT_FILENAME_PREFIX = new char [ strlen(buffer) +1 ] ;
	 			strcpy( OUTPUT_FILENAME_PREFIX , buffer ) ;
				DE_ALLOC_1D( buffer) ;
				param_o_flag = YES ;
				break ;
				
			case 4: // option -F
				if( ( n = CL.getNArgs( "-F") ) < 1 ) 
				{
					cerr << "mustang: -F option takes exactly one arguments.\n" ;
					cerr << "Try \'mustang --help\' for more information.\n" ;
					exit(0) ;
				}
				for( int i = 0 ; i < n ; i++ )
				{
					buffer = CL.getArgs( "-F" , i ) ;
					for( int i = 0 ; i < (signed)strlen(buffer) ; i++ ) buffer[i] = tolower(buffer[i]) ;
					if( !strcmp( buffer , "html"  ) ) param_Fhtml_flag = YES ;
					else if( !strcmp( buffer , "ps"  ) ) param_Fps_flag = YES ;
					else if( !strcmp( buffer , "fasta"  ) ) param_Ffasta_flag = YES ;
					else if( !strcmp( buffer , "pir"  ) ) param_Fpir_flag = YES ;
					else if( !strcmp( buffer , "msf"  ) ) param_Fmsf_flag = YES ;
					else
					{
						cerr << "mustang: \'" << buffer << "\' is an invalid argument for the -F option.\n" ;
						cerr << "       : \'html\', \'ps\', \'fasta\', \'pir\', \'msf\' are the only valid arguments.\n" ;
						cerr << "Try \'mustang --help\' for more information.\n" ;
						exit(0) ;
					}
				}
				
				DE_ALLOC_1D( buffer) ;
				break ;
				
			case 5: // option -D
				if( ( n = CL.getNArgs( "-D") ) >=  2 ) 
				{
					cerr << "mustang: -D option takes either no or just one argument.\n" ;
					cerr << "Try \'mustang --help\' for more information.\n" ;
					exit(0) ;
				}
				if( n == 0 ) CA_CA_DIAMETER = ALIGNMENT_RMSD ;
				else
				{
					buffer = CL.getArgs( "-D" , 0 ) ;
					sscanf( buffer , "%f" , &CA_CA_DIAMETER ) ;
					DE_ALLOC_1D( buffer) ;
				}
				param_D_flag = ON ;
				param_Fhtml_flag = ON ;
				break;	
				
			case 6: // option -s
				if( CL.getNArgs( "-s") != 1 ) 
				{
					cerr << "mustang: -s option takes exactly one argument.\n" ;
					cerr << "Try \'mustang --help\' for more information.\n" ;
					exit(0) ;
				}
				buffer = CL.getArgs( "-s" , 0 ) ;
				for( int i = 0 ; i < (signed)strlen(buffer) ; i++ ) buffer[i] = tolower(buffer[i]) ;
				if( strcmp( buffer , "on"  )&& strcmp(buffer , "off") )
				{
					cerr << "mustang: \"ON\" or \"OFF\" are the only valid arguments for the -s option.\n" ;
					cerr << "Try \'mustang --help\' for more information.\n" ;
					exit(0) ;
				}
				
				if ( !strcmp( buffer , "on"  ) ) param_s_flag = YES ;
				if ( !strcmp( buffer , "off" ) ) param_s_flag = NO  ;
				DE_ALLOC_1D( buffer) ;
				break ;

			case 7:	
				if( CL.getNArgs( "-r") != 1 ) 
				{
					cerr << "mustang: -r option takes exactly one argument -- [ON/OFF]\n" ;
					cerr << "Try \'mustang --help\' for more information.\n" ;
					exit(0) ;
				}
				buffer = CL.getArgs( "-r" , 0 ) ;
				for( int i = 0 ; i < (signed)strlen(buffer) ; i++ ) buffer[i] = tolower(buffer[i]) ;
				if( strcmp( buffer , "on"  )&& strcmp(buffer , "off") )
				{
					cerr << "mustang: \"ON\" or \"OFF\" are the only valid arguments for the -r option.\n" ;
					cerr << "Try \'mustang --help\' for more information.\n" ;
					exit(0) ;
				}
				
				if ( !strcmp( buffer , "on"  ) ) param_r_flag = YES ;
				if ( !strcmp( buffer , "off" ) ) param_r_flag = NO  ;
				DE_ALLOC_1D( buffer) ;
				break ;
			case 8: // option --help
				if( CL.getNArgs( "--help") != 0 ) 
				{
					cerr << "mustang: --help option does not take any arguments.\n" ;
					cerr << "Try \'mustang --help\' for more information.\n" ;
					exit(0) ;
				}
				CMDL_HELP() ;
				exit(0);

			case 9: // option --version
				if( CL.getNArgs( "--version") != 0 ) 
				{
					cerr << "mustang: --version option does not take any arguments.\n" ;
					cerr << "Try \'mustang --help\' for more information.\n" ;
					exit(0) ;
				}
				exit(0);

		}
	}
	DE_ALLOC_2D( options , NOPTS ) ;
	

}

void MUSTANG_CMDLINE_PARSER( int argc, char **argv )
{
	//instaniate a CmdLine object
	CmdLine CL ;

	CL.parseCmdLine( argc , argv ) ;
	CHECK_FOR_UNKNOWN_OPTS( CL ) ;
	CHECK_OPTION_ARGS( CL ) ;
	CHECK_FOR_CONFLICTS( CL ) ;

	if( param_i_flag == YES && param_p_flag == YES )
	{
		char buffer [1000] ;
		struct_paths = new char* [NSTRUCTS] ;
		for( int i = 0 ; i < NSTRUCTS ; i++ )
		{
			strcpy( buffer, STRUCT_PATH ) ;
			strcat( buffer, struct_names[i] ) ;
			struct_paths[i] = new char [ strlen(buffer) +1 ] ;
			strcpy( struct_paths[i] , buffer ) ;
		}
	}
	else if( param_i_flag == YES && param_p_flag == NO )
	{
		char buffer [1000] ;
		struct_paths = new char* [NSTRUCTS] ;
		for( int i = 0 ; i < NSTRUCTS ; i++ )
		{
			strcpy( buffer, struct_names[i] ) ;
			struct_paths[i] = new char [ strlen(buffer) +1 ] ;
			strcpy( struct_paths[i], buffer ) ;
		}
	}
	else if( param_i_flag == NO && param_f_flag == NO )
	{
		cerr << "mustang: Impossible to proceed without input structures!\n" ;
		cerr << "Try \'mustang --help\' for more information.\n" ;
		exit(0);
	}

	//deciding the formats in which output should be generated.
	if( param_Fhtml_flag == OFF && param_Fps_flag == OFF && 
	    param_Fpir_flag == OFF && param_Ffasta_flag == OFF &&
	    param_Fmsf_flag == OFF ) param_Fhtml_flag = ON ;
	

	//the following piece of code makes sure that the struct_names does not contain
	//the absolute/relative path prefixes.

	for( int i = 0  ; i < NSTRUCTS ; i++ )
	{
		char buffer[1000] = "" ;
		for( int j = strlen(struct_names[i])-1 , k = 0 ; j >= 0 ; j-- , k++ )
			if( struct_names[i][j] != '/' )	buffer[k] = struct_names[i][j] ;
			else buffer[k] =  '\0' ;
			
		//reverse
		char temp;	
		for( int j = 0 , k = strlen(buffer)-1 ; j < k ; j++ , k-- )
		{
			temp = buffer[j] ;
			buffer[j] = buffer[k] ;
			buffer[k] = temp;
			//buffer[j] ^= buffer[k] ^= buffer[j] ^= buffer[k] ;
		}

		strcpy( struct_names[i] , buffer ) ;
	}
	
	
}
