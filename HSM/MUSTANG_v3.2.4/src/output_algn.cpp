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

#include "macros.h"
#include "globals.h"
#include "output_algn.h"
#include "de_alloc_routines.h"

//color codes for postscipt printing
char PSCOLORS[5][20]   ={ 
			  { " 0.8  0.08 0.24 " }, //RED
			  { " 0.25 0.41 0.88 " }, //BLUE
			  { " 0.53 0    0.53 " }, //MAGENTA
	                  { " 0.33 0.42 0.19 " }, //GREEN
			  { " 0    0    0    " }  //BLACK
			};
char PSORIGIN[] = " 72 -72 " ;
char PSLINE_INCR[] = " 15 " ;


void PSPRINT( char str[] , int bold_flag,  int color , int lc , ofstream &outfile , int line_cont_flag ) 
{
	if( bold_flag == ON )
		outfile << "/Courier-Bold findfont 10 scalefont setfont\n" ;
	else
		outfile << "/Courier findfont 10 scalefont setfont\n" ;
	outfile << PSCOLORS[color] << " " << "setrgbcolor" << endl ;
	
	if( line_cont_flag == OFF )
		outfile << PSORIGIN << " " << PSLINE_INCR << " " << lc << " mul sub moveto\n" ;
	outfile << "("<< str << ") show\n" ;
	//if( lc >= 34 ) outfile << "\n\nshowpage\n\n" ;
}
void OUTPUT_ALGN_PIR()
{
	if(!meditate) 
		cout << "Output alignment in pir:   ";
	const int LineWidth = 80 ;
	/* IN PIR */
	char tmpstr[500] = "" ;	
	strcpy( tmpstr , OUTPUT_FILENAME_PREFIX ) ;
			
	strcat(tmpstr , ".pir") ;
	if(!meditate)
	{
		if( strlen(tmpstr) < 18 ) {
		     cout<< "         (...writing " 
			 << "\033[31;1m"   
	       	         << tmpstr  
		         << "\033[0m" << ')'
			 << setw(20-strlen(tmpstr))<< " "  ;
		}
		else {
			cout << "\n" ;
			cout << " (...writing "  ;
		        cout << "\033[31;1m"   
	       	             << tmpstr 
		             << "\033[0m" 
			     << ")"
			     << setw(55-strlen(tmpstr))<< " "  ;
		}
	}	
	char pref_str[500]="" ;
	strcat( pref_str , tmpstr ) ;
	ofstream outfile2( pref_str , ios::out ) ;
	if( !outfile2 )
	{
		cerr << "\nError opening outfile:" << tmpstr << endl ;
		exit(1);
	}

	for( int ISTRUC = 0 ;  ISTRUC < NSTRUCTS ;  ISTRUC++ )
	{
		outfile2 << ">P1;" << struct_names[ISTRUC] << endl ;
		outfile2 << struct_names[ISTRUC] << endl ;
		for( int j = 0 ; j < ALGN_LEN ; j++ )
		{
			if( (j+1)%LineWidth ) outfile2 << ALGN[ISTRUC][j] ;
			else outfile2 << ALGN[ISTRUC][j] << endl ;
		}
		outfile2 << "*\n" ;
		if( ISTRUC+1 < NSTRUCTS ) outfile2 << "\n" ;
	}
	outfile2.close();
	if(!meditate)
		cout << "[ \033[32;4mOK\033[0m ]\n" ;

}

void OUTPUT_ALGN_FASTA() 
{
	if(!meditate) 
		cout << "Output alignment in fasta: ";
	const int LineWidth = 80 ;
	/* IN FASTA */
	char tmpstr[500] = "" ;	
	strcpy( tmpstr , OUTPUT_FILENAME_PREFIX ) ;
			
	strcat(tmpstr , ".afasta") ;
	if(!meditate)
	{
		if( strlen(tmpstr) < 18 ) {
		     cout<< "         (...writing " 
			 << "\033[31;1m"   
	       	         << tmpstr  
		         << "\033[0m" << ')'
			 << setw(20-strlen(tmpstr))<< " "  ;
		}
		else {
			cout << "\n" ;
			cout << " (...writing "  ;
		        cout << "\033[31;1m"   
	       	             << tmpstr 
		             << "\033[0m" 
			     << ")"
			     << setw(55-strlen(tmpstr))<< " "  ;
		}
	}	
	char pref_str[500]="" ;
	strcat( pref_str , tmpstr ) ;
	ofstream outfile2( pref_str , ios::out ) ;
	if( !outfile2 )
	{
		cerr << "\nError opening outfile:" << tmpstr << endl ;
		exit(1);
	}

	for( int ISTRUC = 0 ;  ISTRUC < NSTRUCTS ;  ISTRUC++ )
	{
		int j ;
		outfile2 << ">" << struct_names[ISTRUC] << endl ;
		for( j = 0 ; j < ALGN_LEN ; j++ )
		{
			if( (j+1)%LineWidth ) outfile2 << ALGN[ISTRUC][j] ;
			else outfile2 << ALGN[ISTRUC][j] << endl ;
		}
		if( (j+1)%LineWidth ) outfile2 << "\n" ;
		if( ISTRUC+1 < NSTRUCTS ) outfile2 << "\n" ;
	}
	outfile2.close();

	if(!meditate)
		cout << "[ \033[32;4mOK\033[0m ]\n" ;
}

void OUTPUT_ALGN_MSF( time_t rundate )
{
	if(!meditate) 
		cout << "Output alignment in msf:   ";
	/* IN MSF */
	char tmpstr[500] = "" ;	
	strcpy( tmpstr , OUTPUT_FILENAME_PREFIX ) ;
			
	strcat(tmpstr , ".msf") ;
	if(!meditate)
	{
		if( strlen(tmpstr) < 18 ) {
		     cout<< "         (...writing " 
			 << "\033[31;1m"   
	       	         << tmpstr  
		         << "\033[0m" << ')'
			 << setw(20-strlen(tmpstr))<< " "  ;
		}
		else {
			cout << "\n" ;
			cout << " (...writing "  ;
		        cout << "\033[31;1m"   
	       	             << tmpstr 
		             << "\033[0m" 
			     << ")"
			     << setw(55-strlen(tmpstr))<< " "  ;
		}
	}	
	char pref_str[500]="" ;
	strcat( pref_str , tmpstr ) ;
	ofstream outfile2( pref_str , ios::out ) ;
	if( !outfile2 )
	{
		cerr << "\nError opening outfile:" << tmpstr << endl ;
		exit(1);
	}
	// HEADER PART STARTS
	char *buff =  asctime(localtime(&rundate)) ;
	buff[ strlen(buff)-1 ] = '\0' ;
	outfile2 << "MSF of: " <<  pref_str << " from:    1 to:" << setw(5) << ALGN_LEN << endl ; 
	outfile2 << pref_str << "  MSF:" << setw(5) << ALGN_LEN << " Type: P " << buff << " Check: 9999  .." << endl;
	for( int ISTRUC = 0 ;  ISTRUC < NSTRUCTS ;  ISTRUC++ )
	{
		outfile2.setf(ios::left) ;
		outfile2 << "Name: " << setw(15) << struct_names[ISTRUC] ;
		outfile2.unsetf(ios::left) ;
		outfile2 << "   Len: " << setw(7) << ALGN_LEN ;
		outfile2 << "   Check:    9999" << "   Weight:   1.00\n" ;
	}
	outfile2 << "//\n\n" ;
	// HEADER PART ENDS
	int WIDTH = 50 ;
	for(int COL = 0 ; COL < ALGN_LEN ; COL += WIDTH )
	{
	  for( int i = 0 ; i < NSTRUCTS ; i++ )
	  {
		outfile2.setf( ios::left ) ;
		outfile2 << setw(20) << struct_names[i] ;
		outfile2.unsetf( ios::left ) ;
		for( int j = COL , z = 0 ; z < WIDTH && j < ALGN_LEN ; z++ , j++ )
		{
			
			if(ALGN[i][j] == '-' ) 
				outfile2 << "." ;
			else 
				outfile2 << ALGN[i][j] ;
			if( (z+1)%10 == 0 ) outfile2 << " " ;
		}
	  	outfile2 << endl ;
	  }
	  outfile2 << endl ;
	}

	outfile2.close();
	if(!meditate)
		cout << "[ \033[32;4mOK\033[0m ]\n" ;

}

void OUTPUT_ALGN_HTML( time_t rundate )
{
	if(!meditate) 
		cout << "Output alignment in html:  ";

	int n_iden = 0;
	int n_sim = 0;
	struct residue_color_groups
	{
		char color[10] ;
		char group_list[10] ;
	};
	residue_color_groups  col_groups[4] ={ { "Red"     , "AVFPMILW" },
					       { "Blue"    , "DE"	},
					       { "Magenta" , "RK"       },
					       { "Green"   , "STYHCNGQ" }
					     };  
	char *ALGN_markup = new char [ALGN_LEN] ; for(int i=0;i<ALGN_LEN;i++) ALGN_markup[i]=' ' ;
	char *tmp_col = new char [NSTRUCTS] ; for(int i=0;i<NSTRUCTS;i++) tmp_col[i]=' ' ;
	int *tmp_grp = new int [NSTRUCTS] ; for(int i=0;i<NSTRUCTS;i++) tmp_grp[i]=-9 ;
	
	//calc degrees of residue conservation to output in markup line
	char *type_list = new char [NSTRUCTS] ;
	int *type_list_cnt = new int [NSTRUCTS] ;
	int NTYPES ;
	int gap_flag ;
	for( int j = 0 ; j < ALGN_LEN ; j++ )
	{
		for( int i = 0 ; i < NSTRUCTS ; i++ )
			tmp_col[i] = ALGN[i][j] ;

		if( tmp_col[0] != '-' )
		{
			type_list[0] = tmp_col[0] ;
			type_list_cnt[0] = 1 ;
			NTYPES = 1 ;
			gap_flag = OFF ;

			for( int k = 1 ; k < NSTRUCTS ; k++ )
			{
				if( tmp_col[k] == '-' ){ gap_flag = ON ;  break ; }
				
				int found_flag = OFF ;
				for( int l = 0 ; l < NTYPES ; l++ )
				{
					if( type_list[l] == tmp_col[k] )
					{
						found_flag = ON ;
						type_list_cnt[l]++ ;
						break;
					}
				}
				if( found_flag == OFF )
				{
					type_list[NTYPES]  = tmp_col[k] ;
					type_list_cnt[NTYPES++]  = 1 ;
				}
			}
			if( NTYPES == 1 && gap_flag == OFF ){ ALGN_markup[j] = type_list[0] ;  n_iden++ ; }
			if( NTYPES == 2 && gap_flag == OFF ) 
			{
				if( type_list_cnt[0] > type_list_cnt[1] &&  type_list_cnt[1] == 1 )
				{ 
					ALGN_markup[j] = tolower(type_list[0]) ; 
					n_sim++ ;
				}
				else if( type_list_cnt[1] > type_list_cnt[0] && type_list_cnt[0] == 1)
				{ 
					ALGN_markup[j] = tolower(type_list[1]) ; 
					n_sim++ ; 
				}
			}
				
		}
	}

	/* IN HTML */
	char tmpstr[500] = "" ;	
	/*
	for( int i = strlen(argv)-1 , j = 0 ; i>=0; i-- , j++ )
		if( argv[i] != '/' ) tmpstr[j] = argv[i] ;
		else{  tmpstr[j] = '\0' ; break ; }
	//reverse	
	for( int i = strlen(tmpstr)-1 , j = 0 ; j < i ; i-- , j++  )
		tmpstr[i] ^= tmpstr[j] ^= tmpstr[i] ^= tmpstr[j] ;
	for( int i = strlen(tmpstr)-1 ; i >= 0 ; i-- )
		if( tmpstr[i] == '.') { tmpstr[i] = '\0' ; break ; }
       	*/
	strcpy( tmpstr , OUTPUT_FILENAME_PREFIX ) ;
			
	strcat(tmpstr , ".html") ;
	if(!meditate)
	{
		if( strlen(tmpstr) < 18 ) {
		     cout<< "         (...writing " 
			 << "\033[31;1m"   
	       	         << tmpstr  
		         << "\033[0m" << ')'
			 << setw(20-strlen(tmpstr))<< " "  ;
		}
		else {
			cout << "\n" ;
			cout << " (...writing "  ;
		        cout << "\033[31;1m"   
	       	             << tmpstr 
		             << "\033[0m" 
			     << ")"
			     << setw(55-strlen(tmpstr))<< " "  ;
		}
	}	
	char pref_str[500]="" ;
	strcat( pref_str , tmpstr ) ;
	ofstream outfile2( pref_str , ios::out ) ;
	if( !outfile2 )
	{
		cerr << "\nError opening outfile:" << tmpstr << endl ;
		exit(1);
	}

	outfile2 << "<HTML>\n" ;
	outfile2 << "<BODY bgcolor=white>\n" ;
	outfile2 << "<PRE>\n" ;
	// Header starts
	outfile2 << "################################################################################################\n" ;
	outfile2 << "# Program: <font color=\"Red\">MUSTANG " << VERSION << ": A  Multiple structural alignment algorithm</font>\n";
	outfile2 << "# Authors: <font color=\"Red\">A. S. Konagurthu, J. C. Whisstock, and P. J. Stuckey,  A. M. Lesk</font>\n";
	
	outfile2 << "# Rundate: <font color=\"Brown\">" << asctime(localtime(&rundate)) << "</font>";
	outfile2 << "# Report_file: <font color=\"Brown\">" << tmpstr << "</font>\n";
	outfile2 << "################################################################################################\n" ;
	outfile2 << "#====================================\n" ;
	outfile2 << "# Aligned_structures: <B>" << NSTRUCTS << "</B>\n" ;
	for( int i = 0 ; i < NSTRUCTS ; i++ )
		outfile2 << "# " << setw(3) << i+1 << ": <B>" << struct_names[i] << "</B>\n" ;
	outfile2 << "#" << endl ;
	outfile2 << "# Length: <B>"<< setw(10) << ALGN_LEN << "</B>\n" ;
	outfile2.setf(ios::fixed) ;
	outfile2 << "# Identity: <B>" << setw(8)  << n_iden << "</B>" << "/<B>" << setw(3) << ALGN_LEN << "</B> (<B>" ;
	outfile2 << setw(5) << setprecision(1) << (float)((float)(n_iden*100)/(float)(ALGN_LEN)) << "%</B>)  (Calculated as the percentage of conserved columns in the alignment.)"<< endl ;
	outfile2 << "# Similarity: <B>" << setw(6) <<  n_iden+n_sim << "</B>/<B>" << setw(3) << ALGN_LEN << "</B> (<B>" ;
	outfile2 << setw(5) << setprecision(1) << (float)((float)((n_iden+n_sim)*100)/(float)(ALGN_LEN)) << "%</B>)  (Calculated as the percentage of semi-conserved columns in the alignment)"<< endl ;
	//find #gaps
	int n_gaps = 0 ;
	for( int j = 0 ; j < ALGN_LEN ; j++ )
		for( int i = 0 ; i < NSTRUCTS ; i++ )
			if( ALGN[i][j] == '-' )
			{
				n_gaps++ ;
				break;
			}
	outfile2 << "# Gaps: <B>" << setw(12) << n_gaps << "<B>/</B>" << setw(3) << ALGN_LEN << "</B> (<B>" ;
	outfile2 << setw(5) << setprecision(1) << (float)((float)((n_gaps)*100)/(float)(ALGN_LEN)) << "%</B>)  (Calculated as the percentage of columns with atleast one gap.)"<< endl ;
	outfile2.unsetf(ios::fixed) ;
	outfile2 << "<B>\n" ;
	outfile2 << "#===========================================ALIGNMENT START=========================================\n\n\n" ;
	
	//Header Ends	
	char tmp_color[10] = "";
	int WIDTH = 60 ;
	int *indx_cntr = new int [NSTRUCTS] ; for( int i = 0 ; i < NSTRUCTS ; i++ ) indx_cntr[i] = 0 ;	
	for(int COL = 0 ; COL < ALGN_LEN ; COL += WIDTH )
	{
	  for( int i = 0 ; i < NSTRUCTS ; i++ )
	  {
		outfile2 << "<font color=\"Black\">" ;
		outfile2.setf( ios::left ) ;
		outfile2 << setw(20) << struct_names[i] << "</FONT>" ;
		outfile2.unsetf( ios::left ) ;
		//check if the entire row is a run of nulls
		// if yes, then supress the indx_cntr numbering
		int some_flag =ON;
		for( int j = COL , z = 0 ; z < WIDTH && j < ALGN_LEN ; z++ , j++ )
			if(ALGN[i][j] != '-' )
			{
				some_flag = OFF ;
				break ;
			}
		if(some_flag == OFF)	outfile2 << setw(5) << indx_cntr[i]+1 << "  ";
		else 	outfile2 << setw(5) << " " << "  ";
		for( int j = COL , z = 0 ; z < WIDTH && j < ALGN_LEN ; z++ , j++ )
		{
			//find color 
			
			if( ALGN[i][j] == '-' ) strcpy(tmp_color , "Black") ;
			else
			{
			  indx_cntr[i]++ ;
			  for( int k = 0 ; k < 4; k++ )
			  {
				int found_flag = OFF ;
				for( int l = 0 ; l < (signed)strlen( col_groups[k].group_list) ; l++ )
					if( ALGN[i][j] == col_groups[k].group_list[l] )
					{
						//cout << "("<< ALGN[i][j] << " " << k << ") ";
						found_flag = ON ;
						break ;
					}
				if(found_flag == ON )
				{
					strcpy( tmp_color , col_groups[k].color)  ;
					break ;
				}
			  }
			}
			
			
			if( param_D_flag )
			{
				if( ALGN_QUALITY[i][j] || ALGN[i][j] == '-' ) 
				{
					outfile2 << "<font color=\"" << tmp_color << "\">";
					outfile2 << ALGN[i][j] ;
				}
				else 
				{
					outfile2 << "<font color=\"" << tmp_color << "\"";
					outfile2 << "style=\"BACKGROUND-COLOR: #C4C4C4\">" ;
					outfile2 << (char)tolower(ALGN[i][j]) ;
				}
			}
			else
			{
				outfile2 << "<font color=\"" << tmp_color << "\">";
				outfile2 << ALGN[i][j] ;
			}
			outfile2 << "</font>";

		}
		if( some_flag == OFF ) outfile2 << setw(5) << indx_cntr[i] ;
		else outfile2 << setw(5) << " " ;
		outfile2 << endl ;
	  }
	  //markup line 
	  outfile2 << setw(27) << " " ;
	  for( int j = COL , z = 0 ; z < WIDTH && j < ALGN_LEN ; z++ , j++ )
	  {
		  outfile2 << "<font color=\"Black\">" ;
		  outfile2 << ALGN_markup[j] ;
		  outfile2 << "</font>";
	  }

	  outfile2 << endl ;
	  outfile2 << endl ;
	}
	outfile2 << "\n#=========================================ALIGNMENT END=============================================\n" ;
	outfile2 << "</B>#<font color=\"Brown\"><U>LEGEND:</U></font>\n#\n" ;	
	outfile2 << "# Colours indicate the chemical nature of the amino acid;\n" ;
	outfile2 << "# <B><font color=\"Red\">Red</font></B>         = small hydrophobic including aromatic,{<B><font color=\"Red\">A,F,I,L,M,P,V,W</font></B>}\n" ;
        outfile2 << "# <B><font color=\"Blue\">Blue</font></B>        = Acidic,{<B><font color=\"Blue\">D,E</font></B>}\n" ;
        outfile2 << "# <B><font color=\"Magenta\">Magenta</font></B>     = Basic,{<B><font color=\"Magenta\">K,R</font></B>} and\n" ;
	outfile2 << "# <B><font color=\"Green\">Green</font></B>       = Basic amino acids with hydroxyl groups and/or amine groups {<B><font color=\"Green\">C,G,H,N,Q,S,T,Y</font></B>}.\n#\n" ;

	///*
	if(param_D_flag )
	{
		outfile2 << "# Residues in the alignment that are  shown in <B>";
		char buff[100] = "lower case with grey background" ;
		for( int  i = 0 ; i < (signed)strlen(buff) ; i++ )
		{
			strcpy(tmp_color , "Blue") ;
			for( int k = 0 ; k < 4; k++ )
			{
				int found_flag = OFF ;
				for( int l = 0 ; l < (signed)strlen( col_groups[k].group_list) ; l++ )
				{
					char temp = toupper( buff[i] ) ;
					if( temp == col_groups[k].group_list[l] )
					{
						//cout << "("<< ALGN[i][j] << " " << k << ") ";
						found_flag = ON ;
						break ;
					}
				}
				if(found_flag == ON )
				{
					strcpy( tmp_color , col_groups[k].color)  ;
					break ;
				}
			}
			outfile2 << "<font color=\"" << tmp_color << "\"";
			outfile2 << "style=\"BACKGROUND-COLOR: #C4C4C4\">" ;
			outfile2 << (char)buff[i] ;
			outfile2 << "</font>";
		}
		outfile2 << "</B> are those whose aligned\n" ;
		outfile2 << "# CA-CA distances w.r.t. other members in their column of alignment are above the specified threshold.\n#\n" ;
	}
	//*/
	outfile2 << "# The \"markup row\" below each stretch of the multiple alignment is used to mark completely conserved\n" ;
	outfile2 << "# residue (denoted in <B>UPPERCASE</B>) and semi-conserved reside ( denoted in <B>lowercase</B>) in a column of the alignment.\n#" ;
	
	outfile2 << "\n################################################<B>EOF</B>#################################################\n\n\n" ;
	
	outfile2 << "</PRE>\n" ;
	outfile2 << "</BODY>\n" ;
	outfile2 << "</HTML>\n" ;
	outfile2.close();

	DE_ALLOC_1D( ALGN_markup ) ;
	DE_ALLOC_1D( tmp_col ) ;
	DE_ALLOC_1D( tmp_grp ) ;
	DE_ALLOC_1D( type_list ) ;
	DE_ALLOC_1D( type_list_cnt ) ;
	DE_ALLOC_1D( indx_cntr ) ;
	
	if(!meditate)
		cout << "[ \033[32;4mOK\033[0m ]\n" ;

	
}
void OUTPUT_ALGN_POSTSCRIPT( time_t rundate )
{	
	if(!meditate) 
		cout << "Output alignment in postscript" << std::flush ;

	/* IN POSTSCRIPT */
	int n_iden = 0;
	int n_sim = 0;
	struct residue_color_groups
	{
		char color[10] ;
		char group_list[10] ;
	};
	residue_color_groups  col_groups[4] ={ { "Red"     , "AVFPMILW" },
					       { "Blue"    , "DE"	},
					       { "Magenta" , "RK"       },
					       { "Green"   , "STYHCNGQ" }
					     };  
	char *ALGN_markup = new char [ALGN_LEN] ; for(int i=0;i<ALGN_LEN;i++) ALGN_markup[i]=' ' ;
	char *tmp_col = new char [NSTRUCTS] ; for(int i=0;i<NSTRUCTS;i++) tmp_col[i]=' ' ;
	int *tmp_grp = new int [NSTRUCTS] ; for(int i=0;i<NSTRUCTS;i++) tmp_grp[i]=-9 ;
	
	//calc degrees of residue conservation to output in markup line
	char *type_list = new char [NSTRUCTS] ;
	int *type_list_cnt = new int [NSTRUCTS] ;
	//int NTYPES ;
	//int gap_flag ;


	char tmpstr_ps[500] = "" ;	
	/*
	for( int i = strlen(argv)-1 , j = 0 ; i>=0; i-- , j++ )
		if( argv[i] != '/' ) tmpstr_ps[j] = argv[i] ;
		else{  tmpstr_ps[j] = '\0' ; break ; }
	//reverse	
	for( int i = strlen(tmpstr_ps)-1 , j = 0 ; j < i ; i-- , j++  )
		tmpstr_ps[i] ^= tmpstr_ps[j] ^= tmpstr_ps[i] ^= tmpstr_ps[j] ;
	for( int i = strlen(tmpstr_ps)-1 ; i >= 0 ; i-- )
		if( tmpstr_ps[i] == '.') { tmpstr_ps[i] = '\0' ; break ; }
	*/
	strcpy( tmpstr_ps , OUTPUT_FILENAME_PREFIX ) ;
			
	strcat(tmpstr_ps , ".ps" ) ;
	char pref_str_ps[500] = ""  ;
	strcat( pref_str_ps , tmpstr_ps ) ;
	ofstream outfile_ps( pref_str_ps , ios::out ) ;
	if(!meditate)
	{
		cout  << "\033[35m" << setw(24) ;
	       	cout << tmpstr_ps << "\033[0m" << setw(29)<< " "  ;
	}	

	if( !outfile_ps )
	{
		cerr << "\nError opening outfile:" << tmpstr_ps << endl ;
		exit(1);
	}

	outfile_ps << "%%!PS-Adobe-\n" ;
	outfile_ps << "%%%%Creator: MUSTANG " << VERSION << "\n" ;
	outfile_ps << "%%%%CreationDate:" << asctime(localtime(&rundate)) << "\n" ;
	outfile_ps << "%%%%EndComments\n\n\n" ;
	
	outfile_ps << "90 rotate\n" ;


	
	//outfile_ps << "/Courier-Bold findfont 14 scalefont setfont\n" ;
	char strbuffer[1000] = "" ;
	char strbuffer_temp[1000] = "" ;
	int pscolor_indx = 4 ;
	int line_cntr = 0 ;
	int bold_boolean = ON;
	int line_cont = YES ;
	
	
	// Header starts
	pscolor_indx = 4 ;
	bold_boolean = OFF ;
	line_cont = NO ;
	strcpy( strbuffer, "###################################################################" ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;
	
	
	
	pscolor_indx = 4 ;
	bold_boolean = OFF ;
	line_cont = NO ;
	strcpy(strbuffer , "# Program: ") ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
			
	pscolor_indx = 3 ;
	line_cont = YES ;
	sprintf(strbuffer , "MUSTANG %s: A  Multiple structural alignment algorithm",VERSION) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;

	
	pscolor_indx = 4 ;
	line_cont = NO ;
	strcpy(strbuffer , "# Authors: ") ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;


	pscolor_indx = 3 ;
	line_cont = YES ;
	strcpy(strbuffer , "A. S. Konagurthu, A. M. Lesk, J. C. Whisstock, and P. J. Stuckey") ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;
	

	pscolor_indx = 4 ;
	line_cont = NO ;
	strcpy(strbuffer , "# Rundate: ") ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;


	pscolor_indx = 2 ;
	line_cont = YES ;
	strcpy(strbuffer , asctime(localtime(&rundate)) ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;

	pscolor_indx = 4 ;
	line_cont = NO ;
	strcpy(strbuffer , "# Report_file: ") ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;


	pscolor_indx = 2 ;
	line_cont = YES ;
	strcpy(strbuffer , tmpstr_ps ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;

	pscolor_indx = 4 ;
	line_cont = NO ;
	strcpy( strbuffer, "###################################################################" ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;
	
	strcpy( strbuffer, "#====================================" ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;

	
	strcpy(strbuffer , "# Aligned_structures: " ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;

	bold_boolean = ON ;
	line_cont = YES ;
	sprintf(strbuffer , "%d" , NSTRUCTS ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;

	for( int i = 0 ; i < NSTRUCTS ; i++ )
	{
		bold_boolean = OFF ;
		line_cont = NO ;
		sprintf(strbuffer , "# %3d:" , i+1 ) ;
		PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
		
		bold_boolean = ON ;
		line_cont = YES ;
		sprintf(strbuffer , "%s" , struct_names[i] ) ;
		PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;
	}
	bold_boolean = OFF ;
	line_cont = NO ;
	sprintf(strbuffer , "#") ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;

	bold_boolean = OFF ;
	line_cont = NO ;
	sprintf(strbuffer , "# Length: " ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	
	bold_boolean = ON ;
	line_cont = YES ;
	sprintf(strbuffer , "%10d" , ALGN_LEN ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;
	
	outfile_ps.setf(ios::fixed) ;
	
	//Identity
	pscolor_indx = 4 ;
	bold_boolean = OFF ;
	line_cont = NO ;
	sprintf(strbuffer , "# Identity: " ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	
	bold_boolean = ON ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "%8d", n_iden ) ;
	strcpy( strbuffer , strbuffer_temp ) ; //previously was strcat
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;

	bold_boolean = OFF ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "/") ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	
	bold_boolean = ON ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "%3d", ALGN_LEN ) ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;

	bold_boolean = OFF ;
	line_cont = YES ;
	sprintf( strbuffer_temp , " \\(") ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;

	bold_boolean = ON ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "%5.1f%%", (float)((float)(n_iden*100)/(float)(ALGN_LEN)) ) ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	
	bold_boolean = OFF ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "\\)   \\(Percentage of conserved columns.\\)") ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;

	//Similarity
	pscolor_indx = 4 ;
	bold_boolean = OFF ;
	line_cont = NO ;
	sprintf(strbuffer , "# Similarity: " ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	
	bold_boolean = ON ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "%6d", n_iden+n_sim ) ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;

	bold_boolean = OFF ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "/") ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	
	bold_boolean = ON ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "%3d", ALGN_LEN ) ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;

	bold_boolean = OFF ;
	line_cont = YES ;
	sprintf( strbuffer_temp , " \\(") ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;

	bold_boolean = ON ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "%5.1f%%", (float)((float)((n_iden+n_sim)*100)/(float)(ALGN_LEN)) ) ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	
	bold_boolean = OFF ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "\\)  \\(Percentage of conserved and semi-conserved columns.\\)");
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;


	//find #gaps
	int n_gaps = 0 ;
	for( int j = 0 ; j < ALGN_LEN ; j++ )
		for( int i = 0 ; i < NSTRUCTS ; i++ )
			if( ALGN[i][j] == '-' )
			{
				n_gaps++ ;
				break;
			}
	pscolor_indx = 4 ;
	bold_boolean = OFF ;
	line_cont = NO ;
	sprintf(strbuffer , "# Gaps: " ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	
	bold_boolean = ON ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "%12d", n_gaps ) ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;

	bold_boolean = OFF ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "/") ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	
	bold_boolean = ON ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "%3d", ALGN_LEN ) ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;

	bold_boolean = OFF ;
	line_cont = YES ;
	sprintf( strbuffer_temp , " \\(") ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;

	bold_boolean = ON ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "%5.1f%%", (float)((float)((n_gaps)*100)/(float)(ALGN_LEN)) ) ;
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	
	bold_boolean = OFF ;
	line_cont = YES ;
	sprintf( strbuffer_temp , "\\)  \\(Percentage of columns with atleast one gap.\\)");
	strcpy( strbuffer , strbuffer_temp ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;
	outfile_ps.unsetf(ios::fixed) ;
	
	pscolor_indx = 4 ;
	line_cont = NO ;
	strcpy( strbuffer, "===================================================================" ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	line_cntr += 2;
	//Header Ends	
	
	bold_boolean = ON ;
	int tmp_color_code = 4 ;
	int WIDTH = 60 ;
	//reset indx_cntr
	int *indx_cntr = new int [NSTRUCTS] ;
	for( int i = 0 ; i < NSTRUCTS ; i++ ) indx_cntr[i] = 0 ;	
	for(int COL = 0 ; COL < ALGN_LEN ; COL += WIDTH )
	{
	  for( int i = 0 ; i < NSTRUCTS ; i++ )
	  {
		if( line_cntr  >= 34 )  
		{
			line_cntr = 0 ;
			outfile_ps << "\n\nshowpage\n\n90 rotate\n\n" ;
			
		}
		pscolor_indx = 4 ;
		line_cont = NO ;
		sprintf(strbuffer , "%-10s", struct_names[i] ) ;
		PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;

		//check if the entire row is a run of nulls
		// if yes, then supress the indx_cntr numbering
		int some_flag =ON;
		for( int j = COL , z = 0 ; z < WIDTH && j < ALGN_LEN ; z++ , j++ )
			if(ALGN[i][j] != '-' )
			{
				some_flag = OFF ;
				break ;
			}
		if(some_flag == OFF)	
		{
			line_cont = YES ;
			sprintf(strbuffer , "%5d ", indx_cntr[i]+1 ) ;
			PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
		}
		else 	
		{
			line_cont = YES ;
			sprintf(strbuffer , "%5s ", " " ) ;
			PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
		}
		for( int j = COL , z = 0 ; z < WIDTH && j < ALGN_LEN ; z++ , j++ )
		{
			//find color 
			
			if( ALGN[i][j] == '-' ) tmp_color_code = 4 ;
			else
			{
			  indx_cntr[i]++ ;
			  for( int k = 0 ; k < 4; k++ )
			  {
				int found_flag = OFF ;
				for( int l = 0 ; l < (signed)strlen( col_groups[k].group_list) ; l++ )
					if( ALGN[i][j] == col_groups[k].group_list[l] )
					{
						//cout << "("<< ALGN[i][j] << " " << k << ") ";
						found_flag = ON ;
						break ;
					}
				if(found_flag == ON )
				{
					tmp_color_code = k ;
					break ;
				}
			  }
			}
			
			
			pscolor_indx = tmp_color_code ;
			line_cont = YES ;
			sprintf(strbuffer , "%c", ALGN[i][j] ) ;
			PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
		}
		if( some_flag == OFF ) 
		{	
			pscolor_indx = 4 ;
			line_cont = YES ;
			sprintf(strbuffer , "%5d", indx_cntr[i] ) ;
			PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;
		}
		else 
		{
			pscolor_indx = 4 ;
			line_cont = YES ;
			sprintf(strbuffer , "%5s", " " ) ;
			PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr++ , outfile_ps , line_cont ) ;
		}
	  }
	  //markup line 
	  pscolor_indx = 4 ;
	  line_cont = NO ;
	  sprintf(strbuffer , "%16s", " " ) ;
	  PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	  for( int j = COL , z = 0 ; z < WIDTH && j < ALGN_LEN ; z++ , j++ )
	  {
		  pscolor_indx = 4 ;
		  line_cont = YES ;
		  sprintf(strbuffer , "%c", ALGN_markup[j] ) ;
		  PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	  }
	  line_cntr+=2 ;
	}
	pscolor_indx = 4 ;
	line_cont = NO ;
	sprintf(strbuffer , "%s", "#================================EOF=================================" ) ;
	PSPRINT( strbuffer , bold_boolean, pscolor_indx , line_cntr , outfile_ps , line_cont ) ;
	outfile_ps << "\n\n showpage\n" ;
	outfile_ps.close();
	

	DE_ALLOC_1D( ALGN_markup ) ;
	DE_ALLOC_1D( tmp_col ) ;
	DE_ALLOC_1D( tmp_grp ) ;
	DE_ALLOC_1D( type_list ) ;
	DE_ALLOC_1D( type_list_cnt ) ;
	DE_ALLOC_1D( indx_cntr ) ;

	if(!meditate)
		cout << "[ \033[32;4mOK\033[0m ]\n" ;
	
}
