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
#include "globals.h"
#include "macros.h"
char *STRUCT_PATH ;
char **struct_names ;
char **struct_paths ;
char *OUTPUT_IDENTIFIER ;

int param_o_flag = NO ; 
int param_a_flag = NO  ; 
int param_i_flag = NO  ;  
int param_f_flag = NO  ; 
int param_P_flag = NO ; 
int param_H_flag = NO ; 
int param_s_flag = YES ;  
int param_p_flag = NO  ;
int param_r_flag = NO ;  

int param_Fhtml_flag = NO , param_Fps_flag = NO ,  
    param_Fpir_flag = NO, param_Ffasta_flag = NO, 
    param_Fmsf_flag = NO ;
int param_D_flag = NO ;
int NSTRUCTS = -1 ;
char OUTPUT_FILENAME_PREFIX[200] = "results";
struct pdb_info **PROT ;
int **secondary_stuct_identifier ;
struct Complete_pdb_info ***COMPLETE_PDBs ;
int *PROT_SIZES ;
int *PDB_SIZES ;
float ***distance_matrices;
float ***Edge_Weights ;
struct residues res_names[] = {
			{ "ALA" , 'A' }, { "CYS" , 'C' }, { "ASP" , 'D' },
			{ "GLU" , 'E' }, { "PHE" , 'F' }, { "GLY" , 'G' },
			{ "HIS" , 'H' }, { "ILE" , 'I' }, { "LYS" , 'K' },
			{ "LEU" , 'L' }, { "MET" , 'M' }, { "ASN" , 'N' },
			{ "PRO" , 'P' }, { "GLN" , 'Q' }, { "ARG" , 'R' },
			{ "SER" , 'S' }, { "THR" , 'T' }, { "VAL" , 'V' },
			{ "TRP" , 'W' }, { "TYR" , 'Y' }, { "ASX" , 'B' },
			{ "GLX" , 'Z' }, { "UNK" , 'X' }
		};
struct ramaplot efimov_regions[] = {
					{    0.0 , -180.0 ,  50.0 ,  -90.0 , 0 }, //synonymous to region-A 
					{    0.0 , -110.0 , 180.0 ,  100.0 , 1 }, //synonymous to region-B
					{    0.0 , -110.0 , -90.0 , -180.0 , 1 }, //synonymous to region-B
					{ -110.0 , -180.0 , 180.0 ,  100.0 , 2 }, //synonymous to region-C
					{ -110.0 , -180.0 , -90.0 , -180.0 , 2 }, //synonymous to region-C
					{    0.0 , -180.0 , 100.0 ,   50.0 , 3 }, //synonymous to region-D
					{  140.0 ,   20.0 ,  80.0 ,  -40.0 , 4 }  //synonymous to region-E
					//others 5 				  //synonymous to region-F
				   };
struct SecondaryStruct_identifier **SS_identifier ;
struct glob_lib *Global_Library;
struct  loc_lib  **Local_Library;
struct  maximal_fragments  **Max_Frags_Library ;
struct merge_lib **Merged_Library;
float ***Extended_Edge_Weights ;
int **I_ALGN ;
char **ALGN ;
int **ALGN_QUALITY ;
int ALGN_LEN ;
float CA_CA_DIAMETER ;
