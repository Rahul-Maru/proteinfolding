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
#ifndef GLOBALS_H
#define GLOBALS_H
//*********************************
//MACROS
#define NSTRUCTS_ULIMIT 20         
#define CP_WINDOW_SIZE 6        
#define MIN_CP_WINDOW_SIZE 6    // default Minimum Contact Pattern Window size used in length-proportional
				// edge-weights routine.
#define MIN_WINDOW_WIDTH 15     // default size of maximal sub-structure for Quaternion-superposition
#define RMSD_THRESH 0.5   
#define RMSD_THRESH_JOINT 0.5   
#define RMSD_THRESH_MAXIMAL 1.80 
#define JOINT_RMSD_THRESH_MAXIMAL 6.5   
#define RMSD_THRESH_4 0.7   
#define RMSD_THRESH_5 0.9   
#define NINFO_LOCAL_LIBRARY 10   
#define ALIGNMENT_RMSD 8.0
//MACRO ENDS HERE
extern const char INSTALL_DIR[] ;

//STRUCT TEMPLATES 
struct residues
{
	char three_letter_code[4] ;
	char single_letter_code ;
};

struct pdb_info  // Contains only the backbone info
{
	char  res_num[10]  ; //as it appears in PDB
	int  res_indx ; //count from 0
	char res_name ;
	float N_coords[3] ;
	float CA_coords[3] ;
	float C_coords[3] ;
	float O_coords[3] ;
	float occupancy ;
	float backbone_occupancies[4] ;
	float Phi , Psi, Omega ;
	pdb_info()
	{
		res_num[0] = '\0' ; res_indx = 0 ; res_name = 'Z' ;
		for( int i = 0 ; i < 3 ; i++ ) 
		{
			N_coords[i] = -999 ; CA_coords[i] = -999 ; 
			C_coords[i] = -999 ; O_coords[i] = -999 ;
		}
		for( int i = 0 ; i < 4 ; i++ ) backbone_occupancies[i] = -9999 ;
		occupancy = -1 ;
		Phi = Psi = Omega = -999 ;
	}
};
struct Complete_pdb_info // Contains entire info; used to generate superposed pdbs
{
	char atom_name[6];
	char residue[4] ;
	char residue_num[6];
	float coords[3] ;
	float sup_coords[3] ;
	float occupancy ;
	float B_factor ;
	struct Complete_pdb_info *link;
};

struct ramaplot
{
	float Phi_Ulimit ;
	float Phi_Llimit ;
	float Psi_Ulimit ;
	float Psi_Llimit ;
	int  Region_Code ;
} ;

struct SecondaryStruct_identifier
{
	int SS_code ;
	int start ;
	int end ;
	SecondaryStruct_identifier()
	{
		SS_code = start = end = -1 ;
	}
} ;

struct glob_lib
{
	int **mates;
	int sizes[2] ;
};
struct maximal_fragments
{
	int start[2] ;
	int size ;
	float rmsd ;
	struct maximal_fragments *link ;
};

struct loc_lib
{
	int frag_start[2] ;
	int frag_size ;
	float frag_rmsd ;
	struct loc_lib *link ;
};
struct merge_lib
{
	int *mates[2];
	int *res_indices[2] ;
	int sizes[2] ;
	int seq_index[2] ;
	struct merge_lib *link;
};

struct tree_info
{
	int node_num ;
	char leaf_flag;
	int child_node_nums[2] ;
	int node_size ;
	int *node_list ;
	float branch_lengths[2] ;
} ;

struct node_composition_from_Tree
{
	int num ;
	int *list ;
} ;

//STRUCT TEMPLATES ENDS HERE
extern char *STRUCT_PATH ;
extern char  **struct_names;
extern char **struct_paths ;
extern char *OUTPUT_IDENTIFIER ;
extern char cmdlineparam_o[500] ; 

extern int param_o_flag  ; 
extern int param_a_flag  ; 
extern int param_i_flag  ; 
extern int param_f_flag  ; 
extern int param_P_flag  ; 
extern int param_H_flag  ; 
extern int param_s_flag  ; 
extern int param_p_flag  ;
extern int param_r_flag  ; 

extern int param_Fhtml_flag , param_Fps_flag,  
       param_Fpir_flag, param_Ffasta_flag , 
       param_Fmsf_flag;

extern int param_D_flag ;
extern int NSTRUCTS ;
extern char OUTPUT_FILENAME_PREFIX[200];
extern struct pdb_info **PROT ;
extern int **secondary_stuct_identifier ;
extern struct Complete_pdb_info ***COMPLETE_PDBs ;
extern struct ramaplot efimov_regions[] ;
extern struct SecondaryStruct_identifier **SS_identifier ;
extern int *PROT_SIZES ;
extern int *PDB_SIZES ;
extern struct residues res_names[] ;
extern float ***distance_matrices;
extern float ***Edge_Weights ;
extern struct glob_lib *Global_Library ;
extern struct  loc_lib  **Local_Library ;
extern struct  maximal_fragments  **Max_Frags_Library ;
extern struct merge_lib  **Merged_Library ;
extern float ***Extended_Edge_Weights;
extern char **ALGN ;
extern int **I_ALGN ;       // contains alignment of residue indices
extern int **ALGN_QUALITY ; // a boolean mask containing which will help to decide whether or not residues 
 			    //in a column of multiple alignment are correctly aligned.
extern int ALGN_LEN ;
extern float CA_CA_DIAMETER ; 
//extern struct glob_lib Pair_Algn[] ;
//extern float Pair_Algn_Costs[];
#endif /*GLOBALS_H*/

