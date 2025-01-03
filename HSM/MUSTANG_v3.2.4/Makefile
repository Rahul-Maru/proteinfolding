#******************************************************************************
# TARGETS : ./bin/mustang-3.2.4
# AUTHOR  : ARUN S KONAGURTHU
#******************************************************************************
#directories
MUSTANG = .

#local directories
SRC = $(MUSTANG)/src
OBJ = $(MUSTANG)/obj
BIN = $(MUSTANG)/bin

#compiler options
CPP = g++ 
CPPFLAGS = -traditional -Wall -O3 
LDFLAGS =
.cpp.o:
	$(CPP) $(CPPFLAGS) -c -o $@ $<

#macros
OBJECTS = $(OBJ)/globals.o $(OBJ)/CmdLineParser.o $(OBJ)/distmat.o $(OBJ)/sse_RK.o \
	  $(OBJ)/read_structures.o $(OBJ)/primary_lib_gen.o  \
	  $(OBJ)/pairwise_global_structalgn.o $(OBJ)/pdb_ripper.o  \
	  $(OBJ)/merge_global_local_libs.o \
	  $(OBJ)/ew.o $(OBJ)/refine_pairalgn.o $(OBJ)/superpose_weightedRMS.o\
	  $(OBJ)/superpose_2.o $(OBJ)/jacobi.o  $(OBJ)/3D_manip_functions.o \
	  $(OBJ)/extended_lib_gen.o $(OBJ)/progress_align.o $(OBJ)/neighbour_joining.o \
	  $(OBJ)/upgma.o \
	  $(OBJ)/superpose_on_core.o $(OBJ)/multiple_superposition.o \
	  $(OBJ)/output_algn.o $(OBJ)/main.o
ALL = $(BIN)/mustang-3.2.4

#targets
all: $(ALL)


#------------------------------------------------------------------------------
$(BIN)/mustang-3.2.4: $(OBJECTS) 
	$(CPP) $(CPPFLAGS)  $(LDFLAGS) -o $@ $(OBJECTS)

$(OBJ)/globals.o: $(SRC)/globals.h $(SRC)/macros.h $(SRC)/globals.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/globals.cpp -o $(OBJ)/globals.o
		
$(OBJ)/CmdLineParser.o: $(SRC)/CmdLineParser.h $(SRC)/macros.h $(SRC)/globals.h $(SRC)/de_alloc_routines.h $(SRC)/CmdLineParser_2.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/CmdLineParser_2.cpp -o $(OBJ)/CmdLineParser.o
		
$(OBJ)/pdb_ripper.o: $(SRC)/pdb_ripper.h $(SRC)/macros.h $(SRC)/globals.h $(SRC)/alloc_routines.h $(SRC)/init_routines.h  $(SRC)/pdb_ripper_2.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/pdb_ripper_2.cpp -o $(OBJ)/pdb_ripper.o
		
$(OBJ)/read_structures.o: $(SRC)/read_structures.h $(SRC)/macros.h $(SRC)/globals.h $(SRC)/pdb_ripper.h $(SRC)/alloc_routines.h $(SRC)/init_routines.h $(SRC)/de_alloc_routines.h $(SRC)/read_structures.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/read_structures.cpp -o $(OBJ)/read_structures.o

$(OBJ)/3D_manip_functions.o: $(SRC)/3D_manip_functions.h $(SRC)/3D_manip_functions.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/3D_manip_functions.cpp -o $(OBJ)/3D_manip_functions.o
		
$(OBJ)/distmat.o: $(SRC)/macros.h $(SRC)/globals.h $(SRC)/distmat.h $(SRC)/distmat.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/distmat.cpp -o $(OBJ)/distmat.o
		
$(OBJ)/sse_RK.o: $(SRC)/macros.h $(SRC)/globals.h $(SRC)/sse_RK.h $(SRC)/sse_RK.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/sse_RK.cpp -o $(OBJ)/sse_RK.o
		
$(OBJ)/primary_lib_gen.o: $(SRC)/macros.h  $(SRC)/primary_lib_gen.h $(SRC)/globals.h \
			  $(SRC)/pairwise_global_structalgn.h $(SRC)/primary_lib_gen.cpp \
			  $(SRC)/merge_global_local_libs.h
		$(CPP) $(CPPFLAGS) -c $(SRC)/primary_lib_gen.cpp -o $(OBJ)/primary_lib_gen.o
		
$(OBJ)/pairwise_global_structalgn.o: $(SRC)/macros.h $(SRC)/globals.h  $(SRC)/pairwise_global_structalgn.h \
				     $(SRC)/ew.h $(SRC)/refine_pairalgn.h \
				     $(SRC)/pairwise_global_structalgn.cpp 
		$(CPP) $(CPPFLAGS) -c $(SRC)/pairwise_global_structalgn.cpp -o $(OBJ)/pairwise_global_structalgn.o
		
$(OBJ)/refine_pairalgn.o: $(SRC)/macros.h $(SRC)/globals.h  $(SRC)/refine_pairalgn.h \
				     $(SRC)/de_alloc_routines.h $(SRC)/superpose_2.h \
				     $(SRC)/refine_pairalgn_2.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/refine_pairalgn_2.cpp -o $(OBJ)/refine_pairalgn.o

$(OBJ)/superpose_2.o: $(SRC)/jacobi.h $(SRC)/macros.h $(SRC)/superpose_2.h $(SRC)/superpose_2.cpp 
		$(CPP) $(CPPFLAGS) -c $(SRC)/superpose_2.cpp -o $(OBJ)/superpose_2.o
		
$(OBJ)/superpose_weightedRMS.o: $(SRC)/jacobi.h $(SRC)/macros.h $(SRC)/superpose_weightedRMS.h $(SRC)/superpose_weightedRMS.cpp 
		$(CPP) $(CPPFLAGS) -c $(SRC)/superpose_weightedRMS.cpp -o $(OBJ)/superpose_weightedRMS.o
		
$(OBJ)/jacobi.o: $(SRC)/jacobi.h $(SRC)/jacobi.cpp 
		$(CPP) $(CPPFLAGS) -c $(SRC)/jacobi.cpp -o $(OBJ)/jacobi.o
		
$(OBJ)/ew.o: $(SRC)/globals.h $(SRC)/ew.h $(SRC)/superpose_2.h $(SRC)/superpose_weightedRMS.h \
					$(SRC)/ew_2.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/ew_2.cpp -o $(OBJ)/ew.o
		
$(OBJ)/merge_global_local_libs.o: $(SRC)/merge_global_local_libs.h $(SRC)/macros.h $(SRC)/globals.h \
				  $(SRC)/merge_global_local_libs.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/merge_global_local_libs.cpp -o $(OBJ)/merge_global_local_libs.o
		
$(OBJ)/extended_lib_gen.o: $(SRC)/extended_lib_gen.h $(SRC)/globals.h $(SRC)/macros.h  $(SRC)/extended_lib_gen_3.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/extended_lib_gen_3.cpp -o $(OBJ)/extended_lib_gen.o
		
$(OBJ)/progress_align.o: $(SRC)/progress_align.h $(SRC)/globals.h $(SRC)/macros.h  \
			 $(SRC)/neighbour_joining.h $(SRC)/progress_align_3.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/progress_align_3.cpp -o $(OBJ)/progress_align.o
		
$(OBJ)/neighbour_joining.o: $(SRC)/globals.h  $(SRC)/neighbour_joining.h $(SRC)/neighbour_joining.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/neighbour_joining.cpp -o $(OBJ)/neighbour_joining.o
		
$(OBJ)/upgma.o: $(SRC)/globals.h $(SRC)/upgma.h $(SRC)/upgma.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/upgma.cpp -o $(OBJ)/upgma.o

$(OBJ)/superpose_on_core.o: $(SRC)/superpose_on_core.h $(SRC)/superpose_2.h $(SRC)/macros.h $(SRC)/globals.h $(SRC)/alloc_routines.h  $(SRC)/init_routines.h  $(SRC)/de_alloc_routines.h $(SRC)/read_structures.h $(SRC)/superpose_on_core_2.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/superpose_on_core_2.cpp -o $(OBJ)/superpose_on_core.o


$(OBJ)/multiple_superposition.o: $(SRC)/multiple_superposition.h $(SRC)/jacobi.h $(SRC)/alloc_routines.h $(SRC)/multiple_superposition.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/multiple_superposition.cpp -o $(OBJ)/multiple_superposition.o

$(OBJ)/output_algn.o: $(SRC)/output_algn.h $(SRC)/macros.h $(SRC)/globals.h $(SRC)/output_algn.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/output_algn.cpp -o $(OBJ)/output_algn.o

$(OBJ)/main.o: $(SRC)/macros.h  $(SRC)/globals.h $(SRC)/distmat.h $(SRC)/sse_RK.h $(SRC)/primary_lib_gen.h  $(SRC)/main.cpp
		$(CPP) $(CPPFLAGS) -c $(SRC)/main.cpp -o $(OBJ)/main.o

#------------------------------------------------------------------------------
clean:
	rm $(OBJ)/*.o

remove: clean
	rm  $(BIN)/*
#------------------------------------------------------------------------------

#******************************************************************************
