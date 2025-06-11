# sourceFiles := src/ISING.f90 src/rand.f90 src/CubePartition.f90 

# ISING: $(sourceFiles)
# 	gfortran-11 src/ISING.f90 ./obj/rand.o ./obj/CubePartition.o -fcheck=all -Wno-line-truncation -fopenmp -I ./obj -o ./bin/ISING 

# src/rand.f90: obj/rand.o 
# 	gfortran-11 -fopenmp -c $@ -J ./obj -Wno-line-truncation -o $^ 
# src/CubePartition.f90: obj/CubePartition.o 
# 	gfortran-11 -fopenmp -c $@ -J ./obj -Wno-line-truncation -o  $^


# Compiler settings
FC = mpif90
FCFLAGS =  -I/usr/include -Wno-line-truncation -fopenmp -O3 -march=native
MODDIR = -J./obj
INCDIR = -I./obj 
FFTWFLAGS = -lfftw3_omp -lfftw3
# Directories
SRCDIR = src
OBJDIR = obj
BINDIR = bin

# Source files and their corresponding object files
SOURCES = $(SRCDIR)/rand.f90 $(SRCDIR)/CubePartition.f90 $(SRCDIR)/ISING.f90  $(SRCDIR)/atom.f90 $(SRCDIR)/chainMeshCell.f90 $(SRCDIR)/chainMesh.f90 $(SRCDIR)/energyMin.f90 $(SRCDIR)/vecNd.f90 $(SRCDIR)/LLG.f90 $(SRCDIR)/StereographicProjection.f90 $(SRCDIR)/constants.f90 $(SRCDIR)/reciprocalSpaceProcesses.f90 $(SRCDIR)/PT.f90 $(SRCDIR)/PT-Utils.f90 
# Not in use: $(SRCDIR)/RGFlow.f90

OBJECTS = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/%.o,$(SOURCES))

# Make sure directories exist
$(shell mkdir -p $(OBJDIR) $(BINDIR))

# Main target
all: $(BINDIR)/ISING $(BINDIR)/PT

# Linking the final executable
$(BINDIR)/ISING: $(filter-out $(OBJDIR)/PT.o, $(OBJECTS))
	$(FC) $(FCFLAGS) -o $@ $^ $(FFTWFLAGS)

$(BINDIR)/PT: $(filter-out $(OBJDIR)/ISING.o, $(OBJECTS))
	$(FC) $(FCFLAGS) -o $@ $^ $(FFTWFLAGS)
# Generic rule for object files
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) $(MODDIR) -c $< -o $@ $(FFTWFLAGS)

# Dependencies
$(OBJDIR)/ISING.o: $(OBJDIR)/rand.o $(OBJDIR)/CubePartition.o $(OBJDIR)/atom.o $(OBJDIR)/chainMeshCell.o $(OBJDIR)/chainMesh.o $(OBJDIR)/energyMin.o $(OBJDIR)/LLG.o $(OBJDIR)/reciprocalSpaceProcesses.o 
$(OBJDIR)/PT.o: $(OBJDIR)/rand.o $(OBJDIR)/CubePartition.o $(OBJDIR)/atom.o $(OBJDIR)/chainMeshCell.o $(OBJDIR)/chainMesh.o $(OBJDIR)/energyMin.o $(OBJDIR)/LLG.o $(OBJDIR)/reciprocalSpaceProcesses.o 	
$(OBJDIR)/energyMin.o: $(OBJDIR)/chainMesh.o $(OBJDIR)/constants.o $(OBJDIR)/reciprocalSpaceProcesses.o
$(OBJDIR)/chainMesh.o: $(OBJDIR)/atom.o $(OBJDIR)/chainMeshCell.o $(OBJDIR)/vecNd.o $(OBJDIR)/constants.o 
$(OBJDIR)/chainMeshCell.o: $(OBJDIR)/atom.o 
$(OBJDIR)/rand.o: $(OBJDIR)/chainMesh.o 
$(OBJDIR)/CubePartition.o: $(OBJDIR)/chainMesh.o 
$(OBJDIR)/RGFlow.0: $(OBJDIR)/chainMesh.o  
$(OBJDIR)/atom.o:
$(OBJDIR)/vecNd.o: 
$(OBJDIR)/constants.o:
$(OBJDIR)/PT-Utils.o:
$(OBJDIR)/LLG.o: $(OBJDIR)/chainMesh.o $(OBJDIR)/vecNd.o $(OBJDIR)/StereographicProjection.o $(OBJDIR)/reciprocalSpaceProcesses.o	
$(OBJDIR)/StereographicProjection.o: $(OBJDIR)/vecNd.o
$(OBJDIR)/reciprocalSpaceProcesses.0: $(OBJDIR)/chainMesh.o $(OBJDIR)/vecNd.o $(OBJDIR)/constants.o
# Clean target
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(BINDIR)/ISING

.PHONY: all clean
