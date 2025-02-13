# sourceFiles := src/ISING.f90 src/rand.f90 src/CubePartition.f90 

# ISING: $(sourceFiles)
# 	gfortran-11 src/ISING.f90 ./obj/rand.o ./obj/CubePartition.o -fcheck=all -Wno-line-truncation -fopenmp -I ./obj -o ./bin/ISING 

# src/rand.f90: obj/rand.o 
# 	gfortran-11 -fopenmp -c $@ -J ./obj -Wno-line-truncation -o $^ 
# src/CubePartition.f90: obj/CubePartition.o 
# 	gfortran-11 -fopenmp -c $@ -J ./obj -Wno-line-truncation -o  $^


# Compiler settings
FC = gfortran
FCFLAGS =  -Wno-line-truncation -fopenmp -Ofast 
MODDIR = -J./obj
INCDIR = -I./obj

# Directories
SRCDIR = src
OBJDIR = obj
BINDIR = bin

# Source files and their corresponding object files
SOURCES = $(SRCDIR)/rand.f90 $(SRCDIR)/CubePartition.f90 $(SRCDIR)/ISING.f90 $(SRCDIR)/RGFlow.f90 $(SRCDIR)/atom.f90 $(SRCDIR)/chainMeshCell.f90 $(SRCDIR)/chainMesh.f90 $(SRCDIR)/energyMin.f90 $(SRCDIR)/vecNd.f90  



OBJECTS = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/%.o,$(SOURCES))

# Make sure directories exist
$(shell mkdir -p $(OBJDIR) $(BINDIR))

# Main target
all: $(BINDIR)/ISING

# Linking the final executable
$(BINDIR)/ISING: $(OBJECTS)
	$(FC) $(FCFLAGS) -o $@ $^

# Generic rule for object files
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) $(MODDIR) -c $< -o $@

# Dependencies
$(OBJDIR)/ISING.o: $(OBJDIR)/rand.o $(OBJDIR)/CubePartition.o $(OBJDIR)/RGFlow.o $(OBJDIR)/atom.o $(OBJDIR)/chainMeshCell.o $(OBJDIR)/chainMesh.o $(OBJDIR)/energyMin.o 
$(OBJDIR)/energyMin.o: $(OBJDIR)/chainMesh.o 
$(OBJDIR)/chainMesh.o: $(OBJDIR)/atom.o $(OBJDIR)/chainMeshCell.o $(OBJDIR)/vecNd.o 
$(OBJDIR)/chainMeshCell.o: $(OBJDIR)/atom.o 
$(OBJDIR)/rand.o: $(OBJDIR)/chainMesh.o 
$(OBJDIR)/CubePartition.o: $(OBJDIR)/chainMesh.o 
$(OBJDIR)/RGFlow.0: $(OBJDIR)/chainMesh.o  
$(OBJDIR)/atom.o:
$(OBJDIR)/vecNd.o: 

# Clean target
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(BINDIR)/ISING

.PHONY: all clean
