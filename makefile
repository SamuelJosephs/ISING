# ========================
# Compiler settings
# ========================
FC = mpif90
#FCFLAGS = -I/usr/include -Wno-line-truncation -fopenmp -O0 -fcheck=all -fbacktrace -Wpedantic -march=native
FCFLAGS = -I/usr/include -Wno-line-truncation -fopenmp -O3 -Wpedantic -march=native

MODDIR = -J./obj
INCDIR = -I./obj
FFTWFLAGS = -lfftw3_omp -lfftw3

# ========================
# Directories
# ========================
SRCDIR = src
OBJDIR = obj
BINDIR = bin

# ========================
# Source files and objects
# ========================
SOURCES = $(SRCDIR)/rand.f90 \
          $(SRCDIR)/CubePartition.f90 \
          $(SRCDIR)/ISING.f90 \
          $(SRCDIR)/atom.f90 \
          $(SRCDIR)/chainMeshCell.f90 \
          $(SRCDIR)/chainMesh.f90 \
          $(SRCDIR)/energyMin.f90 \
          $(SRCDIR)/vecNd.f90 \
          $(SRCDIR)/LLG.f90 \
          $(SRCDIR)/StereographicProjection.f90 \
          $(SRCDIR)/constants.f90 \
          $(SRCDIR)/reciprocalSpaceProcesses.f90 \
          $(SRCDIR)/PT.f90 \
          $(SRCDIR)/PT-Utils.f90 \
          $(SRCDIR)/algo.f90 \
          $(SRCDIR)/tests.f90 \
          $(SRCDIR)/io.f90 \
          $(SRCDIR)/testWinding.f90 \
	  $(SRCDIR)/unit_cells.f90 \
	  $(SRCDIR)/fft.f90

OBJECTS = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/%.o,$(SOURCES))

# ========================
# Ensure directories exist
# ========================
$(shell mkdir -p $(OBJDIR) $(BINDIR))

# ========================
# Main targets
# ========================
all: $(BINDIR)/ISING $(BINDIR)/PT $(BINDIR)/tests $(BINDIR)/testWinding

# Linking executables
$(BINDIR)/ISING: $(filter-out $(OBJDIR)/PT.o $(OBJDIR)/tests.o $(OBJDIR)/testWinding.o, $(OBJECTS))
	$(FC) $(FCFLAGS) -o $@ $^ $(FFTWFLAGS)

$(BINDIR)/PT: $(filter-out $(OBJDIR)/ISING.o $(OBJDIR)/tests.o $(OBJDIR)/testWinding.o, $(OBJECTS))
	$(FC) $(FCFLAGS) -o $@ $^ $(FFTWFLAGS)

$(BINDIR)/tests: $(filter-out $(OBJDIR)/PT.o  $(OBJDIR)/ISING.o $(OBJDIR)/testWinding.o, $(OBJECTS))
	$(FC) $(FCFLAGS) -o $@ $^ $(FFTWFLAGS)

$(BINDIR)/testWinding: $(filter-out $(OBJDIR)/PT.o $(OBJDIR)/ISING.o $(OBJDIR)/tests.o, $(OBJECTS))
	$(FC) $(FCFLAGS) -o $@ $^ $(FFTWFLAGS)

# ========================
# Compile object files
# ========================
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) $(MODDIR) -c $< -o $@ $(FFTWFLAGS)

# ========================
# Explicit dependencies
# ========================
$(OBJDIR)/ISING.o: $(OBJDIR)/rand.o $(OBJDIR)/CubePartition.o $(OBJDIR)/atom.o $(OBJDIR)/chainMeshCell.o \
                   $(OBJDIR)/chainMesh.o $(OBJDIR)/energyMin.o $(OBJDIR)/LLG.o $(OBJDIR)/reciprocalSpaceProcesses.o \
                   $(OBJDIR)/algo.o $(OBJDIR)/io.o $(OBJDIR)/unit_cells.o

$(OBJDIR)/PT.o: $(OBJDIR)/rand.o $(OBJDIR)/CubePartition.o $(OBJDIR)/atom.o $(OBJDIR)/chainMeshCell.o \
               $(OBJDIR)/chainMesh.o $(OBJDIR)/energyMin.o $(OBJDIR)/LLG.o $(OBJDIR)/reciprocalSpaceProcesses.o \
               $(OBJDIR)/algo.o $(OBJDIR)/io.o $(OBJDIR)/unit_cells.o

$(OBJDIR)/testWinding.o: $(OBJDIR)/rand.o $(OBJDIR)/CubePartition.o $(OBJDIR)/atom.o $(OBJDIR)/chainMeshCell.o \
                        $(OBJDIR)/chainMesh.o $(OBJDIR)/energyMin.o $(OBJDIR)/LLG.o $(OBJDIR)/reciprocalSpaceProcesses.o \
                        $(OBJDIR)/algo.o $(OBJDIR)/io.o $(OBJDIR)/unit_cells.o

$(OBJDIR)/tests.o: $(OBJDIR)/rand.o $(OBJDIR)/CubePartition.o $(OBJDIR)/atom.o $(OBJDIR)/chainMeshCell.o \
                  $(OBJDIR)/chainMesh.o $(OBJDIR)/energyMin.o $(OBJDIR)/LLG.o $(OBJDIR)/reciprocalSpaceProcesses.o \
                  $(OBJDIR)/algo.o $(OBJDIR)/unit_cells.o

$(OBJDIR)/energyMin.o: $(OBJDIR)/chainMesh.o $(OBJDIR)/constants.o $(OBJDIR)/reciprocalSpaceProcesses.o
$(OBJDIR)/chainMesh.o: $(OBJDIR)/atom.o $(OBJDIR)/chainMeshCell.o $(OBJDIR)/vecNd.o $(OBJDIR)/constants.o $(OBJDIR)/algo.o $(OBJDIR)/fft.o
$(OBJDIR)/chainMeshCell.o: $(OBJDIR)/atom.o
$(OBJDIR)/rand.o: $(OBJDIR)/chainMesh.o
$(OBJDIR)/CubePartition.o: $(OBJDIR)/chainMesh.o
$(OBJDIR)/atom.o:
$(OBJDIR)/vecNd.o:
$(OBJDIR)/constants.o:
$(OBJDIR)/PT-Utils.o: $(OBJDIR)/chainMesh.o $(OBJDIR)/constants.o $(OBJDIR)/reciprocalSpaceProcesses.o
$(OBJDIR)/LLG.o: $(OBJDIR)/chainMesh.o $(OBJDIR)/vecNd.o $(OBJDIR)/StereographicProjection.o $(OBJDIR)/reciprocalSpaceProcesses.o
$(OBJDIR)/StereographicProjection.o: $(OBJDIR)/vecNd.o
$(OBJDIR)/reciprocalSpaceProcesses.o: $(OBJDIR)/chainMesh.o $(OBJDIR)/vecNd.o $(OBJDIR)/constants.o $(OBJDIR)/fft.o
$(OBJDIR)/algo.o:
$(OBJDIR)/io.o: $(OBJDIR)/chainMesh.o
$(OBJDIR)/unit_cells.o: $(OBJDIR)/atom.o
$(OBJDIR)/fft.o: 
# ========================
# Clean
# ========================
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(BINDIR)/*

.PHONY: all clean

