#Compilador Tupa
#CC=ftn 

#Compilador 
CC= mpif90 

# Opcoes de compilacao
CFLAG= -c -O3

# Opcoes de otimizacao
CFLAGOPT= -O3

# Para Tupa
# nr= ../nr2f90/mod/

nr= ../nr2f90/mpif90/mod/

VPATH = src/
MODDIR = mod/
BUILDDIR = build/

all: annMPCA annMLP removemod
	
annMPCA: $(BUILDDIR)/newTypes.o $(BUILDDIR)/uniformR8.o $(BUILDDIR)/normalR8.o \
	$(BUILDDIR)/Globais.o $(BUILDDIR)/ModuloRNA.o\
	$(BUILDDIR)/mpcaFunctions.o $(BUILDDIR)/mpca.o
	$(CC) $(CFLAGOPT) -o annMPCA $(BUILDDIR)/uniformR8.o $(BUILDDIR)/normalR8.o $(BUILDDIR)/Globais.o $(BUILDDIR)/ModuloRNA.o $(BUILDDIR)/newTypes.o $(BUILDDIR)/mpcaFunctions.o $(BUILDDIR)/mpca.o ../nr2f90/mpif90/lib/libnr2f90.a

annMLP:	$(BUILDDIR)/Globais.o $(BUILDDIR)/annGeneralization.o $(BUILDDIR)/main_generalization.o
	$(CC) $(CFLAGOPT) -o annMLP $(BUILDDIR)/Globais.o $(BUILDDIR)/annGeneralization.o $(BUILDDIR)/main_generalization.o

$(BUILDDIR)/%.o: $(VPATH)/%.f90
	$(CC) $(CFLAG) $< -o $@

$(BUILDDIR)/mpca.o: mpca.f90
	$(CC) $(CFLAG) src/mpca.f90 -I $(nr) -o $@

$(BUILDDIR)/oppositionFunctions.o: oppositionFunctions.f90
	$(CC) $(CFLAG) src/oppositionFunctions.f90 -I $(nr) -o $@

$(BUILDDIR)/mpcaFunctions.o: mpcaFunctions.f90 
	$(CC) $(CFLAG) src/mpcaFunctions.f90 -I $(nr) -o $@

clean:
	rm -rf *.*~ Makefile~ *.o *.mod annMLP annMPCA output/*.out

removemod:
	rm -f build/*.o *.mod
