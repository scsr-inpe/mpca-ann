#Compilador Tupa
#CC=ftn 

#Compilador 
CC= mpif90 

# Opcoes de compilacao
CFLAG= -c -O3 -g

# Opcoes de otimizacao
CFLAGOPT= -O3 -g

VPATH = src
MODDIR = mod
BUILDDIR = build

# Arquivos objeto do annMPCA
SRCMPCA := $(BUILDDIR)/foul.o \
$(BUILDDIR)/uniformR8.o \
$(BUILDDIR)/newTypes.o \
$(BUILDDIR)/normalR8.o \
$(BUILDDIR)/annTraining.o \
$(BUILDDIR)/mpcaFunctions.o \
$(BUILDDIR)/mpca.o

# Arquivos objeto do annMLP
SRCMLP := $(BUILDDIR)/foul.o \
$(BUILDDIR)/newTypes.o \
$(BUILDDIR)/annGeneralization.o \
$(BUILDDIR)/main_generalization.o

# Arquivos objeto do annActivation
SRCACTIVATION := $(BUILDDIR)/foul.o \
$(BUILDDIR)/newTypes.o \
$(BUILDDIR)/annActivation.o \
$(BUILDDIR)/main_activation.o


all: 	$(BUILDDIR)/foul.o \
	$(BUILDDIR)/newTypes.o \
	$(BUILDDIR)/uniformR8.o \
	$(BUILDDIR)/normalR8.o \
	$(BUILDDIR)/annTraining.o \
	$(BUILDDIR)/mpcaFunctions.o \
	$(BUILDDIR)/mpca.o \
	$(BUILDDIR)/annActivation.o \
	$(BUILDDIR)/annGeneralization.o \
	$(BUILDDIR)/main_generalization.o \
	$(BUILDDIR)/main_activation.o \
	annMPCA \
	annMLP \
	annActivation

annMPCA:
	$(CC) $(CFLAGOPT) -o annMPCA $(SRCMPCA)

annMLP:
	$(CC) $(CFLAGOPT) -o annMLP $(SRCMLP)

annActivation:
	$(CC) $(CFLAGOPT) -o annActivation $(SRCACTIVATION)

$(BUILDDIR)/%.o: $(VPATH)/%.f90
	@mkdir -p $(@D)
	$(CC) $(CFLAG) $< -o $@

clean:	removemod
	rm -rf *.*~ Makefile~ build/*.o *.mod annActivation annMLP annMPCA

removemod:
	rm -f build/*.o *.mod
