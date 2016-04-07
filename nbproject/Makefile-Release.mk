#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/Globais.o \
	${OBJECTDIR}/src/ModuloRNA.o \
	${OBJECTDIR}/src/annGeneralization.o \
	${OBJECTDIR}/src/main_generalization.o \
	${OBJECTDIR}/src/mpca.o \
	${OBJECTDIR}/src/mpcaFunctions.o \
	${OBJECTDIR}/src/newTypes.o \
	${OBJECTDIR}/src/normalR8.o \
	${OBJECTDIR}/src/uniformR8.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpcaann

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpcaann: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.f} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpcaann ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/src/Globais.o: src/Globais.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/Globais.o src/Globais.f90

${OBJECTDIR}/src/ModuloRNA.o: src/ModuloRNA.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/ModuloRNA.o src/ModuloRNA.f90

${OBJECTDIR}/src/annGeneralization.o: src/annGeneralization.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/annGeneralization.o src/annGeneralization.f90

${OBJECTDIR}/src/main_generalization.o: src/main_generalization.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/main_generalization.o src/main_generalization.f90

${OBJECTDIR}/src/mpca.o: src/mpca.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/mpca.o src/mpca.f90

${OBJECTDIR}/src/mpcaFunctions.o: src/mpcaFunctions.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/mpcaFunctions.o src/mpcaFunctions.f90

${OBJECTDIR}/src/newTypes.o: src/newTypes.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/newTypes.o src/newTypes.f90

${OBJECTDIR}/src/normalR8.o: src/normalR8.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/normalR8.o src/normalR8.f90

${OBJECTDIR}/src/uniformR8.o: src/uniformR8.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/uniformR8.o src/uniformR8.f90

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mpcaann
	${RM} *.mod

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
