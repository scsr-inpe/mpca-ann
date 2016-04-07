MODULE newTypes

IMPLICIT NONE

TYPE :: Particle
    REAL*8 :: fitness
    REAL*8, ALLOCATABLE, DIMENSION(:) :: solution
END TYPE Particle

TYPE :: OptionsMPCA
    INTEGER :: nDimensions
    INTEGER :: typeProbability
    INTEGER :: iCycleBlackboard
    INTEGER :: nProcessors
    INTEGER :: iProcessor
    INTEGER :: iExperiment
    INTEGER :: cycleOpposition
    INTEGER :: nParticlesProcessor
    INTEGER*8 :: iterPerturbation
    INTEGER*8 :: maxNFE
    INTEGER*8 :: hookeMax
    REAL*8, ALLOCATABLE, DIMENSION(:) :: lowerBound
    REAL*8, ALLOCATABLE, DIMENSION(:) :: upperBound
    REAL*8 :: lo_small
    REAL*8 :: up_small
    REAL*8 :: fmin
    REAL*8 :: emin
    REAL*8 :: Jr
    REAL*8 :: epsH
    REAL*8 :: rho
    CHARACTER*16 :: functionName
    CHARACTER*16 :: typeOpposition
    LOGICAL :: isOppositionEnabled
    LOGICAL :: saveEvolution
    LOGICAL :: saveSuccess
    LOGICAL :: verbose
    
    !MLP
    INTEGER*8 :: nClasses
    INTEGER*8 :: nClassesValidation
    INTEGER*8 :: nInputs
    INTEGER*8 :: nOutputs
    REAL*8 :: targetError
    INTEGER*8 :: nEpochs
    LOGICAL :: randomWeights
    LOGICAL :: haveValidation
    
    CHARACTER*24 :: path_ent
    CHARACTER*24 :: path_sai
    CHARACTER*24 :: path_sai_valid
    CHARACTER*24 :: path_ent_valid
        
END TYPE OptionsMPCA
    
TYPE :: StatusMPCA
    INTEGER*8 :: NFE
    INTEGER*8 :: it
    INTEGER*8 :: higherNFE
    INTEGER*8 :: lastUpdate
    INTEGER*8 :: totalNFE
    LOGICAL :: flag
    REAL*8, ALLOCATABLE, DIMENSION(:) :: minB
    REAL*8, ALLOCATABLE, DIMENSION(:) :: maxB
    REAL*8 :: bestObjectiveFunction
    LOGICAL :: fileUpdated
END TYPE StatusMPCA

TYPE :: ConfigANN
    REAL*8, ALLOCATABLE, DIMENSION(:,:) :: wh1
    REAL*8, ALLOCATABLE, DIMENSION(:,:) :: wh2
    REAL*8, ALLOCATABLE, DIMENSION(:) :: bh1
    REAL*8, ALLOCATABLE, DIMENSION(:) :: bh2
    REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ws
    REAL*8, ALLOCATABLE, DIMENSION(:) :: bs
    integer :: nClasses
    integer :: nClassesValidation
    integer :: nInputs
    integer :: nOutputs
    integer :: hideLayers
    integer :: hide_neuron(2)
    integer :: activationFunction
END TYPE ConfigANN

END MODULE newTypes