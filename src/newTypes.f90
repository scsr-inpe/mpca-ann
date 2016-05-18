MODULE newTypes

IMPLICIT NONE

TYPE :: Particle
    REAL (kind = 8) :: fitness
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: solution
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
    INTEGER (kind = 8) :: iterPerturbation
    INTEGER (kind = 8) :: maxNFE
    INTEGER (kind = 8) :: hookeMax
    REAL (kind = 8), DIMENSION(6) :: lowerBound
    REAL (kind = 8), DIMENSION(6) :: upperBound
    REAL (kind = 8) :: lo_small
    REAL (kind = 8) :: up_small
    REAL (kind = 8) :: fmin
    REAL (kind = 8) :: emin
    REAL (kind = 8) :: Jr
    REAL (kind = 8) :: epsH
    REAL (kind = 8) :: rho
    CHARACTER*16 :: functionName
    CHARACTER*16 :: typeOpposition
    LOGICAL :: isOppositionEnabled
    LOGICAL :: verbose
    
    !MLP
    INTEGER (kind = 8) :: nClasses
    INTEGER (kind = 8) :: nClassesValidation
    INTEGER (kind = 8) :: nInputs
    INTEGER (kind = 8) :: nOutputs
    REAL (kind = 8) :: targetError
    INTEGER (kind = 8) :: nEpochs
    integer :: loadWeightsBias
    LOGICAL :: haveValidation
    real (8), allocatable, dimension(:,:) :: x, x_valid
    real (8), allocatable, dimension(:,:) :: y, y_valid
    
    CHARACTER*24 :: path_ent
    CHARACTER*24 :: path_sai
    CHARACTER*24 :: path_sai_valid
    CHARACTER*24 :: path_ent_valid
        
END TYPE OptionsMPCA
    
TYPE :: StatusMPCA
    INTEGER (kind = 8) :: NFE
    INTEGER (kind = 8) :: it
    INTEGER (kind = 8) :: higherNFE
    INTEGER (kind = 8) :: lastUpdate
    INTEGER (kind = 8) :: totalNFE
    LOGICAL :: flag
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: minB
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: maxB
    REAL (kind = 8) :: bestObjectiveFunction
    LOGICAL :: fileUpdated
    logical :: doStop
END TYPE StatusMPCA

TYPE :: annConfig
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:,:) :: wh1
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:,:) :: wh2
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:,:) :: ws
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: bh1
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: bh2
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: bs
    integer :: nClasses
    integer :: nClassesValidation
    integer :: nInputs
    integer :: nOutputs
    integer :: hiddenLayers
    integer :: neuronsLayer(2)
    integer :: activationFunction
    real (kind = 8) :: alpha
    real (kind = 8) :: eta
    real (kind = 8) :: MeanSquaredError
    real (kind = 8) :: MeanSquaredErrorValidation
END TYPE annConfig

END MODULE newTypes