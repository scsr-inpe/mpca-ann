!***************************************************************!
! MULTIPLE PARTICLE COLLISION ALGORITHM (MPCA)                  !
!***************************************************************!
! Developed by: Eduardo Favero Pacheco da Luz (CAP/INPE)        !
! Modified by: Reynier Hernandez Torres (CAP/INPE)              !
! Based on PCA (Wagner F. Sacco)                                !
! Updated: 23-Apr-2015                                          !
!***************************************************************!

PROGRAM MPCA
    USE nr
    USE nrtype
    USE ran_state, ONLY: ran_seed
    USE ModuloRNA
    USE Globais
    USE mpcaFunctions
    USE newTypes

    IMPLICIT NONE
    INCLUDE 'mpif.h'

    !**********************
    ! VARIABLES DEFINITION 
    !**********************
    INTEGER :: iSeed, contD, contP, iCycleBlackboard, iError
    INTEGER :: iInitialProcessor, iFinalProcessor
    INTEGER :: nSeed, un, istat, it, nDimensions, paramRNA
    INTEGER, ALLOCATABLE, DIMENSION(:) :: vSeed, vSeed2
    REAL*8, ALLOCATABLE, DIMENSION(:) :: minB, maxB
    REAL*8, ALLOCATABLE, DIMENSION(:) :: realSolution
    REAL(SP), ALLOCATABLE, DIMENSION(:) :: vRandom
    REAL(SP) :: rRandom
    CHARACTER(len = 10) :: string, str0, str1
    character(len = 100) :: filename
    character(len = 100) :: format_string
    logical :: tryFixedInitial

    PARAMETER (paramRNA = 17)
    DOUBLE PRECISION, DIMENSION(paramRNA) :: parametros

    TYPE (Particle), ALLOCATABLE, DIMENSION(:) :: oldParticle
    TYPE (Particle), ALLOCATABLE, DIMENSION(:) :: newParticle
    TYPE (Particle), ALLOCATABLE, DIMENSION(:) :: bestParticle
    TYPE (Particle) :: bestParticleProcessor
    TYPE (OptionsMPCA) :: opMpca
    TYPE (StatusMPCA) :: stMpca

    !******************
    ! INITIALIZING MPI 
    !******************
    CALL MPI_Init(iError)
    CALL MPI_Comm_Size(MPI_COMM_WORLD, opMpca % nProcessors, iError)
    CALL MPI_Comm_Rank(MPI_COMM_WORLD, opMpca % iProcessor, iError)

    !*********************************************
    ! SETTING PARAMETERS                          
    ! 1 - Number of the experiment                
    ! 2 - Maximum number of function evaluations  
    ! 3 - Opposition Enabled                      
    ! 4 - Type of opposition                      
    !*********************************************
    CALL getarg(1, string)
    READ (string, '(I10)') opMpca % iExperiment
    CALL getarg(2, string)
    READ (string, '(I10)') opMpca % maxNFE
    CALL getarg(3, string)
    IF (string == 'true') THEN
        opMpca % isOppositionEnabled = .TRUE.
        CALL getarg(4, opMpca % typeOpposition)
    ELSE
        opMpca % isOppositionEnabled = .FALSE.
        opMpca % typeOpposition = 'MPCA'
    END IF

    opMpca % nDimensions = 6
    opMpca % saveEvolution = .false.
    opMpca % saveSuccess = .false.
    opMpca % typeProbability = 1
    opMpca % iCycleBlackboard = 500
    opMpca % nParticlesProcessor = 4
    opMpca % iterPerturbation = 500
    opMpca % lo_small = 0.7
    opMpca % up_small = 1.1
    opMpca % emin = 1.0e-7
    opMpca % Jr = 0.01
    opMpca % epsH = 1.0e-11
    opMpca % rho = 0.80
    opMpca % hookeMax = 100000
    opMpca % verbose = .false.
    
    opMpca % path_ent = 'data/x.txt'
    opMpca % path_sai = 'data/yd.txt'
    opMpca % path_sai_valid = 'data/yd_valid.txt'
    opMpca % path_ent_valid = 'data/x_valid.txt'

    !Particles by Processor
    iInitialProcessor = (opMpca % iProcessor * opMpca % nParticlesProcessor) + 1
    iFinalProcessor = iInitialProcessor + opMpca % nParticlesProcessor - 1

    !Output files
    IF (opMpca % iProcessor < 10) THEN
        WRITE (str0, '(I1)') opMpca % iProcessor
    ELSE
        WRITE (str0, '(I2)') opMpca % iProcessor
    END IF

    IF (opMpca % iExperiment < 10) THEN
        WRITE (str1, '(I1)') opMpca % iExperiment
    ELSE
        WRITE (str1, '(I2)') opMpca % iExperiment
    END IF

    if (opMpca % saveEvolution .eqv. .true.) then
        OPEN(UNIT = 10, &
        FILE = './output/evolution' // trim(str1) // '_' // trim(str0) // '.out', &
        ACCESS = 'APPEND')
    end if

    if (opMpca % iProcessor == 0) then
        OPEN(UNIT = 20, FILE = './output/final.out', ACCESS = 'APPEND')
    end if

    if (opMpca % saveSuccess .eqv. .true.) then
        OPEN(UNIT = 30, FILE = './output/success' // trim(str0) // '.dat')
    end if

    !Allocating space for dynamic variables
    ALLOCATE(vRandom(opMpca % nDimensions))
    ALLOCATE(minB(opMpca % nDimensions))
    ALLOCATE(maxB(opMpca % nDimensions))
    ALLOCATE(oldParticle(opMpca % nParticlesProcessor))
    ALLOCATE(newParticle(opMpca % nParticlesProcessor))
    ALLOCATE(bestParticle(opMpca % nParticlesProcessor))
    ALLOCATE(opMpca % lowerBound(opMpca % nDimensions))
    ALLOCATE(opMpca % upperBound(opMpca % nDimensions))
    ALLOCATE(stMpca % minB(opMpca % nDimensions))
    ALLOCATE(stMpca % maxB(opMpca % nDimensions))
    DO contP = 1, opMpca % nParticlesProcessor
        ALLOCATE(oldParticle(contP) % solution(opMpca % nDimensions))
        ALLOCATE(newParticle(contP) % solution(opMpca % nDimensions))
        ALLOCATE(bestParticle(contP) % solution(opMpca % nDimensions))
    END DO
    ALLOCATE(realSolution(opMpca % nDimensions))

    ! SET BOUNDARIES
    OPEN(UNIT = 40, FILE = './config/boundsMPCA.in', STATUS = "old")
    DO contD = 1, opMpca % nDimensions
        READ(40, *) opMpca % lowerBound(contD), opMpca % upperBound(contD)
    ENDDO
    CLOSE(40)

    !Le o restante das informações de entrada   
    OPEN(UNIT = 30, FILE = './config/parametersMPCA.in', STATUS = "old")
    READ(30, *) opMpca % nClasses
    READ(30, *) opMpca % nClassesValidation
    READ(30, *) opMpca % nInputs
    READ(30, *) opMpca % nOutputs
    READ(30, *) opMpca % targetError
    READ(30, *) opMpca % nEpochs
    READ(30, *) opMpca % randomWeights
    READ(30, *) opMpca % haveValidation
    CLOSE(30)

    ! RANDOM NUMBER CONFIGURATION
    CALL random_seed(size = nSeed)
    ALLOCATE(vSeed2(nSeed))
    OPEN(un, FILE = "/dev/urandom", ACCESS = "stream", FORM = "unformatted", &
    ACTION = "read", STATUS = "old", iostat = istat)
    READ(un) vSeed2
    CLOSE(un)
    iSeed = vSeed2(1) + vSeed2(2) / opMpca % iExperiment

    ALLOCATE(vSeed(25))
    CALL ran_seed(sequence = iSeed)

    stMpca % NFE = 0
    stMpca % lastUpdate = 0
    stMpca % totalNFE = 0
    stMpca % higherNFE = 0
    stMpca % flag = .false.
    stMpca % bestObjectiveFunction = huge(0.D0)

    DO contP = 1, opMpca % nParticlesProcessor
        bestParticle(contP) % fitness = huge(0.D0)
    end do
    bestParticleProcessor % fitness = huge(0.D0)

    !*****************************
    ! CREATING INITIAL POPULATION 
    !*****************************

    DO contP = 1, opMpca % nParticlesProcessor
        if (contP == 1) then
            open(UNIT = 40, FILE = './config/iniConfigMPCA.in', STATUS = "old")
            read(40, *) tryFixedInitial
            if (tryFixedInitial .eqv. .true.) then
                read(40, *) oldParticle(contP) % solution(1)
                read(40, *) oldParticle(contP) % solution(2)
                read(40, *) oldParticle(contP) % solution(3)
                read(40, *) oldParticle(contP) % solution(4)
                read(40, *) oldParticle(contP) % solution(5)
                read(40, *) oldParticle(contP) % solution(6)
            else
                call ran3(vRandom)
                do contD = 1, opMpca % nDimensions
                    oldParticle(contP) % solution(contD) = (DBLE(vRandom(contD)) &
                    *(opMpca % upperBound(contD) - opMpca % lowerBound(contD))) &
                    +opMpca % lowerBound(contD)
                end do
            end if
            close(40)
        else
            call ran3(vRandom)
            do contD = 1, opMpca % nDimensions
                oldParticle(contP) % solution(contD) = (DBLE(vRandom(contD)) &
                *(opMpca % upperBound(contD) - opMpca % lowerBound(contD))) &
                +opMpca % lowerBound(contD)
            end do
        end if

        oldParticle(contP) % fitness = Rede_Neural_BP(oldParticle(contP) % solution, opMpca, stMpca)
        stMpca % NFE = stMpca % NFE + 1

        IF (oldParticle(contP) % fitness < bestParticle(contP) % fitness) THEN
            bestParticle(contP) = oldParticle(contP)
        END IF

    END DO

    DO contP = 1, opMpca % nParticlesProcessor
        if (bestParticle(contP) % fitness < bestParticleProcessor % fitness) then
            bestParticleProcessor = bestParticle(contP)
        end if
    END DO
    
    if(opMpca % verbose .eqv. .true.) then
        write(*, '(I2)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(1))
        write(*, '(I3)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(2))
        write(*, '(I3)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(3))
        write(*, '(I2)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(4))
        write(*, '(S,ES10.2E2)', ADVANCE = 'NO') bestParticleProcessor % solution(5)
        write(*, '(S,ES10.2E2)', ADVANCE = 'NO') bestParticleProcessor % solution(6)
        write(*, '(S,ES10.2E2)', ADVANCE = 'NO') bestParticleProcessor % fitness
        write(*, *) stMpca % totalNFE
    end if

    !***************************************************************************
    ! PRINCIPAL LOOP                                                            
    !.AND. (abs(bestParticleProcessor % fitness - fmin) .GE. emin))             
    !***************************************************************************
    DO WHILE ((stMpca % higherNFE .LE. opMpca % maxNFE / opMpca % nProcessors) &
        .AND. (stMpca % NFE .LE. opMpca % maxNFE / opMpca % nProcessors) &
        .AND. (stMpca % totalNFE .LE. opMpca % maxNFE) &
        )
        !PARTICLES LOOP (OpenMP)
        DO contP = 1, opMpca % nParticlesProcessor
            CALL Perturbation(oldParticle(contP), newParticle(contP), bestParticle(contP), opMpca, stMpca)
            IF (newParticle(contP) % fitness < oldParticle(contP) % fitness) THEN
                IF (newParticle(contP) % fitness < bestParticle(contP) % fitness) THEN
                    bestParticle(contP) = newParticle(contP)
                END IF
                oldParticle(contP) = newParticle(contP)

                CALL Exploration(oldParticle(contP), newParticle(contP), bestParticle(contP), opMpca, stMpca)
            ELSE
                CALL Scattering(oldParticle(contP), newParticle(contP), bestParticle(contP), opMpca, stMpca)
            END IF
        END DO

        ! BLACKBOARD UPDATE
        IF (((stMpca % NFE - stMpca % lastUpdate) > iCycleBlackboard) &
            .AND. (stMpca % higherNFE < opMpca % maxNFE / opMpca % nProcessors) &
            .AND. ((stMpca % higherNFE + iCycleBlackboard) < opMpca % maxNFE / opMpca % nProcessors)) &
            THEN
            DO contP = 1, opMpca % nParticlesProcessor
                if (bestParticle(contP) % fitness < bestParticleProcessor % fitness) then
                    bestParticleProcessor = bestParticle(contP)
                end if
            END DO

            CALL blackboard(bestParticleProcessor, opMpca, stMpca, &
            stMpca % higherNFE, stMpca % totalNFE)
            DO contP = 1, opMpca % nParticlesProcessor
                bestParticle(contP) = bestParticleProcessor
            END DO
            stMpca % lastUpdate = stMpca % NFE
        END IF

        contD = contD + 1

        ! Showing results
        if(opMpca % verbose .eqv. .true.) then
            write(*, '(I2)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(1))
            write(*, '(I3)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(2))
            write(*, '(I3)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(3))
            write(*, '(I2)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(4))
            write(*, '(S,ES10.2E2)', ADVANCE = 'NO') bestParticleProcessor % solution(5)
            write(*, '(S,ES10.2E2)', ADVANCE = 'NO') bestParticleProcessor % solution(6)
            write(*, '(S,ES10.2E2)', ADVANCE = 'NO') bestParticleProcessor % fitness
            write(*, *) stMpca % totalNFE
        end if
    END DO

    !*************************
    ! FINAL BLACKBOARD UPDATE 
    !************************ 
    DO contP = 1, opMpca % nParticlesProcessor
        if (bestParticle(contP) % fitness < bestParticleProcessor % fitness) then
            bestParticleProcessor = bestParticle(contP)
        end if
    END DO
    CALL blackboard(bestParticleProcessor, opMpca, stMpca, &
    stMpca % higherNFE, stMpca % totalNFE)
    DO contP = 1, opMpca % nParticlesProcessor
        bestParticle(contP) = bestParticleProcessor
    END DO

    if (opMpca % iProcessor == 0) then
        write(20, '(ES14.6E2)', ADVANCE = 'NO') bestParticleProcessor % fitness
        write(20, '(I2)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(1))
        write(20, '(I3)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(2))
        if (nint(bestParticleProcessor % solution(1)) == 2) then
            write(20, '(I3)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(3))
        else
            write(20, '(I3)', ADVANCE = 'NO') 0
        end if
        write(20, '(I2)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(4))
        write(20, '(ES14.6E2)', ADVANCE = 'NO') bestParticleProcessor % solution(5)
        write(20, '(ES14.6E2)') bestParticleProcessor % solution(6)

        if(opMpca % verbose .eqv. .true.) then
            write(*, '(I2)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(1))
            write(*, '(I3)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(2))
            write(*, '(I3)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(3))
            write(*, '(I2)', ADVANCE = 'NO') nint(bestParticleProcessor % solution(4))
            write(*, '(S,ES10.2E2)', ADVANCE = 'NO') bestParticleProcessor % solution(5)
            write(*, '(S,ES10.2E2)', ADVANCE = 'NO') bestParticleProcessor % solution(6)
            write(*, '(S,ES10.2E2)') bestParticleProcessor % fitness
        end if
    end if

    ! FINALIZING

    DEALLOCATE(vRandom)
    DEALLOCATE(oldParticle)
    DEALLOCATE(newParticle)
    DEALLOCATE(minB)
    DEALLOCATE(maxB)
    DEALLOCATE(bestParticle)
    DEALLOCATE(realSolution)


    CALL MPI_Finalize(iError)

    if (opMpca % saveEvolution .eqv. .true.) then
        CLOSE(10)
    end if
    if (opMpca % iProcessor == 0) then
        CLOSE(20)
    end if
    if (opMpca % saveSuccess .eqv. .true.) then
        CLOSE(30)
    end if

END PROGRAM MPCA
