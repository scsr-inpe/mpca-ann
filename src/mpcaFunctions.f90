!***************************************************************!
! MPCA functions						!
!***************************************************************!
! Desenvolvido por: Eduardo Favero Pacheco da Luz (CAP/INPE)	!
! Modificado por: Reynier Hernández Torres (CAP/INPE)		!
! Baseado no PCA por Wagner F. Sacco				!
! Atualizacao: 03-Nov-2014					!
!***************************************************************!
MODULE mpcaFunctions

    USE newTypes
    USE annTraining
    
CONTAINS
    SUBROUTINE Perturbation(oldParticle, newParticle, bestParticle, op, st)
        IMPLICIT NONE
        INTEGER :: contD
        REAL (kind = 8) :: alea
        TYPE (Particle), INTENT(INOUT) :: oldParticle
        TYPE (Particle), INTENT(INOUT) :: newParticle
        TYPE (Particle), INTENT(INOUT) :: bestParticle
        TYPE (OptionsMPCA), INTENT(IN) :: op
        TYPE (StatusMPCA), INTENT(INOUT) :: st

        DO contD = 1, op % nDimensions
            CALL RANDOM_NUMBER(alea)
            newParticle % solution(contD) = oldParticle % solution(contD) &
            +((op % upperBound(contD) - oldParticle % solution(contD)) * alea) &
            -((oldParticle % solution(contD) - op % lowerBound(contD))*(1.0D0 - alea))

            IF (newParticle % solution(contD) .GT. op % upperBound(contD)) THEN
                newParticle % solution(contD) = op % upperBound(contD)
            END IF

            IF (newParticle % solution(contD) .LT. op % lowerBound(contD)) THEN
                newParticle % solution(contD) = op % lowerBound(contD)
            END IF
        END DO
        
        newParticle % fitness = neuralNetworkTraining(newParticle % solution, op, st)
        st % NFE = st % NFE + 1
        
        IF (newParticle % fitness .LT. bestParticle % fitness) THEN
            bestParticle = newParticle
        END IF
        
    END SUBROUTINE

    !********************************************************************
    SUBROUTINE Exploration(oldParticle, newParticle, bestParticle, op, st)
        IMPLICIT NONE
        INTEGER :: l2, nDimensions
        INTEGER (kind = 8) :: iteracoes
        TYPE (Particle), INTENT(INOUT) :: oldParticle, newParticle, bestParticle

        TYPE (OptionsMPCA), INTENT(IN) :: op
        TYPE (StatusMPCA), INTENT(INOUT) :: st

        DO l2 = 1, op % iterPerturbation
            CALL Small_Perturbation(oldParticle, newParticle, bestParticle, op, st)
            
            newParticle % fitness = neuralNetworkTraining(newParticle % solution, op, st)

            IF (newParticle % fitness .LT. oldParticle % fitness) THEN
                oldParticle = newParticle
            END IF
            
            IF (newParticle % fitness .LT. bestParticle % fitness) THEN
                bestParticle = newParticle
            END IF

            if ((bestParticle % fitness < op % emin) .and. (st % flag .eqv. .false.)) then
                st % flag = .true.
            end if
            
            if (st % NFE >= op % maxNFE / op % nProcessors) then
                st % doStop = .true.
                return
            end if
        END DO
    END SUBROUTINE

    !********************************************************************
    SUBROUTINE Small_Perturbation(oldParticle, newParticle, bestParticle, op, st)
        IMPLICIT NONE
        INTEGER :: contD
        REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: inferior
        REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: superior
        REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: alea
        REAL (kind = 8) :: alea2
        TYPE (Particle), INTENT(INOUT) :: oldParticle, newParticle
        TYPE (Particle), INTENT(INOUT) :: bestParticle
        TYPE (OptionsMPCA), INTENT(IN) :: op
        TYPE (StatusMPCA), INTENT(INOUT) :: st

        ALLOCATE(inferior(op % nDimensions))
        ALLOCATE(superior(op % nDimensions))
        ALLOCATE(alea(3 * op % nDimensions))

        CALL RANDOM_NUMBER(alea)
        DO contD = 1, op % nDimensions
            superior(contD) = (DBLE(alea(contD)) &
            *(op % up_small - 1.0D0) + 1.0D0) &
            *oldParticle % solution(contD)

            inferior(contD) = (DBLE(alea(contD + op % nDimensions)) &
            *(1.0D0 - op % lo_small) + op % lo_small) &
            *oldParticle % solution(contD)

            IF (superior(contD) .GT. op % upperBound(contD)) THEN
                superior(contD) = op % upperBound(contD)
            END IF

            IF (inferior(contD) .LT. op % lowerBound(contD)) THEN
                inferior(contD) = op % lowerBound(contD)
            END IF

            newParticle % solution(contD) = oldParticle % solution(contD) &
            +((superior(contD) - oldParticle % solution(contD)) * DBLE(alea(contD + 2 * op % nDimensions))) &
            -((oldParticle % solution(contD) - inferior(contD))*(1.0D0 - DBLE(alea(contD + 2 * op % nDimensions))))

            IF (newParticle % solution(contD) .GT. op % upperBound(contD)) THEN
                newParticle % solution(contD) = op % upperBound(contD)
            END IF

            IF (newParticle % solution(contD) .LT. op % lowerBound(contD)) THEN
                newParticle % solution(contD) = op % lowerBound(contD)
            END IF
        END DO
        
        newParticle % fitness = neuralNetworkTraining(newParticle % solution, op, st)
        st % NFE = st % NFE + 1

        DEALLOCATE(inferior)
        DEALLOCATE(superior)
        DEALLOCATE(alea)
    END SUBROUTINE

    !********************************************************************
    SUBROUTINE Scattering(oldParticle, newParticle, bestParticle, op, st)
        IMPLICIT NONE
        INTEGER :: l2
        REAL (kind = 8) :: p_scat
        REAL (kind = 8) :: alea
        TYPE (Particle), INTENT(INOUT) :: oldParticle, newParticle, bestParticle
        TYPE (OptionsMPCA), INTENT(IN) :: op
        TYPE (StatusMPCA), INTENT(INOUT) :: st
        REAL  (kind = 8), PARAMETER :: pi = 3.1415927

        SELECT CASE (op % typeProbability)
        CASE (1) !Truncated exponential distribution
            p_scat = 1.0D0 - (bestParticle % fitness / newParticle % fitness)
        CASE (2) !Probability - Cauchy PDF
            p_scat = 1.0D0 - (1.0D0/(pi * 1.0D0 * (1.0D0 + ((newParticle % fitness - bestParticle % fitness) / 1.0D0)**2)))
        CASE (3) !Probability - Cauchy PDF complement
            p_scat = 1.0D0/(pi * 1.0D0 * (1.0D0 + ((newParticle % fitness - bestParticle % fitness)/1.0D0)**2))
        END SELECT

        CALL RANDOM_NUMBER(alea)
        IF (alea < p_scat) THEN
            DO l2 = 1, op % nDimensions
                CALL RANDOM_NUMBER(alea)
                oldParticle % solution(l2) = (alea * (op % upperBound(l2) - op % lowerBound(l2))) &
                +op % lowerBound(l2)
            END DO
            
            oldParticle % fitness = neuralNetworkTraining(oldParticle % solution, op, st)
            st % NFE = st % NFE + 1

            IF (oldParticle % fitness .LT. bestParticle % fitness) THEN
                bestParticle = oldParticle
            END IF

            if ((bestParticle % fitness < op % emin) .and. (st % flag .eqv. .false.)) then
                st % flag = .true.
            end if
        END IF
        
        CALL Exploration(oldParticle, newParticle, bestParticle, op, st)

    END SUBROUTINE

    !********************************************************************
    !********************************************************************
    !********************************************************************
    !********************************************************************
    SUBROUTINE blackboard(bo, NFE, higherNFE, totalNFE, doStop, doStopMPCA, op)
    USE newtypes

    implicit none
    INCLUDE 'mpif.h'

    INTEGER :: i, iBest, j, ierr, status(MPI_STATUS_SIZE), stopCount
    INTEGER (kind=8), INTENT(in) :: NFE
    INTEGER (kind=8), INTENT(inout) :: higherNFE, totalNFE
    INTEGER :: world_rank, world_size, nDimensions
    logical, intent(inout) :: doStop, doStopMPCA

    type (Particle), intent(inout) :: bo
    TYPE(OptionsMPCA), intent(in) :: op
    real (kind = 8), allocatable, dimension(:) :: send

    world_rank = op % iProcessor
    world_size = op % nProcessors
    nDimensions = op % nDimensions
    
    if (world_size == 1) then
        call copyFileBest(op, 0)
        higherNFE = NFE
        totalNFE = NFE
        doStopMPCA = doStop
        return
    end if
    
    ALLOCATE(send(nDimensions + 4))

    stopCount = 0
    IF (world_rank .EQ. 0) THEN
        higherNFE = NFE
        totalNFE = NFE

        iBest = 0
        DO i = 1, world_size - 1
            CALL MPI_Recv(send, nDimensions + 4, MPI_REAL8, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            IF (INT(send(nDimensions + 2)) .GT. higherNFE) THEN
                higherNFE = INT8(send(nDimensions + 2))
            END IF

            totalNFE = totalNFE + INT8(send(nDimensions + 2))

            IF (send(nDimensions + 1) < bo % fitness) THEN
                DO j = 1, nDimensions
                    bo % solution(j) = send(j)
                ENDDO
                bo % fitness = send(nDimensions + 1)
                iBest = i
            ENDIF

            if (send(nDimensions + 3) > 0) then
                stopCount = stopCount + 1;
            end if
        ENDDO

        DO i = 1, nDimensions
            send(i) = bo % solution(i)
        END DO

        send(nDimensions + 1) = bo % fitness
        send(nDimensions + 2) = DFLOAT(higherNFE)
        send(nDimensions + 3) = DFLOAT(totalNFE)
        send(nDimensions + 4) = DFLOAT(stopCount)

        CALL MPI_Bcast(send, nDimensions + 4, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        if (send(nDimensions + 4) > 0.9) then
            doStopMPCA = .true.
        else
            doStopMPCA = .false.
        end if
        
        call copyFileBest(op, iBest)

    ELSE
        DO i = 1, nDimensions
            send(i) = bo % solution(i)
        END DO

        send(nDimensions + 1) = bo % fitness
        send(nDimensions + 2) = DFLOAT(NFE)
        
        if (doStop .eqv. .true.) then
            send(nDimensions + 3) = DFLOAT(1);
        else
            send(nDimensions + 3) = DFLOAT(0);
        end if
        
        send(nDimensions + 4) = DFLOAT(0);

        CALL MPI_Send(send, nDimensions + 4, MPI_REAL8, 0, 0, MPI_COMM_WORLD, ierr)

        CALL MPI_Bcast(send, nDimensions + 4, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        bo % solution = send(1:nDimensions)
        bo % fitness = send(nDimensions + 1)
        higherNFE = INT8(send(nDimensions + 2))
        totalNFE = INT8(send(nDimensions + 3))

        if (send(nDimensions + 4) > 0.9) then
        doStopMPCA = .true.
    else
        doStopMPCA = .false.
        end if

    END IF

    deallocate(send)

END SUBROUTINE blackboard

!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
    
    !********************************************************************
        function best_nearby(delta, point, prevbest, nvars, f, funevals, lb, ub)
        
        implicit none

        integer ( kind = 4) nvars

        real ( kind = 8) best_nearby
        real ( kind = 8) delta(nvars)
        real ( kind = 8), external :: f
        real ( kind = 8) ftmp
        integer ( kind = 8) funevals
        integer ( kind = 8) i
        real ( kind = 8) minf
        real ( kind = 8) point(nvars)
        real ( kind = 8) prevbest
        real ( kind = 8) z(nvars)
        real ( kind = 8) lb(nvars)
        real ( kind = 8) ub(nvars)

        minf = prevbest
        z(1:nvars) = point(1:nvars)

        do i = 1, nvars
            z(i) = point(i) + delta(i)

            if (z(i) .lt. lb(i)) then
                z(i) = lb(i)
            elseif (z(i) .gt. ub(i)) then
                z(i) = ub(i)
            end if

            ftmp = f(z, nvars)
            funevals = funevals + 1
            
            
            if (ftmp < minf) then
                minf = ftmp
            else
                delta(i) = -delta(i)
                z(i) = point(i) + delta(i)
                
                if (z(i) .lt. lb(i)) then
                    z(i) = lb(i)
                elseif (z(i) .gt. ub(i)) then
                    z(i) = ub(i)
                end if

                ftmp = f(z, nvars)
                funevals = funevals + 1

                if (ftmp < minf) then
                    minf = ftmp
                else
                    z(i) = point(i)
                end if
            end if
        end do

        point(1:nvars) = z(1:nvars)
        best_nearby = minf

        return
    end function best_nearby
    
    !*****************************************************************
    function hooke(nvars, startpt, endpt, rho, eps, nfemax, f, lb, ub)
        implicit none

        integer ( kind = 4) nvars

        real ( kind = 8) delta(nvars)
        real ( kind = 8) endpt(nvars)
        real ( kind = 8) eps
        real ( kind = 8), external :: f
        real ( kind = 8) fbefore
        integer ( kind = 8) funevals
        integer ( kind = 8) hooke
        integer ( kind = 4) i
        integer ( kind = 8) nfemax
        integer ( kind = 4) iters
        integer ( kind = 4) j
        integer ( kind = 4) keep
        real ( kind = 8) newf
        real ( kind = 8) newx(nvars)
        real ( kind = 8) rho
        real ( kind = 8) startpt(nvars)
        real ( kind = 8) lb(nvars)
        real ( kind = 8) ub(nvars)
        real ( kind = 8) steplength
        real ( kind = 8) tmp
        logical, parameter :: verbose = .false.
        real ( kind = 8) xbefore(nvars)

        newx(1:nvars) = startpt(1:nvars)
        xbefore(1:nvars) = startpt(1:nvars)

        do i = 1, nvars
            if (startpt(i) == 0.0D+00) then
                delta(i) = rho
            else
                delta(i) = rho * abs(startpt(i))
            end if
        end do

        funevals = 0
        steplength = rho
        iters = 0
        fbefore = f(newx)
        funevals = funevals + 1
        newf = fbefore

        do while(funevals < nfemax .and. eps < steplength)

            iters = iters + 1

            if (verbose) then

                write ( *, '(a)') ' '
                write ( *, '(a,i8,a,ES12.2E4)') &
                '  FUNEVALS, = ', funevals, '  F(X) = ', fbefore

                !do j = 1, nvars
                !    write ( *, '(2x,i8,2x,g14.6)') j, xbefore(j)
                !end do
            end if
            !
            !  Find best new point, one coordinate at a time.
            !
            newx(1:nvars) = xbefore(1:nvars)

            newf = best_nearby(delta, newx, fbefore, nvars, f, funevals, lb, ub)
            !
            !  If we made some improvements, pursue that direction.
            !
            keep = 1

            do while (newf < fbefore .and. keep == 1)
                do i = 1, nvars
                    !
                    !  Arrange the sign of DELTA.
                    !
                    if (newx(i) <= xbefore(i)) then
                        delta(i) = -abs(delta(i))
                    else
                        delta(i) = abs(delta(i))
                    end if
                    !
                    !  Now, move further in this direction.
                    !
                    tmp = xbefore(i)
                    xbefore(i) = newx(i)
                    newx(i) = newx(i) + newx(i) - tmp
                    
                    if (newx(i) .lt. lb(i)) then
                        newx(i) = lb(i)
                    elseif (newx(i) .gt. ub(i)) then
                        newx(i) = ub(i)
                    end if
                        
                end do

                fbefore = newf
                newf = best_nearby(delta, newx, fbefore, nvars, f, funevals, lb, ub)
                !
                !  If the further (optimistic) move was bad...
                !
                if (fbefore <= newf) then
                    exit
                end if
                !
                !  Make sure that the differences between the new and the old points
                !  are due to actual displacements; beware of roundoff errors that
                !  might cause NEWF < FBEFORE.
                !
                keep = 0

                do i = 1, nvars
                    if (0.5D+00 * abs(delta(i)) < &
                        abs(newx(i) - xbefore(i))) then
                        keep = 1
                        exit
                    end if
                end do

            end do

            if (eps <= steplength .and. fbefore <= newf) then
                steplength = steplength * rho
                delta(1:nvars) = delta(1:nvars) * rho
            end if

        end do

        endpt(1:nvars) = xbefore(1:nvars)

        hooke = funevals

        return
    end function hooke
    
    
    subroutine copyFileBest(op, iBest)
        implicit none
    
        character (100) :: str0, str1, filename, command
        integer :: iBest
        TYPE(OptionsMPCA), intent(in) :: op
    
        IF (iBest < 9) THEN
            WRITE (str0, '(I1)') iBest + 1
        ELSE IF (iBest < 99) THEN
            WRITE (str0, '(I2)') iBest + 1
        ELSE
            WRITE (str0, '(I3)') iBest + 1
        END IF
    
        IF (op % iExperiment < 10) THEN
            WRITE (str1, '(I1)') op % iExperiment
        ELSE IF (op % iExperiment < 100) THEN
            WRITE (str1, '(I2)') op % iExperiment
        ELSE
            WRITE (str1, '(I3)') op % iExperiment
        END IF

        filename = './output/ann' // trim(str1) // '_' // trim(str0) // '.out'
        command = 'cp ' // TRIM(filename) // ' ./output/ann' // trim(str1) // '.best'
        call system(TRIM(command))
        
    end subroutine copyFileBest
END MODULE mpcaFunctions