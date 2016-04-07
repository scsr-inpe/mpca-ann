!***************************************************************!
! MPCA functions						!
!***************************************************************!
! Desenvolvido por: Eduardo Favero Pacheco da Luz (CAP/INPE)	!
! Modificado por: Reynier Hernández Torres (CAP/INPE)		!
! Baseado no PCA por Wagner F. Sacco				!
! Atualizacao: 03-Nov-2014					!
!***************************************************************!
MODULE mpcaFunctions
    USE nr
    USE nrtype
    USE ran_state, ONLY: ran_seed
    USE newTypes
    USE ModuloRNA
    USE Globais
    
CONTAINS
    SUBROUTINE Perturbation(oldParticle, newParticle, bestParticle, op, st)
        IMPLICIT NONE
        INTEGER :: contD
        REAL(SP) :: alea
        TYPE (Particle), INTENT(INOUT) :: oldParticle
        TYPE (Particle), INTENT(INOUT) :: newParticle
        TYPE (Particle), INTENT(INOUT) :: bestParticle
        TYPE (OptionsMPCA), INTENT(IN) :: op
        TYPE (StatusMPCA), INTENT(INOUT) :: st

        DO contD = 1, op % nDimensions
            CALL ran3(alea)
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
        
        newParticle % fitness = Rede_Neural_BP(newParticle % solution, op, st)
        st % NFE = st % NFE + 1
        
    END SUBROUTINE

    !********************************************************************
    SUBROUTINE Exploration(oldParticle, newParticle, bestParticle, op, st)
        IMPLICIT NONE
        INTEGER :: l2, nDimensions
        INTEGER*8 :: iteracoes
        TYPE (Particle), INTENT(INOUT) :: oldParticle, newParticle, bestParticle

        TYPE (OptionsMPCA), INTENT(IN) :: op
        TYPE (StatusMPCA), INTENT(INOUT) :: st

        DO l2 = 1, op % iterPerturbation
            CALL Small_Perturbation(oldParticle, newParticle, bestParticle, op, st)
            newParticle % fitness = Rede_Neural_BP(newParticle % solution, op, st)

            IF (newParticle % fitness .LT. oldParticle % fitness) THEN
                IF (newParticle % fitness .LT. bestParticle % fitness) THEN
                    bestParticle = newParticle
                    if (op % saveEvolution .eqv. .true.) then
                        write(10, '(S,ES18.10E4)', ADVANCE = 'NO') bestParticle % fitness
                        write(10, '(A8,I6)') ',E,', st % NFE
                    end if
                END IF
                oldParticle = newParticle
            END IF

            if ((bestParticle % fitness < op % emin) .and. (st % flag .eqv. .false.)) then
                if (op % saveSuccess .eqv. .true.) then
                    WRITE(30, '(S,ES18.10E4)', ADVANCE = 'NO') bestParticle % fitness
                    WRITE(30, '(A)', ADVANCE = 'NO') ','
                    WRITE(30, '(I6)', ADVANCE = 'NO') st % NFE
                    WRITE(30, '(A)', ADVANCE = 'NO') ','
                    WRITE(30, '(I6)', ADVANCE = 'NO') st % totalNFE
                    WRITE(30, '(A)', ADVANCE = 'NO') ','
                    WRITE(30, '(A)', ADVANCE = 'NO') trim(op % typeOpposition)
                    WRITE(30, '(A)', ADVANCE = 'NO') ','
                    WRITE(30, '(A)', ADVANCE = 'NO') trim(op % functionName)
                end if
                st % flag = .true.
            end if
        END DO
    END SUBROUTINE

    !********************************************************************
    SUBROUTINE Small_Perturbation(oldParticle, newParticle, bestParticle, op, st)
        IMPLICIT NONE
        INTEGER :: contD
        REAL*8, ALLOCATABLE, DIMENSION(:) :: inferior
        REAL*8, ALLOCATABLE, DIMENSION(:) :: superior
        REAL(SP), ALLOCATABLE, DIMENSION(:) :: alea
        REAL(SP) :: alea2
        TYPE (Particle), INTENT(INOUT) :: oldParticle, newParticle
        TYPE (Particle), INTENT(INOUT) :: bestParticle
        TYPE (OptionsMPCA), INTENT(IN) :: op
        TYPE (StatusMPCA), INTENT(INOUT) :: st

        ALLOCATE(inferior(op % nDimensions))
        ALLOCATE(superior(op % nDimensions))
        ALLOCATE(alea(3 * op % nDimensions))

        CALL ran3(alea)
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

        newParticle % fitness = Rede_Neural_BP(newParticle % solution, op, st)
        st % NFE = st % NFE + 1

        DEALLOCATE(inferior)
        DEALLOCATE(superior)
        DEALLOCATE(alea)
    END SUBROUTINE

    !********************************************************************
    SUBROUTINE Scattering(oldParticle, newParticle, bestParticle, op, st)
        IMPLICIT NONE
        INTEGER :: l2
        REAL*8 :: p_scat
        REAL(SP) :: alea
        REAL(SP), ALLOCATABLE, DIMENSION(:) :: harvest
        TYPE (Particle), INTENT(INOUT) :: oldParticle, newParticle, bestParticle
        TYPE (OptionsMPCA), INTENT(IN) :: op
        TYPE (StatusMPCA), INTENT(INOUT) :: st

        ALLOCATE(harvest(op % nDimensions))

        SELECT CASE (op % typeProbability)
        CASE (1) !Truncated exponential distribution
            p_scat = 1.0D0 - (bestParticle % fitness / newParticle % fitness)
        CASE (2) !Probability - Cauchy PDF
            p_scat = 1.0D0 - (1.0D0/(pi * 1.0D0 * (1.0D0 + ((newParticle % fitness - bestParticle % fitness) / 1.0D0)**2)))
        CASE (3) !Probability - Cauchy PDF complement
            p_scat = 1.0D0/(pi * 1.0D0 * (1.0D0 + ((newParticle % fitness - bestParticle % fitness)/1.0D0)**2))
        END SELECT

        CALL ran3(alea)
        IF (alea < p_scat) THEN
            CALL ran3(harvest)
            DO l2 = 1, op % nDimensions
                oldParticle % solution(l2) = (DBLE(harvest(l2))*(op % upperBound(l2) - op % lowerBound(l2))) &
                +op % lowerBound(l2)
            END DO
            oldParticle % fitness = Rede_Neural_BP(oldParticle % solution, op, st)
            st % NFE = st % NFE + 1

            IF (oldParticle % fitness .LT. bestParticle % fitness) THEN
                bestParticle = oldParticle
                if (op % saveEvolution .eqv. .true.) then
                    write(10, '(S,ES14.6E4)', ADVANCE = 'NO') bestParticle % fitness
                    write(10, *) ',*SCAT,', st % NFE
                end if
            END IF

            if ((bestParticle % fitness < op % emin) .and. (st % flag .eqv. .false.)) then
                if (op % saveSuccess .eqv. .true.) then
                    WRITE(30, '(S,ES18.10E4)', ADVANCE = 'NO') bestParticle % fitness
                    WRITE(30, '(A)', ADVANCE = 'NO') ','
                    WRITE(30, '(I6)', ADVANCE = 'NO') st % NFE
                    WRITE(30, '(A)', ADVANCE = 'NO') ','
                    WRITE(30, '(I6)', ADVANCE = 'NO') st % totalNFE
                    WRITE(30, '(A)', ADVANCE = 'NO') ','
                    WRITE(30, '(A)', ADVANCE = 'NO') trim(op % typeOpposition)
                    WRITE(30, '(A)', ADVANCE = 'NO') ','
                    WRITE(30, '(A)', ADVANCE = 'NO') trim(op % functionName)
                end if
                st % flag = .true.
            end if
        ELSE
            CALL Exploration(oldParticle, newParticle, bestParticle, op, st)
        END IF

        DEALLOCATE(harvest)
    END SUBROUTINE

    !********************************************************************
    SUBROUTINE blackboard(bestParticle, op, st, higherNFE, totalNFE)
        INCLUDE 'mpif.h'
        INTEGER :: ierr, status(MPI_STATUS_SIZE)
        INTEGER*8, INTENT(INOUT) :: higherNFE, totalNFE
        INTEGER :: contD, contP

        REAL*8, ALLOCATABLE, DIMENSION(:) :: sendParticle
        TYPE (Particle), INTENT(INOUT) :: bestParticle
        TYPE (OptionsMPCA), INTENT(IN) :: op
        TYPE (StatusMPCA), INTENT(INOUT) :: st

        ALLOCATE(sendParticle(op % nDimensions + 3))

        IF (op % iProcessor .EQ. 0) THEN
            higherNFE = st % NFE
            totalNFE = st % NFE

            DO contP = 1, op % nProcessors - 1
                CALL MPI_Recv(sendParticle, op % nDimensions + 3, MPI_DOUBLE_PRECISION, &
                MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status, ierr)

                IF (INT(sendParticle(op % nDimensions + 2)) .GT. higherNFE) THEN
                    higherNFE = INT(sendParticle(op % nDimensions + 2))
                END IF
                totalNFE = totalNFE + INT(sendParticle(op % nDimensions + 2))
                IF (sendParticle(op % nDimensions + 1) .LT. bestParticle % fitness) THEN
                    DO j = 1, op % nDimensions
                        bestParticle % solution(j) = sendParticle(j)
                    ENDDO
                    bestParticle % fitness = sendParticle(op % nDimensions + 1)
                ENDIF
            ENDDO

            DO contD = 1, op % nDimensions
                sendParticle(contD) = bestParticle % solution(contD)
            END DO
            sendParticle(op % nDimensions + 1) = bestParticle % fitness
            sendParticle(op % nDimensions + 2) = DBLE(higherNFE)
            sendParticle(op % nDimensions + 3) = DBLE(totalNFE)

            CALL MPI_Bcast(sendParticle, op % nDimensions + 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        ELSE
            DO contD = 1, op % nDimensions
                sendParticle(contD) = bestParticle % solution(contD)
            END DO
            sendParticle(op % nDimensions + 1) = bestParticle % fitness
            sendParticle(op % nDimensions + 2) = DBLE(st % NFE)

            CALL MPI_Send(sendParticle, op % nDimensions + 3, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)

            CALL MPI_Bcast(sendParticle, op % nDimensions + 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

            DO contD = 1, op % nDimensions
                bestParticle % solution(contD) = sendParticle(contD)
            END DO
            bestParticle % fitness = sendParticle(op % nDimensions + 1)
            higherNFE = INT(sendParticle(op % nDimensions + 2))
            totalNFE = INT(sendParticle(op % nDimensions + 3))
        END IF
    END SUBROUTINE
    
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
END MODULE mpcaFunctions