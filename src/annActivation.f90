!*******************************************************************!
! Optimization of the architecture of Artificial Neural Networks    !
! using Multi-Particle Collision Algorithm                          !
!*******************************************************************!
! Developed by: Juliana Anochi and Sabrina Sambatti (CAP/INPE)      !
! Modified by: Reynier Hernandez Torres (CAP/INPE)                  !
! Updated: 11-Mai-2016                                              !
!*******************************************************************!

MODULE annActivation
    use newTypes
    use foul

CONTAINS

    REAL(kind = 8) FUNCTION neuralNetworkActivation(config, op)

        IMPLICIT NONE

        TYPE(OptionsMPCA), intent(in) :: op
        TYPE(annConfig), intent(in) :: config
        
        real (8) :: BestFO
        real (8) :: rNumber
        real (8) :: dOutput
        real (8) :: penaltyObj
        real (8) :: aux
        real (8), allocatable, dimension(:) :: vs
        real (8), allocatable, dimension(:,:) :: error
        real (8), allocatable, dimension(:) :: vh1
        real (8), allocatable, dimension(:) :: vh2
        real (8), allocatable, dimension(:) :: errorClass
        real (8), allocatable, dimension(:) :: yh1
        real (8), allocatable, dimension(:) :: yh2
        real (8), allocatable, dimension(:,:) :: y


        real (8), parameter :: alphaObj = 0.1D0
        real (8), parameter :: betaObj = 1.0D0
        real (8), parameter :: p1 = 5.0e-8
        real (8), parameter :: p2 = 5.0e-5

        character (100) :: str0
        character (100) :: str1
        character (100) :: fString
        character (100) :: dummy

        integer :: i
        integer :: j
        integer :: k

        !Allocating space for error variables
        allocate(error(op % nOutputs, op % nClasses))
        allocate(errorClass(op % nClasses))

        ! Allocating space for v and y variables
        allocate(vh1(config % neuronsLayer(1)))
        allocate(yh1(config % neuronsLayer(1)))
        if (config % hiddenLayers == 2) then
            allocate(vh2(config % neuronsLayer(2)))
            allocate(yh2(config % neuronsLayer(2)))
        else
            allocate(vh2(1))
            allocate(yh2(1))
        end if
        allocate(vs(op % nOutputs))

        neuralNetworkActivation = 0

        !----------------------------------------------------------------------!
        ! ACTIVATION
        !----------------------------------------------------------------------!        
        allocate(y(op % nOutputs, op % nClasses))
        
        fString = '(      F8.5)'
        write(fString(2:7), '(I6)') op % nClasses

        do i = 1, op % nClasses
            ! ACTIVATING HIDDEN LAYER 1
            vh1 = matmul(op % x(:, i), config % wh1) - config % bh1
            fString = '(      F8.5)'
            write(fString(2:7), '(I6)') config % neuronsLayer(1)            
            
            yh1 = activation(vh1, config % activationFunction)

            if (config % hiddenLayers == 1) then
                vs = matmul(yh1, config % ws) - config % bs
            end if

            ! ACTIVATING HIDDEN LAYER 2
            if (config % hiddenLayers == 2) then
                vh2 = matmul(yh1, config % wh2) - config % bh2
                yh2 = activation(vh2, config % activationFunction)
                vs = matmul(yh2, config % ws) - config % bs
            end if

            ! ACTIVATING OUTPUT
            y(:, i) = activation(vs, config % activationFunction)
            error(:, i) = op % y(:, i) - y(:, i)
        end do
        
        neuralNetworkActivation = sum(error) / dfloat(op % nClasses)
    
        fString = '(      F8.5)'
        OPEN (2, file = './output/y_activation.txt')
        DO i = 1, op % nOutputs
            write(2, fString) (y(i, j) , j = 1, op % nClasses)
        END DO
        CLOSE (2)

        deallocate(error)
        deallocate(errorClass)
        deallocate(vh1)
        deallocate(yh1)
        if (config % hiddenLayers == 2) then
            deallocate(vh2)
            deallocate(yh2)
        end if
        deallocate(vs)
        deallocate(y)

    END FUNCTION neuralNetworkActivation


    real(kind = 8) elemental function activation(v, afunction) result(y)
        implicit none

        real (kind = 8), intent(in) :: v
        integer, intent(in) :: afunction

        real (kind = 8), parameter :: a = 1

        select case(afunction)
        case (1) !LOGISTIC
            y = 1.d0/(1.d0 + exp(-a * v))
        case (2) !TANGENT
            y = (1.d0 - exp(-v))/(1.d0 + exp(-v))
        case (3) !GAUSS
            y = exp(-v)
        end select
    end function activation

    real(kind = 8) elemental function derivate(v, y, afunction) result(d)
        implicit none
        real (kind = 8), intent(in) :: v
        real (kind = 8), intent(in) :: y
        integer, intent(in) :: afunction
        real (kind = 8), parameter :: a = 1

        select case(afunction)
        case (1)!LOGISTICA: e^-x / (1+e^-x)^2
            d = ((a * exp(-a * v))/((1.d0 + exp(-a * v))**2.d0))
        case (2)!TANGENTE: 2e^-x / (1+e^-x)^2
            d = (2 * exp(-v)) / ((1 + exp(-v))**2)
        case (3)!GAUSS
            d = -y/a
        end select
    end function derivate

END MODULE annActivation