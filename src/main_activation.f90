PROGRAM MAIN_ACTIVATION
    USE annActivation
    USE newTypes

    implicit none

    TYPE(OptionsMPCA) :: op
    TYPE(annConfig) :: config
    REAL(kind = 8) :: MSE
    character (100) :: fString, dummy
    integer :: i
    integer :: j
    integer :: k

    op % nClasses = 96000
    op % nInputs = 2
    op % nOutputs = 1

    allocate(op % x(op % nInputs, op % nClasses))
    allocate(op % y(op % nOutputs, op % nClasses))

    fString = '(      F8.5)'
    write(fString(2:7), '(I6)') op % nClasses

    OPEN (2, file = './data/x_activation.txt')
    DO i = 1, op % nInputs
        READ(2, fString) (op % x(i, j), j = 1, op % nClasses)
    END DO
    CLOSE (2)
    
    OPEN (2, file = './data/y_activation.txt')
    DO i = 1, op % nInputs
        READ(2, fString) (op % y(i, j), j = 1, op % nClasses)
    END DO
    CLOSE (2)

    !--------------------------------------------------------------------!
    !WEIGHTS AND BIASES
    !--------------------------------------------------------------------!

    open(12, file = trim('./data/ann.best'), STATUS = "old")
    read(12, '(A)') dummy
    read(12, *) config % hiddenLayers
    read(12, *) config % neuronsLayer(1)

    if (config % hiddenLayers == 2) then
        read(12, *) config % neuronsLayer(2)
    end if
    
    ! Allocating space for config
    allocate(config % wh1(op % nInputs, config % neuronsLayer(1)))
    allocate(config % bh1(config % neuronsLayer(1)))
    if (config % hiddenLayers == 2) then
        allocate(config % wh2(config % neuronsLayer(1), config % neuronsLayer(2)))
        allocate(config % bh2(config % neuronsLayer(2)))
        allocate(config % ws(config % neuronsLayer(2), op % nOutputs))
    else
        allocate(config % wh2(1, 1))
        allocate(config % bh2(1))
        allocate(config % ws(config % neuronsLayer(1), op % nOutputs))
    end if
    allocate(config % bs(op % nOutputs))
    
    !Reading file
    read(12, '(A)') dummy
    fString = '(   F11.5)'
    write(fString(2:4), '(I3)') config % neuronsLayer(1)
    
    DO i = 1, op % nInputs
        read(12, fString) (config % wh1(i, k), k = 1, config % neuronsLayer(1))
    ENDDO

    read(12, '(A)') dummy
    read(12, fString) (config % bh1(k), k = 1, config % neuronsLayer(1))

    if (config % hiddenLayers == 2) then
        read(12, '(A)') dummy
        write(fString(2:4), '(I3)') config % neuronsLayer(2)
        DO i = 1, config % neuronsLayer(1)
            read(12, fString) (config % wh2(i, k), k = 1, config % neuronsLayer(2))
        ENDDO

        read(12, '(A)') dummy
        read(12, fString) (config % bh2(k), k = 1, config % neuronsLayer(2))

        write(fString(2:4), '(I3)') op % nOutputs
        read(12, '(A)') dummy
        DO i = 1, config % neuronsLayer(2)
            read(12, fString) (config % ws(i, k), k = 1, op % nOutputs)
        ENDDO
    else
        write(fString(2:4), '(I3)') op % nOutputs
        read(12, '(A)') dummy
        DO i = 1, config % neuronsLayer(1)
            read(12, fString) (config % ws(i, k), k = 1, op % nOutputs)
        ENDDO
    end if

    read(12, '(A)') dummy
    read(12, fString) (config % bs(k), k = 1, op % nOutputs)

    close(12)

    config % activationFunction = 2
    MSE = neuralNetworkActivation(config, op)

    deallocate(config % bh1)
    deallocate(config % bs)
    deallocate(config % wh1)
    deallocate(config % ws)
    if (config % hiddenLayers == 2) then
        deallocate(config % bh2)
        deallocate(config % wh2)
    end if
    
END PROGRAM MAIN_ACTIVATION