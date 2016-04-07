!***************************************************************!
! Optimization of ANN architecture by metaheuristic MPCA        !
!***************************************************************!
! Developed by: Juliana Anochi and Sabrina Sambatti (CAP/INPE)  !
! Modified by: Reynier Hernandez Torres (CAP/INPE)              !
! Updated: 24-Mar-2016                                          !
!***************************************************************!
MODULE ModuloRNA
    USE Globais
    USE newTypes

CONTAINS

    !Funcao para treinar a rede neural com os solution calculados pela heuristica
    DOUBLE PRECISION FUNCTION Rede_Neural_BP(solution, op, st)

        IMPLICIT NONE

        DOUBLE PRECISION :: solution(:)

        TYPE(OptionsMPCA), intent(in) :: op
        TYPE(StatusMPCA), intent(inout) :: st

        IF (nint(solution(1)) .EQ. 1) THEN
            Rede_Neural_BP = Rede_Neural_C1(solution, op, st)
        ELSE IF (nint(solution(1)) .EQ. 2) THEN
            Rede_Neural_BP = Rede_Neural_C2(solution, op, st)
        ENDIF

    END FUNCTION Rede_Neural_BP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DOUBLE PRECISION FUNCTION Rede_Neural_C1(solution, op, st)

        IMPLICIT NONE

        !Parametros recebidos pela funcao (passados pelo algoritmo de otimizacao)
        double precision :: solution(:)
        character*10 :: str0, str1

        !Variaveis usadas para o calculo da funcao objetivo
        double precision :: BestFO
        double precision :: penaltyObj, alphaObj, betaObj, p1, p2
        double precision, allocatable, dimension(:,:) :: x
        double precision, allocatable, dimension(:,:) :: yd
        double precision, allocatable, dimension(:,:) :: x_valid
        double precision, allocatable, dimension(:,:) :: yd_valid
        double precision, allocatable, dimension(:) :: vs
        double precision, allocatable, dimension(:,:) :: ys, ys_menor
        double precision, allocatable, dimension(:,:) :: ws, ws_velho, ws_menor
        double precision, allocatable, dimension(:) :: bs, bs_velho, bs_menor, bs_melhor
        double precision, allocatable, dimension(:,:) :: eqm_menor, eqm_valid_menor
        double precision, allocatable, dimension(:,:) :: erro
        double precision, allocatable, dimension(:) :: erro_pad
        double precision, allocatable, dimension(:,:) :: eqm
        double precision, allocatable, dimension(:,:) :: eqm_valid
        double precision, allocatable, dimension(:) :: grad_saida, grad_hide1
        double precision, allocatable, dimension(:,:) :: deltaw_saida
        double precision, allocatable, dimension(:) :: deltabs, deltabh1
        double precision, allocatable, dimension(:) :: vh1, vh2
        double precision, allocatable, dimension(:,:) :: wh1, wh1_velho, wh1_menor
        double precision, allocatable, dimension(:) :: bh1, bh1_velho, bh1_menor, bh1_melhor
        double precision, allocatable, dimension(:) :: yh1
        double precision, allocatable, dimension(:,:) :: deltaw_h1

        double precision, parameter :: a = 1

        double precision :: dv, soma, erro_menor, alpha
        double precision :: eta ! taxa de aprendizado
        double precision :: ed ! erro desejado
        double precision :: aleat ! gerar pesos e bias aleatoriamente

        integer :: hide_neuron(2) ! numero de neuronios na camada escondida
        integer :: cont
        integer :: l
        integer :: max_it ! numero maximo de iteracoes
        integer :: vetor_entrada ! dimensao do vetor de entrada
        integer :: num_pad ! numero de padroes a serem treinados
        integer :: num_pad_valid
        integer :: vetor_saida ! dimensao do vetor de saida
        integer :: i, j, k, hn ! variavel para controle do programa
        integer :: f_ativa ! controle de qual funcao ativacao
        integer :: f_deriva ! controle de qual funcao derivada
        integer :: hide_camada ! numero de camadas escondidas
        logical :: peso ! controlar como os pesos serao inicializados 
        logical :: validacao ! controlar se terah ou nao validacao

        TYPE(OptionsMPCA), intent(in) :: op
        TYPE(StatusMPCA), intent(inout) :: st

        !Atribuicao empirica dos pesos para a funcao objetivo
        alphaObj = 0.1D0
        betaObj = 1.0D0
        p1 = 5.0000E-08 !5*10^-8 
        p2 = 5.0000e-05 !5*10^-5
        erro_menor = 0.1D0

        hide_camada = nint(solution(1))
        hide_neuron(1) = nint(solution(2))
        hide_neuron(2) = nint(solution(3))
        f_ativa = nint(solution(4))
        alpha = solution(5)
        eta = solution(6)

        num_pad = op % nClasses
        num_pad_valid = op % nClassesValidation
        vetor_entrada = op % nInputs
        vetor_saida = op % nOutputs
        ed = op % targetError
        max_it = op % nEpochs
        peso = op % randomWeights
        validacao = op % haveValidation
        BestFO = st % bestObjectiveFunction
        f_deriva = f_ativa

        !Bias 1. camada oculta
        allocate(bh1(hide_neuron(1)))
        allocate(bh1_menor(hide_neuron(1)))
        allocate(bh1_melhor(hide_neuron(1)))
        allocate(bh1_velho(hide_neuron(1)))
        !Bias camada de saida
        allocate(bs(vetor_saida))
        allocate(bs_menor(vetor_saida))
        allocate(bs_melhor(vetor_saida))
        allocate(bs_velho(vetor_saida))
        !DeltaBias
        allocate(deltabh1(hide_neuron(1)))
        allocate(deltabs(vetor_saida))
        !DeltaW's
        allocate(deltaw_h1(vetor_entrada, hide_neuron(1)))
        allocate(deltaw_saida(hide_neuron(hide_camada), vetor_saida))
        !Erros
        allocate(eqm(1, max_it))
        allocate(eqm_valid(1, max_it))
        allocate(eqm_valid_menor(1, max_it))
        allocate(eqm_menor(1, max_it))
        allocate(erro(vetor_saida, num_pad))
        allocate(erro_pad(num_pad))
        !Gradientes
        allocate(grad_hide1(hide_neuron(1)))
        allocate(grad_saida(vetor_saida))
        !Campo local induzido
        allocate(vh1(hide_neuron(1)))
        allocate(vs(vetor_saida))
        !Pesos 1. camada oculta
        allocate(wh1(vetor_entrada, hide_neuron(1)))
        allocate(wh1_menor(vetor_entrada, hide_neuron(1)))
        allocate(wh1_velho(vetor_entrada, hide_neuron(1)))
        !Pesos camada de saida
        allocate(ws(hide_neuron(hide_camada), vetor_saida))
        allocate(ws_menor(hide_neuron(hide_camada), vetor_saida))
        allocate(ws_velho(hide_neuron(hide_camada), vetor_saida))
        !Dados de entrada e saida
        allocate(x(vetor_entrada, num_pad))
        allocate(x_valid(vetor_entrada, num_pad_valid))
        allocate(yd(vetor_saida, num_pad))
        allocate(yd_valid(vetor_saida, num_pad_valid))
        !Saidas obtidas
        allocate(yh1(hide_neuron(1)))
        allocate(ys(vetor_saida, num_pad))
        allocate(ys_menor(vetor_saida, num_pad))

        !------------------------------------------------------------!
        !LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
        !------------------------------------------------------------!
        OPEN (1, file = op % path_sai)
        DO I = 1, vetor_saida
            READ(1, *) (yd(I, J), J = 1, num_pad)
        END DO
        CLOSE (1)

        OPEN (2, file = op % path_ent)
        DO I = 1, vetor_entrada
            READ(2, *) (x(I, J), J = 1, num_pad)
        END DO
        CLOSE (2)

        IF (validacao .eqv. .true.) THEN
            OPEN (1, file = op % path_sai_valid)
            DO I = 1, vetor_saida
                READ(1, *) (yd_valid(I, J), J = 1, num_pad_valid)
            END DO
            CLOSE (1)

            OPEN (2, file = op % path_ent_valid)
            DO I = 1, vetor_entrada
                READ(2, *) (x_valid(I, J), J = 1, num_pad_valid)
            END DO
            CLOSE (2)
        ENDIF
        !--------------------------------------------------------------------!
        !INICIALIZANDO OS PARAMETROS: wh1, bh1, ws e bs
        !--------------------------------------------------------------------!

        !PRIMEIRA CAMADA OCULTA
        DO l = 1, vetor_entrada
            DO k = 1, hide_neuron(1)
                if (peso .eqv. .false.) then
                    wh1(l, k) = 0.5
                else
                    call random_number(aleat)
                    wh1(l, k) = aleat
                endif
            ENDDO
        ENDDO
        DO k = 1, hide_neuron(1)
            if (peso .eqv. .false.) then
                bh1(k) = 0.5
            else
                call random_number(aleat)
                bh1(k) = aleat
            endif
        ENDDO

        !CAMADA DE SAIDA
        DO l = 1, hide_neuron(hide_camada)
            DO k = 1, vetor_saida
                if (peso .eqv. .false.) then
                    ws(l, k) = 0.5
                else
                    call random_number(aleat)
                    ws(l, k) = aleat
                endif
            ENDDO
        ENDDO
        DO k = 1, vetor_saida
            if (peso .eqv. .false.) then
                bs(k) = 0.5
            else
                call random_number(aleat)
                bs(k) = aleat
            endif
        ENDDO

        !----------------------------------------------------------------------!
        ! INICIO DA REDE: FEEDFORWARD
        !----------------------------------------------------------------------!
        cont = 1 !Contador de epocas
        l = 0

        !DO WHILE ((erro_menor .GE. ed) .AND. (l .LT. max_it))
        DO WHILE (l .LT. max_it)

            l = l + 1

            DO i = 1, num_pad !DO NUMERO TOTAL DE PADROES

                ! ATIVACAO: 1. CAMADA OCULTA
                vh1 = 0.d0
                vh1(:) = matmul(x(:, i), wh1(:,:))
                vh1(:) = vh1(:) - bh1(:);
                select case(f_ativa)
                case (1) !LOGISTICA
                    yh1(:) = 1.d0/(1.d0 + DEXP(-a * vh1(:)))
                case (2) !TANGENTE
                    yh1(:) = (1.d0 - DEXP(-vh1(:)))/(1.d0 + DEXP(-vh1(:)))
                case (3) !GAUSS
                    yh1(:) = DEXP(-vh1(:))
                end select

                ! ATIVACAO: CAMADA DE SAIDA
                vs = 0.d0
                vs(:) = matmul(yh1(:), ws(:,:))
                vs(:) = vs(:) - bs(:)

                select case(f_ativa)
                case (1) !LOGISTICA
                    ys(:, i) = 1.d0/(1.d0 + DEXP(-a * vs(:)))
                case (2) !TANGENTE
                    ys(:, i) = (1.d0 - DEXP(-vs(:)))/(1.d0 + DEXP(-vs(:)))
                case (3) !GAUSS
                    ys(:, i) = DEXP(-vs(:))
                end select

                ! PARA O CALCULO DO NOVO PESO
                wh1_velho = wh1
                bh1_velho = bh1
                ws_velho = ws
                bs_velho = bs

                ! CALCULO ERRO TREINAMENTO
                erro(:, i) = yd(:, i) - ys(:, i)

                !-------------------------------------------------------------------------!
                !                        BACKPROPAGATION
                !-------------------------------------------------------------------------!
                ! TREINAMENTO: CAMADA DE SAIDA
                DO j = 1, vetor_saida

                    select case(f_deriva)
                    case (1)!LOGISTICA: e^-x / (1+e^-x)^2
                        dv = ((a * DEXP(-a * vs(j)))/((1.d0 + DEXP(-a * vs(j)))**2.d0))
                    case (2)!TANGENTE: 2e^-x / (1+e^-x)^2
                        dv = (2 * DEXP(-vs(j)))/ ((1 + DEXP(-vs(j)))**2)
                    case (3)!GAUSS
                        dv = -ys(j, i)/a
                    end select

                    grad_saida(j) = erro(j, i) * dv
                    deltaw_saida(:, j) = eta * grad_saida(j) * yh1(:)
                    ws(:, j) = ws(:, j) + alpha * (ws(:, j) - ws_velho(:, j)) + deltaw_saida(:, j)
                    deltabs(j) = eta * grad_saida(j) * (-1.d0)
                    !deltabs(j) = eta * grad_saida(j) 
                    bs(j) = bs(j) + deltabs(j)

                enddo !CAMADA SAIDA

                ! TREINAMENTO: 1. CAMADA OCULTA
                DO j = 1, hide_neuron(1)

                    soma = 0.d0
                    DO k = 1, vetor_saida
                        soma = soma + (grad_saida(k) * ws_velho(j, k))
                    enddo

                    select case(f_deriva)
                    case (1)!LOGISTICA: e^-x / (1+e^-x)^2
                        dv = ((a * DEXP(-a * vh1(j)))/((1.d0 + DEXP(-a * vh1(j)))**2.d0))
                    case (2)!TANGENTE
                        dv = (2 * DEXP(-vh1(j)))/ ((1 + DEXP(-vh1(j)))**2)
                    case (3)!GAUSS
                        dv = -yh1(j)/a
                    end select

                    grad_hide1(j) = dv * soma
                    deltaw_h1(:, j) = eta * grad_hide1(j) * x(:, i)
                    wh1(:, j) = wh1(:, j) + alpha * (wh1(:, j) - wh1_velho(:, j)) + deltaw_h1(:, j)
                    deltabh1(j) = eta * grad_hide1(j) * (-1.d0)
                    bh1(j) = bh1(j) + deltabh1(j)

                ENDDO ! 1. CAMADA OCULTA

                !CALCULO PADRAO DO ERRO
                erro_pad(i) = sum(erro(:, i), dim = 1)
                erro_pad(i) = 0.5d0 * (erro_pad(i)**2.d0)

            ENDDO !NUMERO PADROES
            eqm(1, l) = sum(erro_pad(:))
            eqm(1, l) = (1.d0/(num_pad)) * eqm(1, l)

            !***********************************************************
            ! VALIDACAO
            !***********************************************************

            DO i = 1, num_pad_valid

                ! ATIVACAO: 1. CAMADA OCULTA
                vh1 = 0.d0
                vh1(:) = matmul(x_valid(:, i), wh1(:,:))
                vh1(:) = vh1(:) - bh1(:);

                select case(f_ativa)
                case (1) !LOGISTICA
                    yh1(:) = 1.d0/(1.d0 + DEXP(-a * vh1(:)))
                case (2) !TANGENTE
                    yh1(:) = (1.d0 - DEXP(-vh1(:)))/(1.d0 + DEXP(-vh1(:)))
                case (3) !GAUSS
                    yh1(:) = DEXP(-vh1(:))
                end select

                ! ATIVACAO: CAMADA DE SAIDA
                vs = 0.d0
                vs(:) = matmul(yh1(:), ws(:,:))
                vs(:) = vs(:) - bs(:)

                select case(f_ativa)
                case (1) !LOGISTICA
                    ys(:, i) = 1.d0/(1.d0 + DEXP(-a * vs(:)))
                case (2) !TANGENTE
                    ys(:, i) = (1.d0 - DEXP(-vs(:)))/(1.d0 + DEXP(-vs(:)))
                case (3) !GAUSS
                    ys(:, i) = DEXP(-vs(:))
                end select

                !CALCULO DO ERRO RNA
                erro(:, i) = yd_valid(:, i) - ys(:, i)

                !CALCULO PADRAO DO ERRO!!!!
                erro_pad(i) = sum(erro(:, i), dim = 1)
                erro_pad(i) = 0.5d0 * (erro_pad(i)**2.d0)

            ENDDO !DO VALIDACAO

            !CALCULO DOS ERROS
            eqm_valid(1, l) = sum(erro_pad(:))
            eqm_valid(1, l) = (1.d0/(num_pad_valid)) * eqm_valid(1, l)

            if (eqm_valid(1, l) < erro_menor) then
                erro_menor = eqm_valid(1, l)
                wh1_menor = wh1
                bh1_menor = bh1
                ys_menor = ys
                ws_menor = ws
                bs_menor = bs
                eqm_menor = eqm
                eqm_valid_menor = eqm_valid
            endif

            if (cont >= 100) then
                cont = 1
                eta = eta * 0.99
            else
                cont = cont + 1
            endif

        ENDDO !DO MAXIMO ITERACAO

        !CALCULO DA FUNCAO OBJETIVO DEFINIDA (Carvalho,2011)
        penaltyObj = 100 * p1 * DEXP(DBLE(hide_neuron(1))) + p2 * DBLE(max_it) + 1
        IF (validacao .eqv. .true.) THEN
            Rede_Neural_C1 = penaltyObj * ((alphaObj * eqm(1, l) + betaObj * &
            eqm_valid(1, l)) / (alphaObj + betaObj))
        ELSE
            Rede_Neural_C1 = penaltyObj * eqm(1, l)
        ENDIF

        !GRAVA ARQUIVOS QUANDO O VALOR DA F.OBJ ATUAL FOR MELHOR QUE A ANTERIOR
        IF (Rede_Neural_C1 .LT. BestFO) THEN
            BestFO = Rede_Neural_C1
            st % bestObjectiveFunction = BestFO
            wh1_melhor = wh1_menor
            ws_melhor = ws_menor
            bh1_melhor = bh1_menor
            bs_melhor = bs_menor

            IF (op % iProcessor < 10) THEN
                WRITE (str0, '(I1)') op % iProcessor + 1
            ELSE IF (op % iProcessor < 100) THEN
                WRITE (str0, '(I2)') op % iProcessor + 1
            ELSE
                WRITE (str0, '(I3)') op % iProcessor + 1
            END IF

            IF (op % iExperiment < 10) THEN
                WRITE (str1, '(I1)') op % iExperiment
            ELSE IF (op % iExperiment < 100) THEN
                WRITE (str1, '(I2)') op % iExperiment
            ELSE
                WRITE (str1, '(I3)') op % iExperiment
            END IF

            OPEN(12, FILE = './output/nn' // trim(str1) // '_' // trim(str0) // '.out')

            WRITE(12, '(A)', ADVANCE = 'NO') 'objFun: '
            WRITE(12, '(ES14.6E2)') Rede_Neural_C1

            WRITE(12, *) 'wh1'
            DO i = 1, vetor_entrada
                WRITE(12, '(10F8.4)') (wh1_melhor(i, j), j = 1, hide_neuron(1))
            ENDDO

            WRITE(12, *) 'ws'
            DO i = 1, hide_neuron(1)
                WRITE(12, '(10F8.4)') (ws_melhor(i, j), j = 1, vetor_saida)
            ENDDO

            WRITE(12, *) 'bh1'
            WRITE(12, '(10F8.4)') (bh1_melhor(j), j = 1, hide_neuron(1))

            WRITE(12, *) 'bs'
            WRITE(12, '(10F8.4)') (bs_melhor(j), j = 1, vetor_saida)

            WRITE(12, *) 'mse'
            WRITE(12, '(15F13.4)') eqm(1, l)

            WRITE(12, *) 'mse_valid'
            WRITE(12, '(15F13.4)') (eqm_valid(1, l))

            WRITE(12, *) 'nfe'
            WRITE(12, *) st % NFE, st % totalNFE

            CLOSE(12)

        ENDIF

        Rede_Neural_C1 = BestFO

        deallocate(x, yd, x_valid, yd_valid, vh1, wh1)
        deallocate(wh1_velho, wh1_menor, bh1)
        deallocate(bh1_velho, bh1_menor, bh1_melhor, vs, ys)
        deallocate(ys_menor, ws, ws_velho)
        deallocate(ws_menor, bs, bs_velho)
        deallocate(bs_menor, bs_melhor, erro, erro_pad)
        deallocate(grad_saida, grad_hide1)
        deallocate(deltabs, deltabh1)
        deallocate(deltaw_saida, deltaw_h1)
        deallocate(eqm, eqm_valid, eqm_menor)
        deallocate(eqm_valid_menor, yh1)

        !FIM DA FUNCAO
    END FUNCTION Rede_Neural_C1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DOUBLE PRECISION FUNCTION Rede_Neural_C2(solution, op, st)

        IMPLICIT NONE

        !Parametros recebidos pela funcao (passados pelo algoritmo de otimizacao)
        double precision :: solution(:)
        character*10, str0, str1

        !Variaveis usadas para o calculo da funcao objetivo
        double precision :: BestFO
        double precision :: penaltyObj, alphaObj, betaObj, p1, p2
        double precision, allocatable, dimension(:,:) :: x
        double precision, allocatable, dimension(:,:) :: yd
        double precision, allocatable, dimension(:,:) :: x_valid
        double precision, allocatable, dimension(:,:) :: yd_valid
        double precision, allocatable, dimension(:) :: vs
        double precision, allocatable, dimension(:,:) :: ys, ys_menor
        double precision, allocatable, dimension(:,:) :: ws, ws_velho, ws_menor, ws_melhor
        double precision, allocatable, dimension(:) :: bs, bs_velho, bs_menor, bs_melhor
        double precision, allocatable, dimension(:,:) :: eqm_menor, eqm_valid_menor
        double precision, allocatable, dimension(:,:) :: erro
        double precision, allocatable, dimension(:) :: erro_pad
        double precision, allocatable, dimension(:,:) :: eqm
        double precision, allocatable, dimension(:,:) :: eqm_valid
        double precision, allocatable, dimension(:) :: grad_saida, grad_hide1, grad_hide2
        double precision, allocatable, dimension(:,:) :: deltaw_saida
        double precision, allocatable, dimension(:) :: deltabs, deltabh1, deltabh2
        double precision, allocatable, dimension(:) :: vh1, vh2
        double precision, allocatable, dimension(:,:) :: wh1, wh1_velho, wh1_menor, wh1_melhor
        double precision, allocatable, dimension(:,:) :: wh2, wh2_velho, wh2_menor, wh2_melhor
        double precision, allocatable, dimension(:) :: bh1, bh1_velho, bh1_menor, bh1_melhor
        double precision, allocatable, dimension(:) :: bh2, bh2_velho, bh2_menor, bh2_melhor
        double precision, allocatable, dimension(:) :: yh1, yh2
        double precision, allocatable, dimension(:,:) :: deltaw_h1, deltaw_h2

        double precision, parameter :: a = 1

        !PARAMETROS RECEBIDOS PELA FUNCAO (PASSADO PELO ALGORITMO DE OTIMIZACAO)
        double precision :: dv, soma, erro_menor, alpha
        double precision :: eta ! taxa de aprendizado
        double precision :: ed ! erro desejado
        double precision :: aleat ! gerar pesos e bias aleatoriamente

        integer :: hide_neuron(2) ! numero de neuronios na camada escondida
        integer :: cont
        integer :: l
        integer :: max_it ! numero maximo de iteracoes
        integer :: vetor_entrada ! dimensao do vetor de entrada
        integer :: num_pad ! numero de padroes a serem treinados
        integer :: num_pad_valid
        integer :: vetor_saida ! dimensao do vetor de saida
        integer :: i, j, k, hn ! variavel para controle do programa
        integer :: f_ativa ! controle de qual funcao ativacao
        integer :: f_deriva ! controle de qual funcao derivada
        integer :: hide_camada ! numero de camadas escondidas
        logical :: peso ! controlar como os pesos serao inicializados 
        logical :: validacao ! controlar se terah ou nao validacao

        TYPE(OptionsMPCA), intent(in) :: op
        TYPE(StatusMPCA), intent(inout) :: st

        !Atribuicao empirica dos pesos para a funcao objetivo
        alphaObj = 1.0D0
        betaObj = 1.0D0
        erro_menor = 0.1D0
        p1 = 5.0000E-08 !5*10^-8 
        p2 = 5.0000e-05 !5*10^-5

        hide_camada = nint(solution(1))
        hide_neuron(1) = nint(solution(2))
        hide_neuron(2) = nint(solution(3))
        f_ativa = nint(solution(4))
        alpha = solution(5)
        eta = solution(6)

        num_pad = op % nClasses
        num_pad_valid = op % nClassesValidation
        vetor_entrada = op % nInputs
        vetor_saida = op % nOutputs
        ed = op % targetError
        max_it = op % nEpochs
        peso = op % randomWeights
        validacao = op % haveValidation
        BestFO = st % bestObjectiveFunction
        f_deriva = f_ativa

        !Bias 1. camada escondida
        allocate(bh1(hide_neuron(1)))
        allocate(bh1_menor(hide_neuron(1)))
        allocate(bh1_velho(hide_neuron(1)))
        allocate(bh1_melhor(hide_neuron(1)))
        !Bias 2. camada escondida
        allocate(bh2(hide_neuron(2)))
        allocate(bh2_menor(hide_neuron(2)))
        allocate(bh2_velho(hide_neuron(2)))
        allocate(bh2_melhor(hide_neuron(2)))
        !Bias camada de saida
        allocate(bs(vetor_saida))
        allocate(bs_menor(vetor_saida))
        allocate(bs_velho(vetor_saida))
        allocate(bs_melhor(vetor_saida))
        !DeltaBias
        allocate(deltabh1(hide_neuron(1)))
        allocate(deltabh2(hide_neuron(2)))
        allocate(deltabs(vetor_saida))
        !DeltaW's
        allocate(deltaw_h1(vetor_entrada, hide_neuron(1)))
        allocate(deltaw_h2(hide_neuron(1), hide_neuron(2)))
        allocate(deltaw_saida(hide_neuron(hide_camada), vetor_saida))
        !Erros
        allocate(eqm(1, max_it))
        allocate(eqm_valid(1, max_it))
        allocate(eqm_valid_menor(1, max_it))
        allocate(eqm_menor(1, max_it))
        allocate(erro(vetor_saida, num_pad))
        allocate(erro_pad(num_pad))
        !Gradientes
        allocate(grad_hide1(hide_neuron(1)))
        allocate(grad_hide2(hide_neuron(2)))
        allocate(grad_saida(vetor_saida))
        !Campo local induzido
        allocate(vh1(hide_neuron(1)))
        allocate(vh2(hide_neuron(2)))
        allocate(vs(vetor_saida))
        !Pesos 1. camada escondida
        allocate(wh1(vetor_entrada, hide_neuron(1)))
        allocate(wh1_menor(vetor_entrada, hide_neuron(1)))
        allocate(wh1_melhor(vetor_entrada, hide_neuron(1)))
        allocate(wh1_velho(vetor_entrada, hide_neuron(1)))
        !Pesos 2. camada escondida
        allocate(wh2(hide_neuron(1), hide_neuron(2)))
        allocate(wh2_menor(hide_neuron(1), hide_neuron(2)))
        allocate(wh2_melhor(hide_neuron(1), hide_neuron(2)))
        allocate(wh2_velho(hide_neuron(1), hide_neuron(2)))
        !Pesos camada de saida
        allocate(ws(hide_neuron(hide_camada), vetor_saida))
        allocate(ws_menor(hide_neuron(hide_camada), vetor_saida))
        allocate(ws_velho(hide_neuron(hide_camada), vetor_saida))
        allocate(ws_melhor(hide_neuron(hide_camada), vetor_saida))
        !Dados de entrada e saida
        allocate(x(vetor_entrada, num_pad))
        allocate(x_valid(vetor_entrada, num_pad_valid))
        allocate(yd(vetor_saida, num_pad))
        allocate(yd_valid(vetor_saida, num_pad_valid))
        !Saidas obtidas
        allocate(yh1(hide_neuron(1)))
        allocate(yh2(hide_neuron(2)))
        allocate(ys(vetor_saida, num_pad))
        allocate(ys_menor(vetor_saida, num_pad))

        !------------------------------------------------------------!
        !LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
        !------------------------------------------------------------!
        OPEN (1, file = op % path_sai)
        DO I = 1, vetor_saida
            READ(1, *) (yd(I, J), J = 1, num_pad)
        END DO
        CLOSE (1)

        OPEN (2, file = op % path_ent)
        DO I = 1, vetor_entrada
            READ(2, *) (x(I, J), J = 1, num_pad)
        END DO
        CLOSE (2)

        IF (validacao .eqv. .true.) THEN
            OPEN (1, file = op % path_sai_valid)
            DO I = 1, vetor_saida
                READ(1, *) (yd_valid(I, J), J = 1, num_pad_valid)
            END DO
            CLOSE (1)

            OPEN (2, file = op % path_ent_valid)
            DO I = 1, vetor_entrada
                READ(2, *) (x_valid(I, J), J = 1, num_pad_valid)
            END DO
            CLOSE (2)
        ENDIF

        !--------------------------------------------------------------------!
        !INICIALIZANDO OS PARAMETROS: wh1, bh1, ws e bs
        !--------------------------------------------------------------------!

        !PRIMEIRA CAMADA OCULTA
        DO l = 1, vetor_entrada
            DO k = 1, hide_neuron(1)
                if (peso .eqv. .false.) then
                    wh1(l, k) = 0.5
                else
                    call random_number(aleat)
                    wh1(l, k) = aleat
                endif
            ENDDO
        ENDDO
        DO k = 1, hide_neuron(1)
            if (peso .eqv. .false.) then
                bh1(k) = 0.5
            else
                call random_number(aleat)
                bh1(k) = aleat
            endif
        ENDDO

        !SEGUNDA CAMADA OCULTA
        DO k = 1, hide_neuron(1)
            DO i = 1, hide_neuron(2)
                if (peso .eqv. .false.) then
                    wh2(k, i) = 0.5
                else
                    call random_number(aleat)
                    wh2(k, i) = aleat
                endif
            ENDDO
        ENDDO
        DO i = 1, hide_neuron(2)
            if (peso .eqv. .false.) then
                bh2(i) = 0.5
            else
                call random_number(aleat)
                bh2(i) = aleat
            endif
        ENDDO

        !CAMADA DE SAIDA
        DO l = 1, hide_neuron(hide_camada)
            DO k = 1, vetor_saida
                if (peso .eqv. .false.) then
                    ws(l, k) = 0.5
                else
                    call random_number(aleat)
                    ws(l, k) = aleat
                endif
            ENDDO
        ENDDO
        DO k = 1, vetor_saida
            if (peso .eqv. .false.) then
                bs(k) = 0.5
            else
                call random_number(aleat)
                bs(k) = aleat
            endif
        ENDDO

        !----------------------------------------------------------------------!
        ! INICIO DA REDE: FEEDFORWARD
        !----------------------------------------------------------------------!
        cont = 1 !Contador de epocas
        l = 0

        !DO WHILE ((erro_menor .GE. ed) .AND. (l .LT. max_it))
        DO WHILE (l .LT. max_it)

            l = l + 1

            DO i = 1, num_pad !DO NUMERO TOTAL DE PADROES

                ! ATIVACAO: 1. CAMADA OCULTA
                vh1 = 0.d0
                vh1(:) = matmul(x(:, i), wh1(:,:))
                vh1(:) = vh1(:) - bh1(:);

                SELECT CASE(f_ativa)
                CASE (1) !LOGISTICA
                    yh1(:) = 1.d0/(1.d0 + DEXP(-a * vh1(:)))
                CASE (2) !TANGENTE
                    yh1(:) = (1.d0 - DEXP(-vh1(:)))/(1.d0 + DEXP(-vh1(:)))
                CASE (3) !GAUSS
                    yh1(:) = DEXP(-vh1(:))
                END SELECT

                ! ATIVACAO: 2. CAMADA OCULTA
                vh2(:) = 0.0
                vh2(:) = matmul(yh1(:), wh2(:,:))
                vh2(:) = vh2(:) - bh2(:)
                SELECT CASE(f_ativa)
                CASE (1) !LOGISTICA
                    yh2(:) = 1.d0/(1.d0 + DEXP(-a * vh2(:)))
                CASE (2) !TANGENTE
                    yh2(:) = (1.d0 - DEXP(-vh2(:)))/(1.d0 + DEXP(-vh2(:)))
                CASE (3) !GAUSS
                    yh2(:) = DEXP(-vh2(:))
                END SELECT

                ! ATIVACAO: CAMADA DE SAIDA
                vs = 0.d0
                vs(:) = MATMUL(yh2(:), ws(:,:))
                vs(:) = vs(:) - bs(:)
                SELECT CASE(f_ativa)
                CASE (1) !LOGISTICA
                    ys(:, i) = 1.d0 / (1.d0 + DEXP(-a * vs(:)))
                CASE (2) !TANGENTE
                    ys(:, i) = (1.d0 - DEXP(-vs(:))) / (1.d0 + DEXP(-vs(:)))
                CASE (3) !GAUSS
                    ys(:, i) = DEXP(-vs(:))
                END SELECT

                ! PARA CALCULO DO NOVO PESO
                wh1_velho = wh1
                bh1_velho = bh1
                wh2_velho = wh2
                bh2_velho = bh2
                ws_velho = ws
                bs_velho = bs

                ! CALCULO ERRO TREINAMENTO
                erro(:, i) = yd(:, i) - ys(:, i)

                !-------------------------------------------------------------------------!
                !                        BACKPROPAGATION
                !-------------------------------------------------------------------------!
                ! TREINAMENTO: CAMADA DE SAIDA
                DO j = 1, vetor_saida

                    SELECT CASE(f_deriva)
                    CASE (1)!LOGISTICA: e^-x / (1+e^-x)^2
                        dv = ((a * DEXP(-a * vs(j))) / ((1.d0 + DEXP(-a * vs(j))) ** 2.d0))
                    CASE (2)!TANGENTE: 2e^-x / (1+e^-x)^2
                        dv = (2 * DEXP(-vs(j))) / ((1 + DEXP(-vs(j))) ** 2)
                    CASE (3)!GAUSS
                        dv = -ys(j, i) / a
                    END SELECT

                    grad_saida(j) = erro(j, i) * dv
                    deltaw_saida(:, j) = eta * grad_saida(j) * yh2(:)
                    ws(:, j) = ws(:, j) + alpha * (ws(:, j) - ws_velho(:, j)) + deltaw_saida(:, j)
                    deltabs(j) = eta * grad_saida(j) * (-1.d0)
                    bs(j) = bs(j) + deltabs(j)

                enddo !CAMADA SAIDA

                !TREINAMENTO: 2. CAMADA OCULTA
                DO j = 1, hide_neuron(2)

                    soma = 0.d0
                    do k = 1, vetor_saida
                        soma = soma + (grad_saida(k) * ws_velho(j, k))
                    enddo

                    select case(f_deriva)
                    case (1)!LOGISTICA: e^-x / (1+e^-x)^2
                        dv = ((a * DEXP(-a * vh2(j)))/((1.d0 + DEXP(-a * vh2(j)))**2.d0))
                    case (2)!TANGENTE
                        dv = (2 * DEXP(-vh2(j)))/ ((1 + DEXP(-vh2(j)))**2)
                    case (3)!GAUSS
                        dv = -yh2(j)/a
                    end select

                    grad_hide2(j) = dv * soma
                    deltaw_h2(:, j) = eta * grad_hide2(j) * yh1(:)
                    wh2(:, j) = wh2(:, j) + alpha * (wh2(:, j) - wh2_velho(:, j)) + deltaw_h2(:, j)
                    deltabh2(j) = eta * grad_hide2(j) * (-1.d0)
                    bh2(j) = bh2(j) + deltabh2(j)

                ENDDO !2. CAMADA OCULTA

                !TREINAMENTO: 1. CAMADA OCULTA
                DO j = 1, hide_neuron(1)

                    soma = 0.d0
                    do k = 1, hide_neuron(2)
                        soma = soma + (grad_hide2(k) * wh2_velho(j, k))
                    enddo

                    select case(f_deriva)
                    case (1)!LOGISTICA: e^-x / (1+e^-x)^2
                        dv = ((a * DEXP(-a * vh1(j))) / ((1.d0 + DEXP(-a * vh1(j))) ** 2.d0))
                    case (2)!TANGENTE
                        dv = (2 * DEXP(-vh1(j))) / ((1 + DEXP(-vh1(j))) ** 2)
                    case (3)!GAUSS
                        dv = -yh1(j) / a
                    end select

                    grad_hide1(j) = dv * soma
                    deltaw_h1(:, j) = eta * grad_hide1(j) * x(:, i)
                    wh1(:, j) = wh1(:, j) + alpha * (wh1(:, j) - wh1_velho(:, j)) + deltaw_h1(:, j)
                    deltabh1(j) = eta * grad_hide1(j) * (-1.d0)
                    bh1(j) = bh1(j) + deltabh1(j)

                ENDDO !1. CAMADA OCULTA

                !CALCULO PADRAO DO ERRO
                erro_pad(i) = sum(erro(:, i), dim = 1)
                erro_pad(i) = 0.5d0 * (erro_pad(i)**2.d0)

            ENDDO !NUMERO PADROES

            eqm(1, l) = SUM(erro_pad(:))
            eqm(1, l) = (1.d0 / (num_pad)) * eqm(1, l)


            !***********************************************************************
            ! VALIDACAO
            !***********************************************************************

            do i = 1, num_pad_valid

                ! ATIVACAO: 1. CAMADA OCULTA
                vh1 = 0.d0
                vh1(:) = matmul(x_valid(:, i), wh1(:,:))
                vh1(:) = vh1(:) - bh1(:);
                select case(f_ativa)
                case (1) !LOGISTICA
                    yh1(:) = 1.d0 / (1.d0 + DEXP(-a * vh1(:)))
                case (2) !TANGENTE
                    yh1(:) = (1.d0 - DEXP(-vh1(:))) / (1.d0 + DEXP(-vh1(:)))
                case (3) !GAUSS
                    yh1(:) = DEXP(-vh1(:))
                end select

                ! ATIVACAO: 2. CAMADA OCULTA
                vh2(:) = 0.0
                vh2(:) = MATMUL(yh1(:), wh2(:,:))
                vh2(:) = vh2(:) - bh2(:)

                SELECT CASE(f_ativa)
                CASE (1) !LOGISTICA
                    yh2(:) = 1.d0/(1.d0 + DEXP(-a * vh2(:)))
                CASE (2) !TANGENTE
                    yh2(:) = (1.d0 - DEXP(-vh2(:)))/(1.d0 + DEXP(-vh2(:)))
                CASE (3) !GAUSS
                    yh2(:) = DEXP(-vh2(:))
                END SELECT

                ! ATIVACAO: CAMADA DE SAIDA
                vs = 0.d0
                vs(:) = MATMUL(yh2(:), ws(:,:))
                vs(:) = vs(:) - bs(:)

                select case(f_ativa)
                CASE (1) !LOGISTICA
                    ys(:, i) = 1.d0 / (1.d0 + DEXP(-a * vs(:)))
                CASE (2) !TANGENTE
                    ys(:, i) = (1.d0 - DEXP(-vs(:))) / (1.d0 + DEXP(-vs(:)))
                CASE (3) !GAUSS
                    ys(:, i) = DEXP(-vs(:))
                END SELECT

                !CALCULO DO ERRO RNA
                erro(:, i) = yd_valid(:, i) - ys(:, i)

                !CALCULO PADRAO DO ERRO
                erro_pad(i) = sum(erro(:, i), dim = 1)
                erro_pad(i) = 0.5d0 * (erro_pad(i)**2.d0)

            ENDDO !DO VALIDACAO

            !CALCULO DOS ERROS
            eqm_valid(1, l) = sum(erro_pad(:))
            eqm_valid(1, l) = (1.d0/(num_pad_valid)) * eqm_valid(1, l)

            if (eqm_valid(1, l) < erro_menor) then
                erro_menor = eqm_valid(1, l)
                wh1_menor = wh1
                wh2_menor = wh2
                bh1_menor = bh1
                bh2_menor = bh2
                ys_menor = ys
                ws_menor = ws
                bs_menor = bs
                eqm_menor = eqm
                eqm_valid_menor = eqm_valid
            endif

            if (cont >= 100) then
                cont = 1
                eta = eta * 0.99
            else
                cont = cont + 1
            endif

        ENDDO !DO MAXIMO ITERAÇÃO DO WHILE (l .LT. max_it)

        !CALCULO DA FUNCAO OBJETIVO DEFINIDA (Carvalho,2011)
        penaltyObj = 100 * p1 * DEXP(DBLE(hide_neuron(1))) + 1e+6 * p1 * DEXP(DBLE(hide_neuron(2))) + p2 * DBLE(max_it) + 1
        IF (validacao .eqv. .true.) THEN
            Rede_Neural_C2 = penaltyObj * ((alphaObj * eqm(1, l) + betaObj * &
            eqm_valid(1, l)) / (alphaObj + betaObj))
        ELSE
            Rede_Neural_C2 = penaltyObj * eqm(1, l)
        ENDIF

        !GRAVA ARQUIVOS QUANDO O VALOR DA F.OBJ ATUAL FOR MELHOR QUE A ANTERIOR
        IF (Rede_Neural_C2 .LT. BestFO) THEN
            BestFO = Rede_Neural_C2
            st % bestObjectiveFunction = BestFO
            wh1_melhor = wh1_menor
            wh2_melhor = wh2_menor
            ws_melhor = ws_menor
            bh1_melhor = bh1_menor
            bh2_melhor = bh2_menor
            bs_melhor = bs_menor

            IF (op % iProcessor < 10) THEN
                WRITE (str0, '(I1)') op % iProcessor + 1
            ELSE IF (op % iProcessor < 100) THEN
                WRITE (str0, '(I2)') op % iProcessor + 1
            ELSE
                WRITE (str0, '(I3)') op % iProcessor + 1
            END IF

            IF (op % iExperiment < 10) THEN
                WRITE (str1, '(I1)') op % iExperiment
            ELSE IF (op % iExperiment < 100) THEN
                WRITE (str1, '(I2)') op % iExperiment
            ELSE
                WRITE (str1, '(I3)') op % iExperiment
            END IF

            OPEN(12, FILE = './output/nn' // trim(str1) // '_' // trim(str0) // '.out')

            WRITE(12, '(A)', ADVANCE = 'NO') 'objFun: '
            WRITE(12, '(ES14.6E2)') Rede_Neural_C2

            WRITE(12, *) 'wh1'
            DO i = 1, vetor_entrada
                WRITE(12, '(10F8.4)') (wh1_melhor(i, j), j = 1, hide_neuron(1))
            ENDDO

            WRITE(12, *) 'wh2'
            DO i = 1, hide_neuron(1)
                WRITE(12, '(10F8.4)') (wh2_melhor(i, j), j = 1, hide_neuron(2))
            ENDDO

            WRITE(12, *) 'ws'
            DO i = 1, hide_neuron(2)
                WRITE(12, '(10F8.4)') (ws_melhor(i, j), j = 1, vetor_saida)
            ENDDO

            WRITE(12, *) 'bh1'
            WRITE(12, '(10F8.4)') (bh1_melhor(j), j = 1, hide_neuron(1))

            WRITE(12, *) 'bh2'
            WRITE(12, '(10F8.4)') (bh2_melhor(j), j = 1, hide_neuron(2))

            WRITE(12, *) 'bs'
            WRITE(12, '(10F8.4)') (bs_melhor(j), j = 1, vetor_saida)

            WRITE(12, *) 'mse'
            WRITE(12, '(15F13.4)') eqm(1, l)

            WRITE(12, *) 'mse_valid'
            WRITE(12, '(15F13.4)') (eqm_valid(1, l))

            WRITE(12, *) 'nfe'
            WRITE(12, *) st % NFE, st % totalNFE

            CLOSE(12)

        ENDIF

        Rede_Neural_C2 = BestFO

        deallocate(x, yd, x_valid, yd_valid, vh1, wh1)
        deallocate(wh1_velho, wh1_menor, wh1_melhor, bh1)
        deallocate(bh1_velho, bh1_menor, bh1_melhor, vs, ys)
        deallocate(ys_menor, ws, ws_velho)
        deallocate(ws_menor, bs, bs_velho)
        deallocate(bs_menor, bs_melhor, erro, erro_pad)
        deallocate(grad_saida, grad_hide1)
        deallocate(deltabs, deltabh1)
        deallocate(deltaw_saida, deltaw_h1)
        deallocate(eqm, eqm_valid, eqm_menor)
        deallocate(eqm_valid_menor, yh1)
        deallocate(vh2, wh2, wh2_velho, wh2_menor, wh2_melhor, bh2)
        deallocate(bh2_velho, bh2_menor, yh2)
        deallocate(grad_hide2, deltabh2, deltaw_h2)

    END FUNCTION Rede_Neural_C2

    !Fim do modulo
END MODULE ModuloRNA
