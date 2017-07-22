program main

    use f90getopt ! Tratar parâmetros de entrada
    use commonUtils, only: checkPoisson, lower, delete_file ! funções/subrotinas
    use help, only: help_options

    implicit none

    real              :: start, finish
    character(len=4), parameter  :: tab = '    '
    character(len=50) :: nfile, Ifile
    character(len=12) :: mode
    character(len=1)  :: var
    double precision  :: I1, I2, I_med, w_exp
    double precision  :: a, b, tol, dq, w1, w2, w_med, h
    integer           :: i , k, k_max
    integer           :: metodo ! 1 - falsa posição; 2 - secante; 
                                ! 3 - gera 'k_max' pontos para plotar o gráfico
    logical           :: fexist
    integer           :: fflag = 0, Iflag = 0
    integer           :: nroot, OS

    ! Varíaveis para mostrar a hora de início
    integer :: t1, t2, count_rate
    character(len=10) :: dateandtime(3)
    integer           :: values(8)
    integer           :: ccount
    double precision  :: tmp, fixedvar, rho

    ! ======= Argumentos =======

    type(option_s) :: opts(5)

    opts(1) = option_s("version", .false., 'v')
    opts(2) = option_s("clean", .false., 'c')
    opts(3) = option_s("help", .false., 'h')
    opts(4) = option_s("include", .true., 'I')
    opts(5) = option_s("file", .true., 'f')

    call system_clock(count_rate = count_rate)

    do
        select case (getopt( "vchI:f:", opts))
            case (char(0))
                ccount = command_argument_count()
                if (ccount == 0 ) then
                    write(*,'(a)') "Nenhuma entrada encontrada"
                    stop
                endif
                exit
            case ('v')
                write(*,'(a)') 'Version 3.1415'
                stop
            case ('c')
                call delete_file("fort.3")
                call delete_file("jacobi.res")
                call delete_file("Modal.res")
                stop
            case ('h')
                call help_options()
            case ('I')
                Iflag = 1
                Ifile = optarg
            case ('f')
                fflag = 1
                nfile = optarg
            case ('?')
                stop
            case default
                call help_options()
        end select
    enddo

    if (fflag == 0 ) call get_command_argument(1, nfile)

    inquire(FIlE=trim(nfile), EXIST=fexist)

    if (.not. fexist) then
        write(*, '(a)') 'Input file does not exist'
        call exit(1)
    endif

    if ( Iflag == 1 ) then
        inquire(FIlE=trim(Ifile), EXIST=fexist)
        if (.not. fexist) then
            write(*, '(a)') 'Included file does not exist'
            call exit(1)
        endif
    endif

    ! ==================== Terminado tratamento de argumentos ===============

    call get_OS(OS)
    call date_and_time(dateandtime(1), dateandtime(2), dateandtime(3), values)
    write(*,'(a,3(i0,a))') 'Start: ', values(5),'h', values(6),'min',values(7),'s'
    print*,

    write(*,'(a)') 'INPUT'
    write(*,'(3a)') tab, 'file: ', trim(nfile)

    if (Iflag == 1 ) then
        write(*,'(3a)') tab, 'Included file: ', trim(Ifile)
    endif

    open(unit=17, file=trim(nfile), status='old')

    read(17,*) mode

    read(17,*) metodo, I1, I2, w_exp, nroot, tol, k_max
    read(17,*) tmp, tmp, tmp, fixedvar, rho

    close(17)

    select case(lower(trim(mode)))
        case('find_elastic')
            var = 'E'
            write(*,'(2a)') tab, "Calculate Young's Modulus"
        case('find_poisson')
            var = 'v'
            write(*,'(2a)') tab, "Calculate Poisson's Ratio"
    end select

    call cpu_time(start)
    
    dq = tol + 1e4 ! Só para garantir que irá entrar pelo menos no primeiro loop...
    k = 0

    select case (metodo)

        case(1) ! Metodo da falsa posição
            write(*,'(2a)') tab, "False position method"

            write(*,'(2a,f0.7)') tab, "rho   = ", rho

            select case(lower(trim(mode)))
                case('find_elastic')
                    write(*,'(2a,f0.7)') tab,   "v     = ", fixedvar
                case('find_poisson')
                    write(*,'(2a,en13.4)') tab, "E     = ", fixedvar
            end select

            write(*,'(2a,f0.7)') tab, "w_exp = ", w_exp

            open(unit=15,file='iterations_fpos.res') ! Criar arquivo para armazenar
                                                     ! os dados obtidos durante a iteração
            write(15,*)
            close(15)

            write(*,'(a,2(a,en15.4))') tab, "Initial values: ", I1, tab, I2
            write(*,'(2a,es12.3)') tab, 'Tolerance: ', tol
            write(*,'(2a,i0)') tab, 'Maximum iterations: ',k_max

            open(15, file='iterations_fpos.res', status='old', position='append', action='write')
            write(15,'(a,2(en15.4))') 'Initial values: ', I1, I2
            close(15)

            print*, 
            write(*,'(a)') 'STARTING COMPUTATIONS'
            write(*,'(2a)') tab,"CALCULATING INITIAL FREQUENCY VALUES"

            if ( Iflag == 1 ) then
                open(unit=33,file=trim(Ifile), status='old',action='read')
                read(33,*) w1, w2
                close(33)
                write(*,'(2a)') tab, 'USING INCLUDED FILE'
                write(*,'(2a,f0.7)') tab, 'w_val1 = ', w1
                write(*,'(2a,f0.7)') tab, 'w_val2 = ', w2
            else

            select case(lower(trim(mode)))
                case('find_elastic')

                    call system_clock(count=t1)
                    call firstfreqE(nfile, I1, w1)
                    call system_clock(count=t2)

                    write(*,'(2a,f0.5,a,f0.5,a)') tab, &
                        'Time per comput.: ',real(t2-t1)/real(count_rate),'s, ',&
                        real(t2-t1)/(60.*real(count_rate)),'min'

                    write(*,'(2a,f0.7)') tab, 'w_val1 = ', w1

                    call firstfreqE(nfile, I2, w2)

                    write(*,'(2a,f0.7)') tab, 'w_val2 = ', w2

                case('find_poisson')

                    call system_clock(count=t1)
                    call firstfreqPoisson(nfile, I1, w1)
                    call system_clock(count=t2)

                    write(*,'(2a,f0.5,a,f0.5,a)') tab, &
                        'Time per comput.: ',real(t2-t1)/real(count_rate),'s, ',&
                        real(t2-t1)/(60.*real(count_rate)),'min'

                    write(*,'(2a,f0.7)') tab, 'w_val1 = ', w1

                    call firstfreqPoisson(nfile, I2, w2)

                    write(*,'(2a,f0.7)') tab, 'w_val2 = ', w2

            end select
            endif

            print*,
            write(*,'(2a)') tab,'ITERATIONS'

            do while(dq > tol .and. k<k_max)
                k = k + 1 !  número de iteracoes

                write(*,'(a,i0,a)') tab, k, ': '

                open(15, file='iterations_fpos.res', status='old', position='append', action='write')
                write(15,*) 'k =', k

                write(15,'(2a,en15.4)') var,'1 = ', I1
                write(15,'(2a,en15.4)') var,'2 = ', I2

                a = w_exp - w1
                write(15,*) 'w1 = ', w1
                write(15,*) 'a = ',a

                b = w_exp - w2
                write(15,*) 'w2 = ', w2
                write(15,*) 'b = ',b

                I_med = (I1*b - I2*a)/(b - a)
                write(15,'(2a,en17.8)') var,'_med = ', I_med
                close(15)

                select case(lower(trim(mode)))
                    case('find_elastic')
                        if (k==1 .and. Iflag == 1) then
                            call system_clock(count=t1)
                            call firstfreqE(nfile, I_med, w_med)
                            call system_clock(count=t2)

                            write(*,'(2a,f0.5,a,f0.5,a)') tab, &
                                'Time per comput.: ',real(t2-t1)/real(count_rate),'s, ',&
                                real(t2-t1)/(60.*real(count_rate)),'min'
                        else
                            call firstfreqE(nfile, I_med, w_med)
                        endif
                    case('find_poisson')

                        if (.not. checkPoisson(I_med)) then
                            write(*,'(2a,f0.7)') tab, 'STOP: v_med = ', I_med
                            stop
                        endif

                        if (k==1 .and. Iflag == 1) then
                            call system_clock(count=t1)
                            call firstfreqPoisson(nfile, I_med, w_med)
                            call system_clock(count=t2)

                            write(*,'(2a,f0.5,a,f0.5,a)') tab, &
                                'Time per comput.: ',real(t2-t1)/real(count_rate),'s, ',&
                                real(t2-t1)/(60.*real(count_rate)),'min'
                        else
                            call firstfreqPoisson(nfile, I_med, w_med)
                        endif
                end select

                open(15, file='iterations_fpos.res', status='old', position='append', action='write')
                write(15,*) 'w_med = ', w_med

                dq = dabs(w_exp - w_med)
                write(15,*) 'dq = ', dq
                write(15,*)
                close(15)

                write(*,fmt="(a1)", advance="no") achar(13)

                select case(lower(trim(mode)))
                    case('find_elastic')
                        write(*,'(3a,f0.7,2a,en16.5,a)') &
                            tab, tab, 'w_med = ',w_med, tab, '(E_med = ', I_med, ')'
                    case('find_poisson')
                        write(*,'(3a,f0.7,2a,f0.7,a)') &
                            tab, tab, 'w_med = ',w_med, tab, '(v_med = ', I_med, ')'
                end select

                if (a*(w_exp - w_med) > 0) then ! Teste para ver onde a função troca de sinal
                    I1 = I_med
                    w1 = w_med
                else
                    I2 = I_med
                    w2 = w_med
                endif

            enddo

            write(*,fmt="(a1,t20,a)", advance="no") achar(13), ' '
            print*,

            write(*,'(a)') 'OUTPUT'

            select case(lower(trim(mode)))
                case('find_elastic')
                    write(*,'(3a,en17.8)') tab, var, ' = ', I_med
                case('find_poisson')
                    write(*,'(3a,f0.8)') tab, var, ' = ', I_med
            end select

            write(*,'(2a,f0.7)') tab, 'w_med = ', w_med
            write(*,'(2a,es10.4)') tab,'Relative error: ', dabs(w_exp - w_med)/w_exp
            write(*,'(2a,i0)') tab, 'Iterations: ', k

        case(2) ! método da secante

            write(*,'(2a)') tab, "Secant method"

            write(*,'(2a,f0.7)') tab, "rho   = ", rho

            select case(lower(trim(mode)))
                case('find_elastic')
                    write(*,'(2a,f0.7)') tab,   "v     = ", fixedvar
                case('find_poisson')
                    write(*,'(2a,en13.4)') tab, "E     = ", fixedvar
            end select

            write(*,'(2a,f0.7)') tab, "w_exp = ", w_exp

            open(unit=15,file='iterations_sec.res')
            write(15,*)
            close(15)

            write(*,'(a,2(a,en15.4))') tab, "Initial values: ", I1, tab, I2
            write(*,'(2a,es12.3)') tab, 'Tolerance: ', tol
            write(*,'(2a,i0)') tab, 'Maximum iterations: ',k_max

            open(15, file='iterations_sec.res', status='old', position='append', action='write')
            write(15,'(a,2(en15.4))') 'Initial values: ', I1, I2
            close(15)

            print*, 
            write(*,'(a)') 'STARTING COMPUTATIONS'
            write(*,'(2a)') tab,"CALCULATING INITIAL FREQUENCY VALUES"

            if ( Iflag == 1 ) then
                open(unit=33,file=trim(Ifile), status='old',action='read')
                read(33,*) w1, w2
                close(33)
                write(*,'(2a)') tab, 'USING INCLUDED FILE'
                write(*,'(2a,f0.7)') tab, 'w_val1 = ', w1
                write(*,'(2a,f0.7)') tab, 'w_val2 = ', w2
            else

            select case(lower(trim(mode)))
                case('find_elastic')

                    call system_clock(count=t1)
                    call firstfreqE(nfile, I1, w1)
                    call system_clock(count=t2)

                    write(*,'(2a,f0.5,a,f0.5,a)') tab, &
                        'Time per comput.: ',real(t2-t1)/real(count_rate),'s, ',&
                        real(t2-t1)/(60.*real(count_rate)),'min'

                    write(*,'(2a,f0.7)') tab, 'w_val1 = ', w1

                    call firstfreqE(nfile, I2, w2)

                    write(*,'(2a,f0.7)') tab, 'w_val2 = ', w2

                case('find_poisson')

                    call system_clock(count=t1)
                    call firstfreqPoisson(nfile, I1, w1)
                    call system_clock(count=t2)

                    write(*,'(2a,f0.5,a,f0.5,a)') tab, &
                        'Time per comput.: ',real(t2-t1)/real(count_rate),'s, ',&
                        real(t2-t1)/(60.*real(count_rate)),'min'

                    write(*,'(2a,f0.7)') tab, 'w_val1 = ', w1

                    call firstfreqPoisson(nfile, I2, w2)

                    write(*,'(2a,f0.7)') tab, 'w_val2 = ', w2

            end select
            endif

            print*,
            write(*,'(2a)') tab,'ITERATIONS'

            do while (k < k_max .and. dabs(dq)>tol)
                k = k + 1

                write(*,'(a,i0,a)') tab, k, ': '

                open(15, file='iterations_sec.res', status='old', position='append', action='write')
                write(15,*) 'k =', k

                write(15,'(2a,en15.4)') var,'1 = ', I1
                write(15,'(2a,en15.4)') var,'2 = ', I2

                a = w_exp - w1
                write(15,*) 'w1 = ', w1
                write(15,*) 'a = ',a

                b = w_exp - w2
                write(15,*) 'w2 = ', w2
                write(15,*) 'b = ',b

                I_med = I2 - (I2 - I1)*b/(b - a)
                write(15,'(2a,en17.8)') var,'_med = ', I_med
                close(15)

                select case(lower(trim(mode)))
                    case('find_elastic')
                        if (k==1 .and. Iflag == 1) then
                            call system_clock(count=t1)
                            call firstfreqE(nfile, I_med, w_med)
                            call system_clock(count=t2)

                            write(*,'(2a,f0.5,a,f0.5,a)') tab, &
                                'Time per comput.: ',real(t2-t1)/real(count_rate),'s, ',&
                                real(t2-t1)/(60.*real(count_rate)),'min'
                        else
                            call firstfreqE(nfile, I_med, w_med)
                        endif
                    case('find_poisson')

                        if (.not. checkPoisson(I_med)) then
                            write(*,'(2a,f0.7)') tab, 'STOP: v_med = ', I_med
                            stop
                        endif

                        if (k==1 .and. Iflag == 1) then
                            call system_clock(count=t1)
                            call firstfreqPoisson(nfile, I_med, w_med)
                            call system_clock(count=t2)

                            write(*,'(2a,f0.5,a,f0.5,a)') tab, &
                                'Time per comput.: ',real(t2-t1)/real(count_rate),'s, ',&
                                real(t2-t1)/(60.*real(count_rate)),'min'
                        else
                            call firstfreqPoisson(nfile, I_med, w_med)
                        endif
                end select

                open(15, file='iterations_sec.res', status='old', position='append', action='write')
                write(15,*) 'w_med = ', w_med

                dq = dabs(w_exp - w_med)
                write(15,*) 'dq = ', dq

                write(*,fmt="(a1)", advance="no") achar(13)

                select case(lower(trim(mode)))
                    case('find_elastic')
                        write(*,'(3a,f0.7,2a,en13.5,a)') &
                            tab, tab, 'w_med = ',w_med, tab, '(E_med = ', I_med, ')'
                    case('find_poisson')
                        write(*,'(3a,f0.7,2a,f0.7,a)') &
                            tab, tab, 'w_med = ',w_med, tab, '(v_med = ', I_med, ')'
                end select

                I1 = I2
                w1 = w2
                I2 = I_med
                w2 = w_med

                write(15,*)
                close(15)

            end do

            write(*,fmt="(a1,t20,a)", advance="no") achar(13), ' '
            print*,

            write(*,'(a)') 'OUTPUT'

            select case(lower(trim(mode)))
                case('find_elastic')
                    write(*,'(3a,en17.8)') tab, var, ' = ', I_med
                case('find_poisson')
                    write(*,'(3a,f0.8)') tab, var, ' = ', I_med
            end select

            write(*,'(2a,f0.7)') tab, 'w_med = ', w_med
            write(*,'(2a,es10.4)') tab,'Relative error: ', dabs(w_exp - w_med)/w_exp
            write(*,'(2a,i0)') tab, 'Iterations: ', k


        case(3) ! Pontos pora plotar função desvio-elasticidade
            ! Precisa de revisão

            h = (I2-I1)/(k_max-1)

            open(unit=15, file = 'points.dat')
            write(15,*)
            close(15)
            write(*,'(a)') "Generating points in 'points.dat' file"

            do i=0, k_max - 1
                I_med = I1 + h*i
                select case(lower(trim(mode)))
                    case('find_elastic')
                        call firstfreqE(nfile, I_med, w_med)
                    case('find_poisson')
                        call firstfreqPoisson(nfile, I_med, w_med)
                end select
                open(15, file='points.dat', status='old', position='append', action='write')
                write(15,*) I_med, w_med, w_exp - w_med
                close(15)
            enddo

    end select

    call cpu_time(finish)
    call date_and_time(dateandtime(1), dateandtime(2), dateandtime(3), values)

    print*,
    write(*,'(a,3(i0,a))') 'End: ', values(5),'h', values(6),'min',values(7),'s'

    call delete_file("jacobi.res")
    call delete_file("Modal.res")
    ! TODO: arrumar esse 'rm fort.3'
    close(3)
    call delete_file("fort.3")
    !call execute_command_line("rm fort.3")

end program main
