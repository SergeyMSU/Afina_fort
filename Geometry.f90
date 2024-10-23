module GEOMETRY
    use STORAGE 
    implicit none 


    contains

    subroutine Init_Setka(SS, N)        ! Выделение памяти под все массивы сетки
        TYPE (Setka), intent(in out) :: SS
        INTEGER, INTENT(IN) :: N   ! Сколько граней будет в сетке

        SS%N = N

        allocate(SS%gran(2, N))
        allocate(SS%if_u_0(N))
        allocate(SS%if_un_0(N))
        allocate(SS%gran_info(N))

        allocate(SS%center_gran_x(N))
        allocate(SS%center_gran_y(N))
        allocate(SS%normal_gran_x(N))
        allocate(SS%normal_gran_y(N))
        allocate(SS%u(N))
        allocate(SS%if_u(N))
        allocate(SS%un(N))
        allocate(SS%if_un(N))
        allocate(SS%len_gran(N))
        allocate(SS%MM(N, N))
        allocate(SS%BB(N))

    end subroutine Init_Setka


    subroutine Dell_Setka(SS)        ! Удаляем сетку - очищаем память
        TYPE (Setka), intent(in out) :: SS

        deallocate(SS%gran)
        deallocate(SS%if_u_0)
        deallocate(SS%if_un_0)

        deallocate(SS%center_gran_x)
        deallocate(SS%center_gran_y)
        deallocate(SS%normal_gran_x)
        deallocate(SS%normal_gran_y)
        deallocate(SS%gran_info)
        deallocate(SS%u)
        deallocate(SS%if_u)
        deallocate(SS%un)
        deallocate(SS%if_un)
        deallocate(SS%len_gran)
        deallocate(SS%MM)
        deallocate(SS%BB)
    end subroutine Dell_Setka


    subroutine ALL_konstruct(N)   ! Заполняем всю сетку и готовим её к вычислениям
        INTEGER, INTENT(IN) :: N   ! Сколько точек будет на сфере (на каждой)
        INTEGER :: i, j
        real(8) :: the, r


        ALLOCATE(yzel_x(2 * N))
        ALLOCATE(yzel_y(2 * N))
        ALLOCATE(vel_gran(N))
        ALLOCATE(vel_yzel_x(N))
        ALLOCATE(vel_yzel_y(N))


        do i = 1, N
            the = (i - 1) * 2.0 * par_pi/N
            r = par_R0

            yzel_x(i) = r * cos(the)
            yzel_y(i) = r * sin(the)
        end do


        do i = N + 1, 2 * N
            the = (i - N - 1) * 2.0 * par_pi/N
            r = par_R1

            yzel_x(i) = r * cos(the)
            yzel_y(i) = r * sin(the)
        end do

        ! Работаем с первой подсеткой

        if (ALLOCATED(gl_S_in%gran) == .False.) then
            print*, "ERROR  GEOMETRY 93u8y4r8t38rg28uh2pr8hf92uhp82ygfpgg"
            STOP
        end if

        do i = 1, N  ! Заполняем грани на внутренней сфере
            j = i + 1
            if(j > N) j = 1
            gl_S_in%gran(1, i) = j
            gl_S_in%gran(2, i) = i
            gl_S_in%if_u_0(i) = .True.
            gl_S_in%if_un_0(i) = .False.
            gl_S_in%gran_info(i) = 1
        end do

        do i = N + 1, 2 * N  ! Заполняем грани на внешней сфере
            j = i + 1
            if(j > 2 * N) j = N + 1
            gl_S_in%gran(1, i) = i
            gl_S_in%gran(2, i) = j
            gl_S_in%if_u_0(i) = .False.
            gl_S_in%if_un_0(i) = .True.
            gl_S_in%gran_info(i) = 2
        end do

        gl_S_in%if_u = gl_S_in%if_u_0
        gl_S_in%if_un = gl_S_in%if_un_0

        call Calc_grans(gl_S_in)

        ! Работаем со второй подсеткой

        if (ALLOCATED(gl_S_out%gran) == .False.) then
            print*, "ERROR  GEOMETRY 93u8y4r8t38rg28uh2pr8hf92uhp82ygfpgg"
            STOP
        end if

        do i = N + 1, 2 * N  ! Заполняем грани на внешней сфере
            j = i + 1
            if(j > 2 * N) j = N + 1
            gl_S_out%gran(1, i - N) = j
            gl_S_out%gran(2, i - N) = i
            gl_S_out%if_u_0(i - N) = .False.
            gl_S_out%if_un_0(i - N) = .True.
            gl_S_out%gran_info(i - N) = 2
        end do

        gl_S_out%if_u = gl_S_out%if_u_0
        gl_S_out%if_un = gl_S_out%if_un_0

        call Calc_grans(gl_S_out)

    end subroutine ALL_konstruct


    subroutine Calc_grans(SS)   ! Считаем нормали и длины
        TYPE (Setka), intent(in out) :: SS  
        real(8) :: x1, y1, x2, y2, l1, l2, l, the, r
        INTEGER :: i, k1, k2

        do i = 1, SS%N
            k1 = SS%gran(1, i)
            k2 = SS%gran(2, i)

            x1 = yzel_x(k1)
            y1 = yzel_y(k1)

            x2 = yzel_x(k2)
            y2 = yzel_y(k2)

            l1 = x2 - x1
            l2 = y2 - y1

            l = sqrt(l1 * l1 + l2 * l2)
            SS%len_gran(i) = l

            ! Задаём нормаль, не забываем нормировать, чтобы была единичной
            SS%normal_gran_x(i) = l2/l
            SS%normal_gran_y(i) = -l1/l

            ! Центр грани
            SS%center_gran_x(i) = (x1 + x2)/2.0
            SS%center_gran_y(i) = (y1 + y2)/2.0

            r = sqrt(((x1 + x2)/2.0)**2 + ((y1 + y2)/2.0)**2)
            the = atan((y1 + y2)/2.0, (x1 + x2)/2.0)

            if(SS%gran_info(i) == 1) then
                SS%un(i) = 0.0
                SS%u(i) = par_B0 * sin(the) !! Задаём граничные условия
                ! SS%u(i) = 3551.22 * sin(2.0 * the)  !! Задаём граничные условия
                SS%if_u(i) = .True.
                SS%if_un(i) = .False.
            else
                if(SS%area == 1) then
                    SS%un(i) = 0.0
                    SS%u(i) = 0.0
                    SS%if_u(i) = .False.
                    SS%if_un(i) = .True.
                else
                    SS%un(i) = -par_Uinf * SS%normal_gran_x(i)
                    SS%u(i) = 0.0
                    SS%if_u(i) = .False.
                    SS%if_un(i) = .True.
                end if
            end if
            

        end do


    end subroutine Calc_grans


    real(8) pure function H_func(i, j, SS)
        TYPE (Setka), intent(in) :: SS  
        integer, intent(in) :: i, j
        real(8) :: xi, yi, x1, y1, x2, y2, x3, y3, r, n1, n2, r1, r2, rr, l
        integer :: k1, k2

        if(i /= j) then
            xi = SS%center_gran_x(i)
            yi = SS%center_gran_y(i)
            k1 = SS%gran(1, j)
            k2 = SS%gran(2, j)
            x1 = yzel_x(k1)
            y1 = yzel_y(k1)
            x2 = yzel_x(k2)
            y2 = yzel_y(k2)
            x3 = SS%center_gran_x(j)
            y3 = SS%center_gran_y(j)
            r = sqrt(x3**2 + y3**2)
            n1 = SS%normal_gran_x(j)
            n2 = SS%normal_gran_y(j)
            r1 = x3 - xi
            r2 = y3 - yi
            rr = r1**2 + r2**2
            l = SS%len_gran(j)
            H_func = (r1 * n1 + r2 * n2)/(rr * 2.0 * par_pi) * l
        else
            H_func = -0.5
        end if
    end function H_func


    real(8) pure function G_func(i, j, SS)
        TYPE (Setka), intent(in) :: SS  
        integer, intent(in) :: i, j
        real(8) :: xi, yi, x1, x2, y1, y2, l, r1, r2
        integer :: k1, k2

        xi = SS%center_gran_x(i)
        yi = SS%center_gran_y(i)
        k1 = SS%gran(1, j)
        k2 = SS%gran(2, j)
        x1 = yzel_x(k1)
        y1 = yzel_y(k1)
        x2 = yzel_x(k2)
        y2 = yzel_y(k2)
        l = SS%len_gran(j)
        r1 = sqrt((xi - x1)**2 + (yi - y1)**2)
        r2 = sqrt((x2 - xi)**2 + (y2 - yi)**2)

        if (i == j) then
            G_func = l * (log(l/2.0) - 1.0)/(2.0 * par_pi)
        else
            G_func = (l/(4.0 * par_pi)) * (log(r1) + log(r2))
        end if

    end function G_func


    subroutine Set_Matrix(SS)   ! Заполняем матрицу
        TYPE (Setka), intent(in out) :: SS  
        integer :: i, j
        real(8) :: S

        SS%MM = 0.0_8

        ! print*, "++", SS%if_u
        ! print*, "++", SS%if_un
        ! print*, "==", SS%gran_info

        ! Считаем MM
        do i = 1, SS%N
            do j = 1, SS%N
                if(SS%if_u(j) == .False.) then
                    SS%MM(i, j) = H_func(i, j, SS)
                else
                    SS%MM(i, j) = -G_func(i, j, SS)
                end if
            end do
        end do

        ! Считаем BB
        do i = 1, SS%N
            S = 0.0
            do j = 1, SS%N
                if( SS%if_u(j) /= .False.) then
                    S = S - H_func(i, j, SS) * SS%u(j)
                else
                    S = S + G_func(i, j, SS) * SS%un(j)
                end if
            end do
            SS%BB(i) = S
        end do

    end subroutine Set_Matrix


    subroutine Culc_equ(SS)
        TYPE (Setka), intent(in out) :: SS  
        external :: dgesv
        integer :: rc, i
        real(8) :: B(SS%N)
        real(8) :: pivot(SS%N)

        B = SS%BB


        call dgesv(SS%N, 1, SS%MM, SS%N, pivot, B, SS%N, rc)
        ! call dgesv(2, 1, M, 2, pivot, B, 2, rc)

        ! print*, "rc = ", rc
        ! print*, "B = ", B

        do i = 1, SS%N
            if (SS%if_u(i) == .False.) then
                SS%u(i) = B(i)
            else 
                SS%un(i) = B(i)
            end if
        end do

    end subroutine Culc_equ


    real(8) pure function my_sign(x)
        real(8), intent(in) :: x
        if(x > 0.0000000000001) then
            my_sign = 1.0_8
        else if (x < 0.0000000000001) then
            my_sign = -1.0_8
        else
            my_sign = 0.0_8
        end if
    end function my_sign


    subroutine Move_surface() ! Двигаем поверхность-магнитопаузу согласно давлению с двух сторон
        integer :: i, j, N, jj, o1, o2, o3, o4
        real(8) :: n1, n2, x0, y0, x1, y1, r, nn, u1, u2, ex, ey
        real(8) :: p1, p2, ux, uy, f1, f2, x2, y2, d, the, pp1, pp2, x3, y3

        ! print*, "Move_surface"
        vel_gran = 0.0_8;
        N = par_N

        open(1, file = 'test1.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'the', 'p1', 'p2'"

        open(2, file = 'test2.txt')
        write(2,*) "TITLE = 'HP'  VARIABLES = 'the', 'p1', 'p2', 'ux', 'uy', 'u1', 'u2'"

        ! Расчитываем скорость грани
        do i = N + 1, 2 * N
            j = i + 1
            if(j > 2 * N) j = N + 1
            jj = i - 1
            if(jj == N) jj = 2 * N
            x0 = gl_S_in%center_gran_x(i)
            y0 = gl_S_in%center_gran_y(i)

            ! print*, "0 ", sqrt(x0 * x0 + y0 * y0)

            n1 = gl_S_in%normal_gran_x(i)
            n2 = gl_S_in%normal_gran_y(i)

            x1 = x0 - par_otstup * n1
            y1 = y0 - par_otstup * n2
            ux = Get_v1(x1, y1, gl_S_in)
            uy = Get_v2(x1, y1, gl_S_in)
            pp1 = (ux**2 + uy**2)/(8.0 * par_pi)

            ! Новый вариант вычисления давления на границе
            p1 = 0.0_8
            do o1 = -3, 3, 1
                o2 = i + o1
                o3 = o2 + 1
                o4 = o2 - 1
                o2 = N + 1 + mod(o2 - N - 1, N)
                o3 = N + 1 + mod(o3 - N - 1, N)
                o4 = N + 1 + mod(o4 - N - 1, N)

                x1 = gl_S_in%center_gran_x(j)
                y1 = gl_S_in%center_gran_y(j)
                f1 = gl_S_in%u(j)

                x2 = gl_S_in%center_gran_x(jj)
                y2 = gl_S_in%center_gran_y(jj)
                f2 = gl_S_in%u(jj)

                ux = (f2 - f1)/sqrt((x2 - x1)**2 + (y2 - y1)**2)
                p1 = (ux**2)/(8.0 * par_pi)
            end do

            the = atan(y0, x0)
            write(1,*) the, pp1, p1

            x1 = x0 + par_otstup * n1
            y1 = y0 + par_otstup * n2
            ux = Get_v1(x1, y1, gl_S_out)
            uy = Get_v2(x1, y1, gl_S_out)
            pp2 = par_Bernully - ((ux + par_Uinf)**2 + uy**2)/2.0

            ! Новый вариант вычисления давления на границе
            x1 = gl_S_out%center_gran_x(j - N)
            y1 = gl_S_out%center_gran_y(j - N)
            f1 = gl_S_out%u(j - N)

            x2 = gl_S_out%center_gran_x(jj - N)
            y2 = gl_S_out%center_gran_y(jj - N)
            f2 = gl_S_out%u(jj - N)

            d = sqrt((x2 - x1)**2 + (y2 - y1)**2)
            u1 = (f2 - f1)/d
            u2 = -par_Uinf * (-n1)
            ex = (x2 - x1)/d
            ey = (y2 - y1)/d

            ux = u1 * ex + u2 * (-n1)
            uy = u1 * ey + u2 * (-n2)

            p2 = par_Bernully - ((ux + par_Uinf)**2 + uy**2)/2.0

            write(2,*) the, pp2, p2, ux + par_Uinf, uy, u1, u2
            
            vel_gran(i - N) = vel_gran(i - N) + par_move * sqrt(abs(p1 - p2))/2.0 * my_sign(p1 - p2)

            ! if(i == N + N/4) then
            !     print*, "i =", i, "vel_gran = ", vel_gran(i - N)
            ! end if

            ! if(i == N + 3 * N/4) then
            !     print*, "i =", i, "vel_gran = ", vel_gran(i - N)
            ! end if


            ! if(vel_gran(i - N) > 10.0) print*, "ERROR orijhfiuffr", vel_gran(i - N)
            ! print*, (i - N), vel_gran(i - N), p1, p2
        end do

        close(1)
        close(2)

        ! Расчитываем скорость каждого узла
        vel_yzel_x = 0.0_8
        vel_yzel_y = 0.0_8

        do i = N + 1, 2 * N   ! Номера узлов
            j = i - 1
            if(j == N) j = 2 * N

            n1 = gl_S_in%normal_gran_x(i)
            n2 = gl_S_in%normal_gran_y(i)

            vel_yzel_x(i - N) = vel_yzel_x(i - N) + n1 * vel_gran(i - N)/2.0
            vel_yzel_y(i - N) = vel_yzel_y(i - N) + n2 * vel_gran(i - N)/2.0

            vel_yzel_x(i - N) = vel_yzel_x(i - N) + n1 * vel_gran(j - N)/2.0
            vel_yzel_y(i - N) = vel_yzel_y(i - N) + n2 * vel_gran(j - N)/2.0

            ! if( abs(vel_yzel_y(i - N)) > 10.0) then
            !     print*, "Vel", i - N, vel_yzel_x(i - N), vel_yzel_y(i - N), " ___ ", vel_gran(i - N), vel_gran(j - N)
            !     print*, n1, n2
            ! end if
        end do

        ! Расчитываем новое положение узла
        do i = N + 1, 2 * N   ! Номера узлов  !! ДОДЕЛАТЬ
            
            x0 = yzel_x(i)
            y0 = yzel_y(i)
            r = sqrt(x0 * x0 + y0 * y0)
            nn = vel_yzel_x(i - N) * x0/r + vel_yzel_y(i - N) * y0/r
            yzel_x(i) = x0 + x0/r * nn
            yzel_y(i) = y0 + y0/r * nn
        end do


        ! Сглаживание
        do i = N + 1, 2 * N   ! Номера узлов  !! ДОДЕЛАТЬ
            j = i + 1
            if(j > 2 * N) j = N + 1
            jj = i - 1
            if(jj == N) jj = 2 * N

            x0 = yzel_x(i)
            y0 = yzel_y(i)
            x1 = yzel_x(jj)
            y1 = yzel_y(jj)
            x2 = yzel_x(j)
            y2 = yzel_y(j)

            x3 = (x1 + x2)/2.0
            y3 = (y1 + y2)/2.0

            yzel_x(i) = x0 + (x3 - x0) * par_sglag
            yzel_y(i) = y0 + (y3 - y0) * par_sglag
        end do

    end subroutine Move_surface


    real(8) pure function Get_f(x, y, SS)
        TYPE (Setka), intent(in) :: SS  
        real(8), intent(in) :: x, y
        real(8) :: S, xi, yi, x1, x2, x3, y1, y2, y3, r, n1, n2, rr1, rr2, rr, l
        INTEGER :: k1, k2, i

        S = 0.0_8

        do i = 1, SS%N
            xi = x
            yi = y
            k1 = SS%gran(1, i)
            k2 = SS%gran(2, i)
            x1 = yzel_x(k1)
            y1 = yzel_y(k1)
            x2 = yzel_x(k2)
            y2 = yzel_y(k2)
            
            x3 = SS%center_gran_x(i)
            y3 = SS%center_gran_y(i)
            r = sqrt(x3**2 + y3**2)
            n1 = SS%normal_gran_x(i)
            n2 = SS%normal_gran_y(i)
            
            rr1 = x3 - xi
            rr2 = y3 - yi
            rr = rr1**2 + rr2**2
            l = SS%len_gran(i)
            S = S + (rr1 * n1 + rr2 * n2)/(rr * 2.0 * par_pi) * l * SS%u(i)
            
            rr1 = sqrt((xi - x1)**2 + (yi - y1)**2)
            rr2 = sqrt((x2 - xi)**2 + (y2 - yi)**2)
            S = S - (l/(4.0 * par_pi)) * (log(rr1) + log(rr2)) * SS%un(i)
        end do
        Get_f = S
    end function Get_f


    real(8) pure function Get_v1(x, y, SS)
        TYPE (Setka), intent(in) :: SS  
        real(8), intent(in) :: x, y
        real(8) :: S, xi, yi, x1, x2, x3, y1, y2, y3, r, n1, n2, r1, r2, rr, l, fn1
        INTEGER :: k1, k2, i, j
        real(8) :: ksi(5), wk(5)

        ksi(1) = -0.906179845938664_8
        ksi(2) = -0.538469310105683_8
        ksi(3) = 0.0_8
        ksi(4) = 0.538469310105683_8
        ksi(5) = 0.906179845938664_8

        wk(1) = 0.236926885056189
        wk(2) = 0.478628670499366
        wk(3) = 0.568888888888889
        wk(4) = 0.478628670499366
        wk(5) = 0.236926885056189

        ! ksi(1) = -0.774596669241483_8
        ! ksi(2) = 0.0_8
        ! ksi(3) = 0.774596669241483_8

        ! wk(1) = 0.555555555555556
        ! wk(2) = 0.888888888888889
        ! wk(3) = 0.555555555555556

        S = 0.0_8

        do i = 1, SS%N
            xi = x
            yi = y
            k1 = SS%gran(1, i)
            k2 = SS%gran(2, i)
            x1 = yzel_x(k1)
            y1 = yzel_y(k1)
            x2 = yzel_x(k2)
            y2 = yzel_y(k2)
            
            n1 = SS%normal_gran_x(i)
            n2 = SS%normal_gran_y(i)
            l = SS%len_gran(i)

            do j = 1, 5
                x3 = (x2 + x1)/2.0 + ksi(j) * (x2 - x1)/2.0
                y3 = (y2 + y1)/2.0 + ksi(j) * (y2 - y1)/2.0
                ! x3 = SS%center_gran_x(i)
                ! y3 = SS%center_gran_y(i)
                r = sqrt(x3**2 + y3**2)
                
                r1 = x3 - xi
                r2 = y3 - yi
                rr = r1**2 + r2**2
                
                S = S + r1 * l * SS%un(i)/2.0/par_pi/rr * wk(j)/2.0

                fn1 = n1 * (x3 + y3 - xi - yi) * (x3 - y3 - xi + yi)/(2.0 * par_pi * rr * rr) +  &
                n2 * r1 * r2/(par_pi * rr * rr) 

                S = S + fn1 * l * SS%u(i) * wk(j)/2.0
            end do
        end do
        Get_v1 = S
    end function Get_v1


    real(8) pure function Get_v2(x, y, SS)
        TYPE (Setka), intent(in) :: SS  
        real(8), intent(in) :: x, y
        real(8) :: S, xi, yi, x1, x2, x3, y1, y2, y3, r, n1, n2, r1, r2, rr, l, fn2
        INTEGER :: k1, k2, i, j
        real(8) :: ksi(5), wk(5)

        ksi(1) = -0.906179845938664_8
        ksi(2) = -0.538469310105683_8
        ksi(3) = 0.0_8
        ksi(4) = 0.538469310105683_8
        ksi(5) = 0.906179845938664_8

        wk(1) = 0.236926885056189
        wk(2) = 0.478628670499366
        wk(3) = 0.568888888888889
        wk(4) = 0.478628670499366
        wk(5) = 0.236926885056189

        ! ksi(1) = -0.774596669241483_8
        ! ksi(2) = 0.0_8
        ! ksi(3) = 0.774596669241483_8

        ! wk(1) = 0.555555555555556
        ! wk(2) = 0.888888888888889
        ! wk(3) = 0.555555555555556

        S = 0.0_8

        do i = 1, SS%N
            xi = x
            yi = y
            k1 = SS%gran(1, i)
            k2 = SS%gran(2, i)
            x1 = yzel_x(k1)
            y1 = yzel_y(k1)
            x2 = yzel_x(k2)
            y2 = yzel_y(k2)
            
            do j = 1, 5
                x3 = (x2 + x1)/2.0 + ksi(j) * (x2 - x1)/2.0
                y3 = (y2 + y1)/2.0 + ksi(j) * (y2 - y1)/2.0
                ! x3 = SS%center_gran_x(i)
                ! y3 = SS%center_gran_y(i)
                r = sqrt(x3**2 + y3**2)
                n1 = SS%normal_gran_x(i)
                n2 = SS%normal_gran_y(i)
                
                r1 = x3 - xi
                r2 = y3 - yi
                rr = r1**2 + r2**2
                l = SS%len_gran(i)
                S = S + r2 * l * SS%un(i)/2.0/par_pi/rr * wk(j)/2.0

                fn2 = -n2 * (x3 + y3 - xi - yi) * (x3 - y3 - xi + yi)/(2.0 * par_pi * rr * rr) +  &
                n1 * r1 * r2/(par_pi * rr * rr)
                
                S = S + fn2 * l * SS%u(i) * wk(j)/2.0
            end do
        end do
        Get_v2 = S
    end function Get_v2

    subroutine Print_test()
        integer :: i
        real(8) :: x, y, f1, f2, d, r, the
        open(1, file = 'test.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'the', 'Fanalit', 'F'"

        do i = par_N + 1, 2 * par_N
            x = gl_S_in%center_gran_x(i)
            y = gl_S_in%center_gran_y(i)
            r = sqrt(x**2 + y**2)
            the = atan(y, x)
            d = par_R0 * (par_R1 * par_R1/par_R0 + par_R0)/(par_R0**2 + par_R1**2) 
            f1 = par_R0 * sin(the) * (par_R1 * par_R1/r + r)/(par_R0**2 + par_R1**2) /d * par_B0
            f2 = gl_S_in%u(i)
            write(1,*) the, f1, f2

        end do

        close(1)

    end subroutine Print_test

    subroutine Print_Solution()
        INTEGER :: i, NNN, j
        real(8) :: r, the, x, y, ux, uy, x1, y1

        NNN = 150
        open(1, file = 'Solution_in.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'F', 'Bx', 'By', 'BB'"

        do i = par_N + 1, 2 * par_N
            x = yzel_x(i)
            y = yzel_y(i)
            do j = 1, NNN
                x1 = x * (1.0 - 1.0 * j/(NNN + 1.0))
                y1 = y * (1.0 - 1.0 * j/(NNN + 1.0))
                if(sqrt(x1 * x1 + y1 * y1) < par_R0 + 0.05) then
                    write(1,*) x1, y1, 0.0, 0.0, 0.0, 0.0
                else
                    ux = Get_v1(x1, y1, gl_S_in)
                    uy = Get_v2(x1, y1, gl_S_in)
                    write(1,*) x1, y1, Get_f(x1, y1, gl_S_in), ux, uy, (ux**2 + uy**2)/(8.0 * par_pi)
                end if
            end do
        end do

        ! do i = 1, NNN
        !     do j = 1, NNN
        !         r = par_R0 + (par_R1 - par_R0) * (i)/(NNN + 1)
        !         the = j * 2.0_8 * par_pi/NNN
        !         x = r * cos(the)
        !         y = r * sin(the)
        !         ux = Get_v1(x, y, gl_S_in)
        !         uy = Get_v2(x, y, gl_S_in)
        !         write(1,*) x, y, Get_f(x, y, gl_S_in), ux, uy, (ux**2 + uy**2)/(8.0 * par_pi)
        !     end do
        ! end do

        close(1)


        open(1, file = 'Solution_out.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'F', 'Ux', 'Uy', 'p'"

        do i = par_N + 1, 2 * par_N
            x = yzel_x(i)
            y = yzel_y(i)
            do j = 1, NNN
                x1 = x * (1.0 + 3.0 * j/(NNN + 1.0))
                y1 = y * (1.0 + 3.0 * j/(NNN + 1.0))
                if(sqrt(x1 * x1 + y1 * y1) < par_R0 + 0.05) then
                    write(1,*) x1, y1, 0.0, 0.0, 0.0, 0.0
                else
                    ux = Get_v1(x1, y1, gl_S_out)
                    uy = Get_v2(x1, y1, gl_S_out)
                    write(1,*) x1, y1, Get_f(x1, y1, gl_S_out) + par_Uinf * x1, ux + par_Uinf, uy, par_Bernully - ((ux+ par_Uinf)**2 + uy**2)/2.0
                end if
            end do
        end do

        ! do i = 1, NNN
        !     do j = 1, NNN
        !         r = par_R1 + (3.0 * par_R1 - par_R1) * (i)/(NNN + 1)
        !         the = j * 2.0_8 * par_pi/NNN
        !         x = r * cos(the)
        !         y = r * sin(the)
        !         ux = Get_v1(x, y, gl_S_out)
        !         uy = Get_v2(x, y, gl_S_out)
        !         write(1,*) x, y, Get_f(x, y, gl_S_out) + par_Uinf * x, ux + par_Uinf, uy, par_Bernully - ((ux+ par_Uinf)**2 + uy**2)/2.0
        !     end do
        ! end do

        close(1)

    end subroutine Print_Solution


    subroutine Print_yzel_all()
        INTEGER :: N1, i

        open(1, file = '_All_yzel.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'"

        N1 = size(yzel_x(:))

        do i = 1, N1
            write(1,*) yzel_x(i), yzel_y(i)
        end do

        close(1)

    end subroutine Print_yzel_all


    subroutine Print_gran_setka(SS)
        TYPE (Setka), intent(in) :: SS
        INTEGER :: i, N, k1, k2, j
        character(len=3) :: name


        write(unit=name,fmt='(i3.3)') SS%area


        N = SS%N
        open(1, file = name // '_Grans.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', ZONE T= 'HP', N= ", 2 * N, ", E =  ", N , ", F=FEPOINT, ET=LINESEG"

        do i = 1, N
            k1 = SS%gran(1, i)
            k2 = SS%gran(2, i)

            write(1,*) yzel_x(k1), yzel_y(k1)
            write(1,*) yzel_x(k2), yzel_y(k2)
        end do

        do j = 0, N - 1
            write(1,*) 2 * j + 1, 2 * j + 2
        end do

        close(1)

    end subroutine Print_gran_setka

    subroutine Print_matrix_real(A)
        real(8), intent(in) :: A(:, :)   !! Matrix
        integer(4) :: mi, ni

        do mi = 1, size(A(:, 1))
            do ni = 1, size(A(1, :))
                write(*,"(F8.2,$)") A(mi,ni)
            end do
            write (*,*) ''
        end do
    end subroutine Print_matrix_real


end module GEOMETRY