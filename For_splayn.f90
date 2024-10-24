module For_splayn
    use STORAGE 
    implicit none 

    contains

    subroutine Init_Splayn(SS, N_, l, r)
        TYPE (Splayn), intent(in out) :: SS
        INTEGER, INTENT(IN) :: N_, l, r   ! Сколько граней будет в сетке

        SS%N = N_

        if(ALLOCATED(SS%x) == .True.) return

        ALLOCATE(SS%x(N_))
        ALLOCATE(SS%f(N_))
        ALLOCATE(SS%a(3, N_ - 1))
        ALLOCATE(SS%M( 3 *(N_ - 1), 3 *(N_ - 1) ))
        ALLOCATE(SS%B( 3 *(N_ - 1) ))

        SS%left = l    
        SS%right = r    
    end subroutine Init_Splayn


    subroutine Splayn_Set_Matrix(SS)   ! Заполняем матрицу
        TYPE (Splayn), intent(in out) :: SS  
        integer :: i, j

        SS%M = 0.0_8
        SS%B = 0.0_8

        ! Заполняем матрицу - значение на правом конце сплайна
        do i = 1, SS%N - 1
            SS%M(i, 3 * (i - 1) + 1) = SS%x(i + 1) - SS%x(i)
            SS%M(i, 3 * (i - 1) + 2) = (SS%x(i + 1) - SS%x(i))**2
            SS%M(i, 3 * (i - 1) + 3) = (SS%x(i + 1) - SS%x(i))**3
            SS%B(i) = SS%f(i + 1) - SS%f(i)
        end do

        ! Приравниваем значения первой производной на правом конце текущего сплайна
        do i = 1, SS%N - 2
            SS%M(i + SS%N - 1, 3 * (i + 1 - 1) + 1) = -1.0_8

            SS%M(i + SS%N - 1, 3 * (i - 1) + 1) = 1.0_8
            SS%M(i + SS%N - 1, 3 * (i - 1) + 2) = 2.0 * (SS%x(i + 1) - SS%x(i))
            SS%M(i + SS%N - 1, 3 * (i - 1) + 3) = 3.0 * (SS%x(i + 1) - SS%x(i))**2

            SS%B(i + SS%N - 1) = 0.0_8
        end do

        ! Приравниваем значения второй производной на правом конце текущего сплайна
        do i = 1, SS%N - 2
            SS%M(i + 2 * SS%N - 3, 3 * (i + 1 - 1) + 2) = -2.0_8
            
            SS%M(i + 2 * SS%N - 3, 3 * (i - 1) + 1) = 0.0_8
            SS%M(i + 2 * SS%N - 3, 3 * (i - 1) + 2) = 2.0
            SS%M(i + 2 * SS%N - 3, 3 * (i - 1) + 3) = 6.0 * (SS%x(i + 1) - SS%x(i))

            SS%B(i + 2 * SS%N - 3) = 0.0_8
        end do

        if(SS%left == 0) then
            SS%M(3 * SS%N - 4, 1) = 1.0_8
        else
            SS%M(3 * SS%N - 4, 2) = 2.0_8
        end if

        if(SS%right == 0) then
            SS%M(3 * SS%N - 3, 3 * (SS%N - 1) - 2) = 1.0_8
            SS%M(3 * SS%N - 3, 3 * (SS%N - 1) - 1) = 2.0 * (SS%x(SS%N) - SS%x(SS%N - 1))
            SS%M(3 * SS%N - 3, 3 * (SS%N - 1)) = 3.0 * (SS%x(SS%N) - SS%x(SS%N - 1))**2
        else
            SS%M(3 * SS%N - 3, 3 * (SS%N - 1) - 2) = 0.0_8
            SS%M(3 * SS%N - 3, 3 * (SS%N - 1) - 1) = 2.0
            SS%M(3 * SS%N - 3, 3 * (SS%N - 1)) = 6.0 * (SS%x(SS%N) - SS%x(SS%N - 1))
        end if

    end subroutine Splayn_Set_Matrix


    subroutine Splayn_Culc_equ(SS)
        TYPE (Splayn), intent(in out) :: SS  
        external :: dgesv
        integer :: rc, i
        real(8) :: B(3 * (SS%N - 1))
        real(8) :: pivot(3 * (SS%N - 1))

        B = SS%B


        call dgesv(3 * (SS%N - 1), 1, SS%M, 3 * (SS%N - 1), pivot, B, 3 * (SS%N - 1), rc)
        ! call dgesv(2, 1, M, 2, pivot, B, 2, rc)

        ! print*, "rc = ", rc
        ! print*, "B = ", B

        do i = 1, SS%N - 1
            SS%a(1, i) = B(3 * (i - 1) + 1)
            SS%a(2, i) = B(3 * (i - 1) + 2)
            SS%a(3, i) = B(3 * (i - 1) + 3)
        end do

    end subroutine Splayn_Culc_equ

  ! pure
    real(8) function Splayn_Get(SS, x) ! Двигаем поверхность-магнитопаузу согласно давлению с двух сторон
        TYPE (Splayn), intent(in) :: SS  
        real(8), intent(in) :: x
        integer :: i, k

        do i = 2, SS%N
            k = i - 1
            if(SS%x(i) > x) EXIT
        end do

        !print*, "k = ", k

        Splayn_Get = SS%f(k) + SS%a(1, k) * (x - SS%x(k)) + SS%a(2, k) * (x - SS%x(k))**2 + &
                + SS%a(3, k) * (x - SS%x(k))**3
    end function Splayn_Get

    subroutine Print_Splayn(SS, num) ! Двигаем поверхность-магнитопаузу согласно давлению с двух сторон
        TYPE (Splayn), intent(in) :: SS 
        integer, intent(in) :: num
        CHARACTER(len = 3) :: name
        integer :: i
        real(8) :: x

        write(unit=name,fmt='(i3.3)') num

        open(2, file = name // '_splayn.txt')
        write(2,*) "TITLE = 'HP'  VARIABLES = 'x', 'f'"

        ! print*, "x1  -  x2", SS%x(1), SS%x(SS%N)
        ! print*, "f1  -  f2", SS%f(1), SS%f(SS%N)

        do i = 1, 10 * SS%N
            x = SS%x(1) + (i - 1) * (SS%x(SS%N) - SS%x(1))/(10 * SS%N - 1)
            write(2,*) x, Splayn_Get(SS, x)
        end do

        ! print*, "Proverka = ", Splayn_Get(SS, SS%x(1)), SS%f(1) 

        close(2)

    end subroutine Print_Splayn

    subroutine Test_Splayn()
        integer :: i, j, N, jj, o1, o2, o3, o4, k_sqlag
        real(8) :: n1, n2, x0, y0, x1, y1, r, nn, u1, u2, ex, ey
        real(8) :: p1, p2, ux, uy, f1, f2, x2, y2, d, the, pp1, pp2, x3, y3

        N = par_N

        ! Сглаживание
        do i = N + 1, 2 * N   
            x0 = yzel_x(i)
            y0 = yzel_y(i)
            r = sqrt(x0 * x0 + y0 * y0)
            the = atan(y0, x0)
            if(the > 0) then
                sur_the(i - N) = the
            else
                sur_the(i - N) = the + 2 * par_pi
            end if
            sur_r(i - N) = r
        end do
        sur_the(1) = 0.0

        o1 = 25
        call Init_Splayn(Splayn_2, N/2 / o1 + 1, 1, 1)
        do i = 1, N/2 / o1 + 1
            Splayn_2%x(i) = sur_the((i - 1) * o1 + N/4 + 1)
            Splayn_2%f(i) = sur_r((i - 1) * o1 + N/4 + 1)
        end do

        o2 = 25
        call Init_Splayn(Splayn_1, N/4 / o2 + 1, 0, 1)
        do i = 1, N/4 / o2 + 1
            Splayn_1%x(i) = sur_the((i - 1) * o2 + 1)
            Splayn_1%f(i) = sur_r((i - 1) * o2 + 1)
        end do

        ! print*, "111 ", sur_the(1), sur_the(N/4 + 1)

        o3 = 25
        call Init_Splayn(Splayn_3, N/4 / o3 + 1, 1, 0)
        do i = 1, N/4 / o3
            Splayn_3%x(i) = sur_the((i - 1) * o3 + 3 * N/4 + 1)
            Splayn_3%f(i) = sur_r((i - 1) * o3 + 3 * N/4 + 1)
        end do

        Splayn_3%x(N/4 / o3 + 1) = 2 * par_pi
        Splayn_3%f(N/4 / o3 + 1) = sur_r(1)



        call Splayn_Set_Matrix(Splayn_1)
        call Splayn_Culc_equ(Splayn_1)
        call Print_Splayn(Splayn_1, 1)

        call Splayn_Set_Matrix(Splayn_2)
        call Splayn_Culc_equ(Splayn_2)
        call Print_Splayn(Splayn_2, 2)

        call Splayn_Set_Matrix(Splayn_3)
        call Splayn_Culc_equ(Splayn_3)
        call Print_Splayn(Splayn_3, 3)



    end subroutine Test_Splayn

end module For_splayn