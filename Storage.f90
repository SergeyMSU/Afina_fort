module STORAGE
    implicit none
    public
    real(8), parameter :: par_pi = acos(-1.0_8) 
    real(8), parameter :: par_sqrtpi = sqrt(par_pi)


    real(8), parameter :: par_R0 = 1.0_8  ! Внутренняя сфера (Земля)
    real(8), parameter :: par_R1 = 20.0_8  ! Внешняя сфера (Магнитосфера)
    real(8), parameter :: par_Uinf = -1.0_8 
    real(8), parameter :: par_Bernully = 15.5_8 
    real(8), parameter :: par_move = 0.01_8 
    real(8), parameter :: par_otstup = 0.5_8 
    real(8) :: par_sglag = 0.01_8 !0.3_8    ! Сглаживание
    real(8), parameter :: par_B0 = 3551.22_8 


    integer :: par_N
    real(8), allocatable :: yzel_x(:)  ! Все узлы в программе
    real(8), allocatable :: yzel_y(:)  ! Все узлы в программе

    ! Далее для магнитопаузы, т.к. мы её двигаем
    ! От N + 1 до 2N - поверхность Магнитопаузы
    real(8), allocatable :: vel_gran(:)  ! Скорость грани
    real(8), allocatable :: vel_yzel_x(:)  ! Скорость узла
    real(8), allocatable :: vel_yzel_y(:)  ! Скорость узла

    real(8), allocatable :: sur_r(:)  ! Все узлы в программе
    real(8), allocatable :: sur_r2(:)  ! Все узлы в программе
    real(8), allocatable :: sur_the(:)  ! Все узлы в программе


    TYPE :: Setka
        ! Данные, которые не меняются
        integer :: area
        ! 1 - внутренняя область
        ! 2 - внешняя область

        integer :: N  ! Сколько граней в подсетке
        integer, allocatable :: gran(:, :)  ! (2, :)  Грани данной подсетки - подзадачи

        integer, allocatable :: gran_info(:)  ! (:)  Грани данной подсетки - подзадачи
        ! 1 - Внутренняя сфера - Земля
        ! 2 - Внешняя сфера - Магнитосфера

        LOGICAL, allocatable :: if_u_0(:)  ! Задаётся ли значение функции на грани изначально
        LOGICAL, allocatable :: if_un_0(:)  ! (:)  Задаётся ли значение нормальной скорости на грани изначально

        ! Данные, которые нужно пересчитывать при движении сетки
        real(8), allocatable :: center_gran_x(:)  ! (:)  Грани данной подсетки - подзадачи
        real(8), allocatable :: center_gran_y(:)  ! (:)  Грани данной подсетки - подзадачи
        real(8), allocatable :: normal_gran_x(:)  ! (:)  Грани данной подсетки - подзадачи
        real(8), allocatable :: normal_gran_y(:)  ! (:)  Грани данной подсетки - подзадачи
        real(8), allocatable :: len_gran(:)  ! (:)  Длина грани
        real(8), public, allocatable :: u(:)  ! Значение функции на грани
        LOGICAL, allocatable :: if_u(:)  ! Значение функции на грани
        real(8), allocatable :: un(:)  ! (:)  Значение нормальной скорости на грани
        LOGICAL, allocatable :: if_un(:)  ! (:)  Значение нормальной скорости на грани

        ! Матрица
        real(8), allocatable :: MM(:, :)
        real(8), allocatable :: BB(:)

    end TYPE


    TYPE (Setka):: gl_S_in    ! Подседка для внутренней задачи
    TYPE (Setka):: gl_S_out   ! Подсетка для внешней задачи


end module STORAGE