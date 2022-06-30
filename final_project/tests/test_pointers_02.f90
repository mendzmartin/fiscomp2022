! make clean && make test_pointers_02.o && ./test_pointers_02.o
program test_pointers_02
    use module_precision
    implicit none
    type person
        integer(sp)           :: edad
        character (len=50)    :: nomb
        type(person), pointer :: hacia ! apunta a un tipo person
    end type person
    type (person), pointer    :: inic,sig
    type (person)             :: temp
    integer(sp)               :: ios

    temp%edad = 0_sp
    temp%nomb = ' '
    nullify(temp%hacia)
    allocate(inic) ! reserva memoria para el primer elemento
    if (.NOT. ASSOCIATED(inic)) stop 'error de memoria'
    ! lectura de lista
    sig => inic ! puntero => destino. El puntero es un alias
    ! para su destino
    do
        sig = temp ! inicializa sig , puede suprimirse
        read(*,*,iostat=ios) sig%edad,sig%nomb ! lee sig
        if (ios < 0_sp) exit
        if (sig%edad < 0_sp) exit
        allocate(sig%hacia)
        if (.NOT. ASSOCIATED(sig%hacia)) stop 'error de memoria'
        sig => sig%hacia
    end do
    ! escribe lista
    sig => inic
    do
        if (.NOT. ASSOCIATED(sig%hacia)) exit
        write (*,*) sig%edad,sig%nomb
        sig => sig%hacia
    end do
end program test_pointers_02