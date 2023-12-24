subroutine diagonalize_matrix(matrix, eigenvalues)
    implicit none
    integer, parameter :: n = 3 ! Tamaño de la matriz (3x3)
    real, intent(in) :: matrix(n,n)
    real, intent(out) :: eigenvalues(n)
    ! Declaraciones para la diagonalización de la matriz (pueden variar según tu compilador)
    ! ...

    ! Aquí deberías colocar el código para diagonalizar la matriz
    ! ...

    ! Supongamos que 'eigenvalues' contiene los valores propios
    ! y la matriz 'matrix' ha sido diagonalizada

    eigenvalues = [1.0, 2.0, 3.0] ! Valores propios de ejemplo
end subroutine diagonalize_matrix
