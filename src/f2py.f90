python module f2py
    interface
        subroutine suma() bind(c)
        end subroutine suma
    end interface
end python module f2py

python module variables ! Nombre del módulo en Python
    real :: num1 ! Declara las variables de módulo para acceder a ellas desde Python
    real :: num2
    real :: resultado
end python module variables


