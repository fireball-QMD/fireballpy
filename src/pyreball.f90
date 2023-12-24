program pyreball
    use variables    

    open(unit=10, file='fireball.in', status='old')
    read(10, *) num1, num2
    close(10)

    call suma()
    !call load_fdata()
    print *, "La suma es:", resultado
end program pyreball

