program pyreball

    !use variables    
    use M_fdata
    
    fdatalocation='/home/dani/FB/git/create/coutput'
    fdatalocation='/home/dani/Fdata_HC-new'
    call load_fdata()

    call getenergy()
    
    !call diagonalize_matrix()

end program pyreball

