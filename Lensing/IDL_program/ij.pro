;ij are the indices identified with the WHERE function
;N is the number of elements in one side of the array

FUNCTION IJ, ij, N

    N = LONG(N)
    
	j = LONG(FIX(DOUBLE(ij/N)))
	i = LONG(ROUND(DOUBLE(ij - j*N)))
    return, [[i],[j]]
END
