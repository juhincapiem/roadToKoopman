function A = matrixA(A,I,alfa,beta,up)
    A(1,1) = 1;
    A(end,end) = 1;
    for j = I
        for i = I
            if j == i
                A(j,i-1) = -(alfa+beta*up(i));
                A(j,i) = (1+2*alfa+beta*up(i));
                A(j,i+1) = -alfa;
            end
        end
    end
end