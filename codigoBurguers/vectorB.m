function b = vectorB(b,I,alfa,beta,up)
    b(1,1) = up(1);
    b(end,1) = up(end);
    for i=I
        b(i) = (alfa+beta*up(i))*up(i-1)...
              +(1-2*alfa-beta*up(i))*up(i)...
              +alfa*up(i+1);
    end

end