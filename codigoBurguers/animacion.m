function [] = animacion(campoVel,x,t)
    [~,cols] = size(campoVel);
    obj = VideoWriter("burgers1DUpwind","MPEG-4");
    open(obj);
    for i = 1:cols
        plot(x,campoVel(:,i),"--r");
        title(sprintf("Estado del sistema dinámico en t = %0.2f",t(i)))
        axis([0 2 0 3.5]);
        grid on;
        xlabel("x[m]")
        ylabel("u[m/s]")
        writeVideo(obj,getframe);
        pause(0.2);
    end
    close(obj)
end