function y = ExplicitSchema(N, M, l, T, gamma, u, phi)

%yOld, y - массивы из N элементов. Добавить Граничные условия

h = l/N;
tau = T/M;

x = linspace(0,l,N);
t = linspace(0,T,M);

%нач условия
yOld = zeros(1);
y = zeros(1);


for i = 1:N
    y(i) = u(x(i), 0);
end


for j = 2:M
    yOld = y;
    y(1) = u(0, t(j));
    y(N) = u(l, t(j));
    for i = 2:N-1
        RaznUp = RaznOp(yOld(i+1), yOld(i),h);
        RaznDown = RaznOp(yOld(i), yOld(i - 1),h);

        first = tau/(2*h) * RaznUp * (PsiFunc(yOld(i+1), gamma, l) + PsiFunc(yOld(i), gamma, l));
        second = (-1)*(tau/(2*h) * RaznDown * (PsiFunc(yOld(i-1), gamma, l) + PsiFunc(yOld(i), gamma, l)));
        third = + phi(i, j) * tau + yOld(i);
        
        y(i) = first + second + third;

        if(sum(isnan(y)) > 0)
            return;
        end

%         y(i) = tau/(2*h) * RaznOp(yOld(i+1), yOld(i),h) * (PsiFunc(yOld(i+1), gamma, l) + PsiFunc(yOld(i), gamma, l))...
%             -(tau/(2*h) * RaznOp(yOld(i), yOld(i - 1),h) * (PsiFunc(yOld(i-1), gamma, l) + PsiFunc(yOld(i), gamma, l)))...
%             + phi(i, j) * tau ...
%             + yOld(i);
       



    end

    if(j == M/4)
        f = u(x, t(j));
        plot(x,f,x,y);
        legend("Точное", "Приближенное");
        grid on;
        xlabel('x');
        ylabel('u');
        err = max(abs(y-f));
        title("T/4 err = ", num2str(err))
        disp("a");


    elseif(j == M/2)
        f = u(x, t(j));
        plot(x,f,x,y);
        legend("Точное", "Приближенное");
        grid on;
        xlabel('x');
        ylabel('u');
        err = max(abs(y-f));
        title("T/2 err = ", num2str(err))
                disp("a");

    elseif(j == 3/4 * M)
        f = u(x, t(j));
        plot(x,f,x,y);
        legend("Точное", "Приближенное");

        grid on;
        xlabel('x');
        ylabel('u');
        err = max(abs(y-f));
        title("3/4 T err = ", num2str(err))
                disp("a");

    end
end
end