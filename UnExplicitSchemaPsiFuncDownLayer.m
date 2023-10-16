function [y,iter] = UnExplicitSchemaPsiFuncDownLayer(N, M, l, T, gamma, u, phi)
%Итерационный процесс с нелинейностью на верхнем слое
%Функция пси вычисляется на предыдущем слое
h = l/N;
tau = T/M;

x = linspace(0,l,N);
t = linspace(0,T,M);

%нач условия
yOld = zeros(1, N + 1);
y = zeros(1, N + 1);

for i = 1:N
    y(i) = u(x(i), 0);
end


iter = 0;

for j = 2:M-1    
    a = zeros(1,N-2);
    b = zeros(1,N-2);
    c = zeros(1,N-2);

    F = zeros(N,1);
        
    for i = 2:N-1
        a(i-1) = tau/(2*h^2) * (...
            (PsiFunc(y(i+1),gamma,l) + ...
            PsiFunc(y(i),gamma,l)) / RaznOp(y(i+1), y(i), h));

        b(i-1) = -tau/(2*h^2) * (...
           (PsiFunc(y(i+1),gamma,l) + ...
            PsiFunc(y(i),gamma,l)) / RaznOp(y(i+1), y(i), h) + ...
           (PsiFunc(y(i),gamma,l) + ...
            PsiFunc(y(i-1),gamma,l)) / RaznOp(y(i), y(i-1), h)) - 1;

        c(i-1) = tau/(2*h^2) * (...
           (PsiFunc(y(i),gamma,l) + ...
            PsiFunc(y(i-1),gamma,l)) / RaznOp(y(i), y(i-1), h));
            
        F(i) = -tau*phi(i,j) - y(i);
    end
        
        
        
    A = zeros(N,N);

    for i = 1:N
        if i == 1
            A(1,1) = 1;
        elseif i == N
            A(N,N) = 1;
        else
            A(i, i) = b(i-1);
            A(i, i-1) = c(i-1);
            A(i, i+1) = a(i-1);
        end
    end

    F(1) = u(0, t(j));
    F(end) = u(l, t(j));
     
    y = A \ F;
    
    if(j == M/4)
        f = u(x, t(j));
        plot(x,f',x,y);
        legend("Точное", "Приближенное");
        grid on;
        xlabel('x');
        ylabel('u');
        err = max(abs(y-f'));
        title("T/4 err = ", num2str(err))
        disp("a");


    elseif(j == M/2)
        f = u(x, t(j));
        plot(x,f',x,y);
        legend("Точное", "Приближенное");
        grid on;
        xlabel('x');
        ylabel('u');
        err = max(abs(y-f'));
        title("T/2 err = ", num2str(err))
                disp("a");

    elseif(j == 3/4 * M)
        f = u(x, t(j));
        plot(x,f',x,y);
        legend("Точное", "Приближенное");

        grid on;
        xlabel('x');
        ylabel('u');
        err = max(abs(y-f'));
        title("3/4 T err = ", num2str(err))
                disp("a");

    end

end
end
