function [y,iter] = UnExplicitSchema2(N, M, l, T, gamma, u, phi)
%Итерационный процесс с нелинейностью на верхнем слое
h = l/N;
tau = T/M;

x = linspace(0,1,N);
t = linspace(0, T, M);

%нач условия
yOld = zeros(1, N);
y = zeros(1, N);

for i = 1:N
    y(i) = u(x(i), 0);
end

iter = 0;

for j = 2:M
    yP = y;
    for k = 1:10 %Итерационный цикл
        yOld = y;
   
        a = zeros(1,N-2);
        b = zeros(1,N-2);
        c = zeros(1,N-2);

        F = zeros(N,1);
        
        for i = 2:N-1
            a(i-1) = tau/(2*h^2) * (...
                (PsiFunc(yOld(i+1),gamma,l) + ...
                PsiFunc(yOld(i),gamma,l)) / RaznOp(yOld(i+1), yOld(i), h));

            b(i-1) = -tau/(2*h^2) * (...
                (PsiFunc(yOld(i+1),gamma,l) + ...
                PsiFunc(yOld(i),gamma,l)) / RaznOp(yOld(i+1), yOld(i), h) + ...
                (PsiFunc(yOld(i),gamma,l) + ...
                PsiFunc(yOld(i-1),gamma,l)) / RaznOp(yOld(i), yOld(i-1), h)) - 1;

            c(i-1) = tau/(2*h^2) * (...
               (PsiFunc(yOld(i),gamma,l) + ...
                PsiFunc(yOld(i-1),gamma,l)) / RaznOp(yOld(i), yOld(i-1), h));
            
            F(i) = -tau*phi(i,j) - yP(i);
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
        
%         alpha = zeros(1,N);
%         beta = zeros(1,N);
% 
%         a1 = zeros(1.N);
% 
%         alpha(1) = c(1)/b(1);
%         beta(1) = F(1)/b(1);
% 
%         for i = 2:N
%             alpha(i) = c(i)/(b(i) - a(i)*alpha(i-1));
%             beta(i) = (F(i) - a(i)*beta(i-1))/(b(i) - a(i)*alpha(i));
%         end
%         
%         y(N) = beta(N);
%         for i = N-1:-1:1
%             y(i) = alpha(i) * y(i+1) + beta(i);
%         end
% 

        y = A \ F;

        if(max(abs(y) - abs(yOld)) < 0.000001)
            %disp(j);
            iter = iter + 1;
            break;
        end
        clearvars A F;
        iter = iter + 1;
    end

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