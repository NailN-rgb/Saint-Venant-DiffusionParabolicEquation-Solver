function [y,iter] = UnExplicitSchema(N, M, l, T, gamma, u, phi)
%Итерационный процесс с нелинейностью на верхнем слое
h = l/N;
tau = T/M;

x = 0:h:l;
t = 0:tau:T;

%нач условия
yOld = zeros(1, N + 1);
y = zeros(1, N + 1);

for i = 1:N
    y(i) = u(x(i), 0);
end

iter = 0;

for j = 2:M
    y(1) = u(0, t(j));
    y(N+1) = u(l, t(j));
    yP = y;
    for k = 1:10 %Итерационный цикл
        yOld = y;
        for i = 2:N 
            first = tau/(2*h) * ((PsiFunc(yOld(i+1), gamma, l) + PsiFunc(yOld(i), gamma, l)) ...
                / RaznOp(yOld(i+1), yOld(i),h)) * RaznOp(y(i+1), y(i),h);
            
            second = tau/(2*h) * ( (PsiFunc(yOld(i), gamma, l) + PsiFunc(yOld(i-1), gamma, l)) ...
                / RaznOp(yOld(i), yOld(i-1),h)) * RaznOp(y(i), y(i-1),h);

            y(i) = first - second + phi(i,j) * tau + yP(i);
        end

        
        if(max(abs(y - yOld)) < 0.000001)
            %disp(j);
            break
        end
    
        iter = iter + 1;
    end

end


end