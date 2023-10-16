function [y,iter] = UnExplicitSchemaModifiedPsiFunc(N, M, l, T, gamma, u, phi)
%Итерационный процесс с нелинейностью на верхнем слое
%Модифицируем способ вычисления ПСИ
h = l/N;
tau = T/M;

x = 0:h:l;
t = 0:tau:T;

PsiNew = @(y1,y2) 1/2 * (PsiFunc(y1,gamma,l) + PsiFunc(y2,gamma,l));
%Если в исходной формуле ф-я пси вычисляется в точке i, то в PsiNew
%    необходимо подать i-1 and i


%нач условия
y = zeros(1, N + 1);
yOld = zeros(1, N + 1);


for i = 1:N
    y(i) = u(x(i), 0);
end


iter = 0;

for j = 2:M
    y(1) = u(0, t(j));
    y(N+1) = u(l, t(j));
    yP = y;
    yOld = yP; 
    for k = 1:10 %Итерационный цикл
        for i = 2:N 
            first = tau/(2*h) * ( (PsiNew(yOld(i),yOld(i+1)) + PsiNew(yOld(i-1), yOld(i)) ) ...
                / RaznOp(yOld(i+1), yOld(i),h)) * RaznOp(y(i+1), y(i),h);
            
            second = tau/(2*h) * ( (PsiNew(yOld(i-1), yOld(i)) + PsiNew(yOld(i-1), yOld(i)) ) ... %?
                / RaznOp(yOld(i), yOld(i-1),h)) * RaznOp(y(i), y(i-1),h);

            y(i) = first - second + phi(i,j) * tau + yP(i);
        end
        
        yOld = y; 

        if(max(abs(y - yOld)) < 0.000001)
            %disp(j);
            break
        end
        
        iter = iter + 1;
    end

end
end