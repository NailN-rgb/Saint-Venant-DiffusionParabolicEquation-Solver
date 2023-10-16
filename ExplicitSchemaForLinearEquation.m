%Функция вычисляющая реешение уравнения без учета функции пси
function res = ExplicitSchemaForLinearEquation(N, M, l, T, u, phi)

h = l/N;
tau = T/M;

x = 0:h:l;
t = 0:tau:T;

yOld = zeros(1, N + 1);
y = zeros(1, N + 1);

cor = zeros(N + 1,M + 1);

for i = 1:N
    y(i) = u(x(i), 0);
end


for j = 2:M+1
    yOld = y;
    y(1) = u(0, t(j));
    y(N+1) = u(l, t(j));
    for i = 2:N
        y(i) = tau/(2*h) * (-RaznOp(yOld(i), yOld(i-1), h) + ...
            RaznOp(yOld(i+1), yOld(i),h) - ...
            RaznOp(yOld(i), yOld(i-1),h) + ...
            RaznOp(yOld(i+1),yOld(i),h)) + ...
            tau * phi(i,j) + yOld(i);
    end
    cor(:,j) = y(:);
end

res = y;
end