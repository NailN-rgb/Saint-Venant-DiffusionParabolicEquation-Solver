clc, clear all;
%Определим константы, порядок разбиения и
%               ширину рассматриваемого канала.


l = 1;
gamma = 0.1; %при 50 и 0.001 разница между итерационным
    % и  методами
T = 1;

N = 20;
M = 500;

u = @(x,t) (t+0.01)*((l-x).*(x+1)/l);

h = l/N;
tau = T/M;

x = 0:h:l;
t = 0:tau:T;


mu1 = @(t) 1 * t + 0.01;
mu2 = @(x) ((l-x).*(x+1))/l;

% [x,t] = meshgrid(0:h:l, 0:tau:T);
% 
% z = t.*((l-x).*(x+1)/l);
% mesh(x,t,z);


dmu1 = @(t) 1;
dmu2 = @(x) 1 - 2*x/l - 1/l;
d2mu2 = @(x) -2/l;



phi = zeros(N,M);

x = 0:h:l;
t = 0:tau:T;

for i = 1:N+1
    for j = 1:M+1
         phi(i,j) = FCalc(mu1, mu2, dmu1,dmu2,d2mu2, gamma,l, x(i), t(j));
    end
end

result2 = ExplicitSchema(N, M, l, T, gamma, u, phi);


x = linspace(0,l,N);
t = 0:tau:T;

f = u(x,T);

plot(x, result2, 'y');
errExpl = max(abs((result2 - f)));
disp(errExpl);
