function f = FCalc(mu1, mu2, dmu1, dmu2, d2mu2, gamma,l, x, t)

m1 = mu1(t);
m2 = mu2(x);

dm1 = dmu1(t);
dm2 = dmu2(x);
d2m2 = d2mu2(x);

func = @(x,t) dm1 .* m2;
     

func2 = @(x,t) (abs(m1 .* dm2)).^(1/2) .* sign(m1 .* dm2) ...
        .* dPsidx(m1 .* m2, m1 .* dm2, gamma, l);

func3 = @(x,t) PsiFunc(m1 .* m2,gamma, l) ...
     .* 1/2 .* abs(m1 .* dm2).^(-1/2) .* m1 .* d2m2;

f1 = func(x,t);
f2 = func2(x,t);
f3 = func3(x,t);
f = f1 - f2 - f3;

%Правая часть для линейного уравнения psi == 0

% func = @(x,t) dm1 * m2 - 0.5 * (abs(m1 * dm2))^(-0.5) * m1 * d2m2;
% f = func(x,t);

end