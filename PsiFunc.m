function psi = PsiFunc(y, gamma, l)


first = ((l * y).^5).^(1/3);
second = (2 * y + l).^(2/3);

psi = gamma * (first./second);
% Тестовая функция

% psi = abs(y);
end