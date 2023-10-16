%Функция численно вычисляет производную функции пси, используется 
%                           для вычисления правой части
function dpsi = dPsidx(y, dydx, gamma, l)

first = (5/3) * ((l* y).^2).^(1/3) .* dydx .* (2 * y + l).^(2/3);

second = (4/3) * (2 * y + l).^(-1/3) .* dydx .* ((l * y).^5).^(1/3);

third = (2 * y + l).^(4/3);

dpsi = gamma * ((first - second)./third);

% dpsi = sign(y);

end