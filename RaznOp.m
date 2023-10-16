function y = RaznOp(yUp, yDown,h)

r = (yUp - yDown) / h;
%y = sqrt((abs(r))) * sign(r);
y = sqrt((abs(r)));

end