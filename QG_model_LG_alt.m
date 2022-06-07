function dydt = QG_model_LG_alt(~,y,d,bslope,dslope,h2,V,to_slope,cr,Rcull)
dydt = zeros(size(y));

% variables
R = y(1);
b = y(2);

fit_grad = 1 - 2*to_slope*b;
d = to_slope*b^2;

dydt(1) = (b - bslope*R)*R - (d + dslope*R)*R - cr*(max((R-Rcull),0));
dydt(2) = V*h2*fit_grad;