function dydt = LGalt_model(~,y,b_max,d_min,bslope,dslope,cr,Rcull)
dydt = zeros(size(y));

% variables
R = y(1);

dydt(1) = (b_max - bslope*R)*R - (d_min + dslope*R)*R - cr*(max((R-Rcull),0));
