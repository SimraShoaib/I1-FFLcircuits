function dx= TX_TL_TL_ODE(t,x,p,par)

dx = zeros(8,1);

%%%%molecule species

STAR = x(1);  
THS = x(2);       
Py_act = x(3);   
RBP = x(4);
Z = x(5);   
Z_act = x(6);
Z_rep  = x(7);
GFP = x(8);

Py = par.P_y - Py_act;

dx(1) = par.P_x*p(1) - p(10)*STAR*Py - p(5)*STAR + p(11)*Py_act; 

dx(2) = par.P_x*p(1) - p(12)*THS*Z - p(6)*THS;

dx(3) =  p(10)*STAR*Py -  p(11)*Py_act; 

dx(4) = p(2)*Py_act - p(13)*RBP*Z - p(7)*RBP + p(14)*Z_rep +  p(8)*Z_rep - p(13)*RBP*Z_act;

dx(5) = p(3)*par.P_z - p(13)*RBP*Z - p(12)*THS*Z - p(8)*Z + p(14)*Z_rep +  p(7)*Z_rep ;

dx(6) = p(12)*THS*Z- p(17)*Z_act - p(13)*RBP*Z_act; % Z_act

dx(7) =  p(13)*RBP*(Z+Z_act) - p(16)*Z_rep - p(14)*Z_rep - p(7)*Z_rep - p(8)*Z_rep; 

dx(8) = p(4)*Z_act - p(9)*GFP;
