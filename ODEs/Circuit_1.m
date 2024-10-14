function dx= Circuit_1(t,x,p,par)

dx = zeros(10,1);

%%%%molecule species

STAR = x(1);  
THS = x(2);       
TetR = x(3);
aTc = x(4);
Pz_rep = x(5);     
Z = x(6);
Z_act = x(7);
GFP = x(8);
Py_act = x(9);
aTc_TetR = x(10);


Pz = par.P_z - Pz_rep;
Py = par.P_y - Py_act;

dx(1) = par.P_x*p(1) - p(10)*STAR*Py - p(5)*STAR + p(11)*Py_act; 

dx(2) = par.P_x*p(1) - p(12)*THS*Z - p(8)*THS;

dx(3) = p(2)*Py_act - p(13)*TetR*Pz - p(14)*TetR*aTc - p(6)*TetR + p(15)*Pz_rep + p(16)*aTc_TetR;

dx(4) = -p(14)*TetR*aTc + p(16)*aTc_TetR;

dx(5) = p(13)*TetR*Pz - p(15)*Pz_rep;

dx(6) = p(3)*Pz - p(12)*THS*Z - p(7)*Z;

dx(7) = p(12)*THS*Z - p(7)*Z_act; 

dx(8) = p(4)*Z_act - p(9)*GFP;

dx(9) = p(10)*STAR*Py - p(11)*Py_act;

dx(10) = p(14)*TetR*aTc - p(16)*aTc_TetR;

