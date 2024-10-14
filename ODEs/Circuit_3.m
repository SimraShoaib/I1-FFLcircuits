function dx= Circuit_3(t,x,p,par)

dx = zeros(9,1);

%%%%molecule species

STAR = x(1);  
THS = x(2);       
Y = x(3);   
Y_act = x(4);
RBP = x(5);
Pz_act  = x(6);
Z = x(7);   
GFP = x(8);
Z_rep = x(9);

Z_free = Z - Z_rep;

Pz = par.P_z  - Pz_act;

dx(1) = par.P_x*p(1) - p(12)*STAR*Pz - p(6)*STAR + p(13)*Pz_act; 

dx(2) = par.P_x*p(1) - p(14)*THS*Y - p(7)*THS;

dx(3) = p(2)*par.P_y - p(14)*THS*Y - p(8)*Y;

dx(4) = p(14)*THS*Y - p(8)*Y_act; 

dx(5) = p(3)*Y_act - p(15)*RBP*Z - p(9)*RBP + p(17)*Z_rep + p(10)*Z_rep;

dx(6) = p(12)*STAR*Pz - p(13)*Pz_act;

dx(7) = p(4)*Pz_act - p(15)*RBP*Z - p(10)*Z +  p(9)*Z_rep; 

dx(8) = p(5)*Z_free - p(11)*GFP;

dx(9) =  p(15)*RBP*Z  - p(18)*Z_rep - p(9)*Z_rep - p(10)*Z_rep; 
