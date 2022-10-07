function qddot = fun_qddotsspsinglecontact(x,u)

global A B Nx Nu pert MI L m  nx ny tx ty g r lam vars misc alp alpval indic kc lamall xdata lamx lamy val af acal fx fy Mmat2 invM phi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1 = m(1);
m2 = m(2);
m3 = m(3);

L2 = L(1);
L3 = L(2);

r2 = r(2);
r3 = r(3);

MI1 = MI(1);
MI2 = MI(2);
MI3 = MI(3);

%%%%%%%%%%%%%%%%%%%%%%
tht2 = x(3);
tht3 = x(4);
x1   = x(1);
y1   = x(2);
omg2 = x(7);
omg3 = x(8);
vh1 = x(5);
vh1 = x(6);


F1 = u(3);
F2 = u(4);
T1  =  u(1);
T2  =  u(2);

c = -1.6;
%T1  = u(1);
%{
T2  = c*omg2*r2 +u(1);
T3  = c*omg3*r3 +u(2);
T4  = c*omg4*r4 +u(3);
T5  = c*omg5*r5 +u(4);
T6  = c*omg6*r6 +u(5);
T7  = c*omg7*r7 +u(6);
%}
%lam1 = lamy;
%lam2 =lamx;
%{
T2  = -c*omg2*r2 +u(1);
T3  = -c*omg3*r3 +u(2);
T4  = -MI4*alp4  -c*omg4*r4 +u(3);
T5  = -c*omg5*r5 +u(4);
T6  = -c*omg6*r6 +u(5);
T7  = -c*omg7*r7 +u(6);

%}

%T4  =  c*omg4 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mmat=[m2+m3,0,L2.*m3.*sin(tht2)+m2.*r2.*sin(tht2),m3.*r3.*sin(tht3);0, ...
  m2+m3,(-1).*L2.*m3.*cos(tht2)+(-1).*m2.*r2.*cos(tht2),(-1).*m3.* ...
  r3.*cos(tht3);L2.*m3.*sin(tht2)+m2.*r2.*sin(tht2),(-1).*L2.*m3.* ...
  cos(tht2)+(-1).*m2.*r2.*cos(tht2),L2.^2.*m3+MI2+m2.*r2.^2,L2.*m3.* ...
  r3.*cos(tht2+(-1).*tht3);m3.*r3.*sin(tht3),(-1).*m3.*r3.*cos(tht3) ...
  ,L2.*m3.*r3.*cos(tht2+(-1).*tht3),MI3+m3.*r3.^2]






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mmat'
Mmat - Mmat'
issymmetric(Mmat)
eig(Mmat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% phi is for finding qddot with ext lambda calculation 
phi=[F1+(-1).*L2.*m3.*omg2.^2.*cos(tht2)+(-1).*m2.*omg2.^2.*r2.*cos( ...
  tht2)+(-1).*m3.*omg3.^2.*r3.*cos(tht3),F2+(-1).*g.*m1+(-1).*g.*m2+ ...
  (-1).*g.*m3+(-1).*L2.*m3.*omg2.^2.*sin(tht2)+(-1).*m2.*omg2.^2.* ...
  r2.*sin(tht2)+(-1).*m3.*omg3.^2.*r3.*sin(tht3),T1+(-1).*T2+g.*L2.* ...
  m3.*cos(tht2)+g.*m2.*r2.*cos(tht2)+(-1).*L2.*m3.*omg3.^2.*r3.*sin( ...
  tht2+(-1).*tht3),T2+g.*m3.*r3.*cos(tht3)+L2.*m3.*omg2.^2.*r3.*sin( ...
  tht2+(-1).*tht3)]



invM = inv(Mmat)
qddot_invdy = invM*phi'
af = qddot_invdy ;
qddot = qddot_invdy;
 

end






 






