function xdot = fun_xdot(x,u,dt)

global A B Nx Nu pert MI L m  nx ny tx ty g ra lam vars alp alpval indic kc lamall xdata lamx lamy

xdot = [x(5:8); fun_qddotsspsinglecontact(x,u)];
%xdot = fun_qddotssp11(x,u,dt);
