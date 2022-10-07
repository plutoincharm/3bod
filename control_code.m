clear; clc; close all

%n = 4; % dimensions of system
m1=2;m2=3;m3=2;m4=3;m5=2;m6=3;m7=2;m8=3;m9=3;
L1=2;L2=3;L3=2;L4=3;L5=2;L6=3;L7=2;L8=3;L9=3;
g = 9.8; % gravity
Nx = 18;
Nu  = 9;
Tf = 1;
dt = 0.1;
Nt = round(Tf/dt)+1;



pert = 0.001;
A = zeros(18,18);



  qdot = x0;
adc = (dt/2)*{x(10);x(11);x(12);x(13);x(14);x(15);x(16);x(17);x(18);f21(x,u,dt);f22(x,u,dt);f23(x,u,dt);f24(x,u,dt);f25(x,u,dt);f26(x,u,dt);f27(x,u,dt);f28(x,u,dt);f29(x,u,dt)};
%%% For finding A matrix 
v=Nx;h=Nx;%%% v is row and h is column
for i = 1:(Nx/2)
        h = Nx;
        for j = 1:(Nx/2)
            xnew = x + adc; % is this constant ??
            xplus = x;
            xplus(j) =  x(j) + pert; %%% pert if iside diffferentiations
            xplusr = xnew;
            xplusr(h) =  xnew(h) + pert; %%% pert if iside diffferentiations
            xplusr(j) =  xnew(j) + pert*dt/2;
            %f = xnew(10);
            f = xnew((Nx/2)+j);
            F = x(j) + f*dt;
            qddot = fun_qddot(xplus,u,dt)
            fplus = x((Nx/2)+j) + qddot(j)*dt/2;  %  f22 (xplus,u,dt)*dt/2;%%%% select function
            Fplus = xplus(j) + fplus*dt; 

            fr = qddot(xnew,u,dt);
            Fr = x(h) + fr*dt;
            fplusr =qddot(xplusr,u,dt); 
            Fplusr = x(h)+pert + fplusr(h)*dt; 

            A(i,j) = (Fplus - F)/pert;
            A(g,h) = (Fplusr - Fr)/pert;
            xplus =0;xplusr =0;
            h = h-1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
            


%{
            f = {xnew(10:18);f21(xnew,u,dt);f22(xnew,u,dt);f23(xnew,u,dt);f24(xnew,u,dt);f25(xnew,u,dt);f26(xnew,u,dt);f27(xnew,u,dt);f28(xnew,u,dt);f29(xnew,u,dt)};
            F =  qdot  + f*dt;
            fplus = {xplus(10:18);f21(xplus,u,dt);f22(xplus,u,dt);f23(xplus,u,dt);f24(xplus,u,dt);f25(xplus,u,dt);f26(xplus,u,dt);f27(xplus,u,dt);f28(xplus,u,dt);f29(xplus,u,dt)};
            Fplus = qdot  + fplus*dt; 
            A(i,j) = (Fplus - F)/pert;
            xplus =0;

%}
            
        end  

        v = v-1;

end

xdot = @(x,u,dt)[x(10:18)'; fun_qddot(x,u,dt)]

        
dynamics_midpoint = @(x,u,dt) x + fc(x + fc(x,u)*dt/2,u)*dt;
%{
A_midpoint = @(x,dt) [(1 + g*cos(x(1))*(dt^2)/(2*l)) (dt - mu*(dt^2)/(2*m*l^2));
                      (g*cos(x(1) + x(2)*dt/2)*dt/l - mu*g*cos(x(1))*(dt^2)/(2*m*l^3)) (1 + g*cos(x(1) + x(2)*dt/2)*(dt^2)/(2*l) - mu*dt/(m*l^2) + (mu^2)*(dt^2)/(2*(m^2)*l^4))];

B_midpoint = @(x,dt) [(dt^2)/(2*m*l^2); 
                      (-mu*(dt^2)/(2*(m^2)*l^4) + dt/(m*l^2))];

%}




Bmat = @(x,dt) [0,0;0,0;(-1).*L2.^2.*m2.*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1).* ...
  L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1),( ...
  L2.^2.*m2+L1.*L2.*m2.*cos(x(2))).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1) ...
  .*L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1);( ...
  L2.^2.*m2+L1.*L2.*m2.*cos(x(2))).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1) ...
  .*L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1),(( ...
  -1).*L1.^2.*m1+(-1).*L1.^2.*m2+(-1).*L2.^2.*m2+(-2).*L1.*L2.*m2.* ...
  cos(x(2))).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+ ...
  L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1)];


% initial conditions
tht1=deg2rad(20);tht2=deg2rad(35);tht3=deg2rad(20);tht4=deg2rad(35);tht5=deg2rad(20);tht6=deg2rad(35);tht7=deg2rad(35);hx=0.5;hy=1.5;
omg1 = 1;omg2 = 1;omg3 = 1;omg4 = 1;omg5 = 1;omg6 = 1;omg7 = 1;vhx = 1;vhy = 1;
x0 = [tht1;tht2;tht3;tht4;tht5;tht6;tht7;hx;hy;omg1;omg2;omg3;omg4;omg5;omg6;omg7;vhx;vhy];

% goal
thtf1=deg2rad(20);thtf2=deg2rad(35);thtf3=deg2rad(20);thtf4=deg2rad(35);thtf5=deg2rad(20);thtf6=deg2rad(35);thtf7=deg2rad(35);hfx=0.5;hfy=1.5;
omgf1 = 1;omgf2 = 1;omgf3 = 1;omgf4 = 1;omgf5 = 1;omgf6 = 1;omgf7 = 1;vhx = 1;vhy = 1;
xf = [thtf1;thtf2;thtf3;thtf4;thtf5;thtf6;thtf7;hfx;hfy;omgf1;omgf2;omgf3;omgf4;omgf5;omgf6;omgf7;vhfx;vhfy];

% costs
%Q = 1e-5*eye(4);
Q =  1e-5*eye(4);
Qf = 15*eye(4);
R = 5*1e-4*eye(2);
I = eye(4);

e_dJ = 1e-12;

% simulation
%dt = 0.1;
%tf = 1;
%N = floor(tf/dt);
%t = linspace(0,tf,N);
%iterations = 100;

% initialization
u = zeros(Nu,Nt-1);
x = zeros(Nx,Nt);
x_prev = zeros(Nx,Nt);
x(:,1) = x0;


% first roll-out
for k = 2:Nt-1
        x(:,k) = dynamics_midpoint(x(:,k-1),u(:,k-1),dt);
        %x(:,k-1)
        %fc(x(:,k-1) + fc(x(:,k-1),u(:,k-1))*dt/2,u(:,k-1))
        % fc() 
end

% original cost
J = 0;
for k = 1:Nt-1
    J = J + (x(:,k)-xf)'*Q*( x(:,k)-xf) + u(:,k)'*R*u(:,k);
end
disp('Original cost:')
J = 0.5*(J + (x(:,Nt)-xf)'*Qf*(x(:,Nt)-xf))



%%%%%%%%%%%%%%%% ILQR Algorithm  %%%%%%%%%%%%%%%%%%%%
p = ones(Nx,Nt);
P = zeros(Nx,Nx,Nt);
%d = ones(Nu,Nu,Nt-1);
d = ones(Nu,Nt-1);
K = zeros(Nu,Nx,Nt-1);
%pdim = ones(Nx,Nu,Nt);
dJ = 1.0;  % change in cost

xn = zeros(Nx,Nt);
un = zeros(Nu,Nt-1);
% func g(dx,du) is perturbation of val func
% grad- g/ hessian-G of change in value fun
gx = zeros(Nx);
gu = zeros(Nu);
Gxx = zeros(Nx,Nx);
Guu = zeros(Nu,Nu);
Gxu = zeros(Nx,Nu);
Gux = zeros(Nu,Nx);

iter = 0;
while max(abs(d(:))) >  1e-3
    
    iter = iter +  1 

 %%%%% Backward Pass %%%%%
    dJ = 0.0;
    p(:,Nt) = Qf*(x(:,Nt)-xf);     %%% P is vx
    P(:,:,Nt) = Qf;                %%% P is vxx
    mu_reg = 0;
    for k = (Nt-1):-1:1
          %Calculate derivatives
           q = Q*( x(:,k)-xf);     % lx
           r = R*u(:,k);           % lu
            
            A = Amat(x(:,k), u(:,k))
            B = Bmat(x(:,k), u(:,k))

           %gradient of change in val fn
            gx = q + A'*p(:,k+1);% gx = dg/dx  
            gu = r + B'*p(:,k+1);% gu = dg/du
    
          %iLQR (Gauss-Newton) version
          %Hessian
             Gxx = Q + A'*(P(:,:,k+1))*A;
             Guu = R + B'*(P(:,:,k+1)+ mu_reg*eye(Nx))*B;
             Gxu = A'*P(:,:,k+1)*B;
             Gux = B'*(P(:,:,k+1) + mu_reg*eye(Nx))*A;     
             
             %beta = 0.1;
             log = issymmetric([Guu]);
             eigv = eig([Guu]);

          if any(eig(Guu)<0)
            mu_reg = mu_reg + 1;
            k = Nt-1;
            disp('regularized')
          end
        %{
              while (log==0) || all(eigv < 0) 
                    Gxx = Gxx + A'*beta*I*A
                    Guu = Guu + B'*beta*I*B
                    Gxu = Gxu + A'*beta*I*B
                    Gux = Gux + B'*beta*I*A
                    beta = 2*beta
                    %display("regularizing G")
                    display(beta)
                    log = issymmetric([Gxx Gxu; Gux Guu]);
                    eigv = eig([Gxx Gxu; Gux Guu]);
              end
         %}
            d(:,k) = Guu\gu;  % feedforward term
            K(:,:,k) = Guu\Gux; % feedback gain term
    
             p(:,k) = gx - K(:,:,k)'*gu + K(:,:,k)'*Guu*d(:,k) - Gxu*d(:,k);
             P(:,:,k) = Gxx + K(:,:,k)'*Guu*K(:,:,k) - Gxu*K(:,:,k) - K(:,:,k)'*Gux;
             dJ = dJ +  gu'*d(:,k);

          
    end
    disp('EOBP')
  
%%%% End of Backward Pass %%%%%
      
    %Forward rollout with line search
    xn(:,1) = x(:,1);
    alpha = 1.0;
   for k = 1:(Nt-1)
        un(:,k) = u(:,k) - alpha*d(:,k) - (K(:,:,k)*(xn(:,k)-x(:,k)));
        xn(:,k+1) = dynamics_midpoint(xn(:,k),un(:,k),dt);
    end
    disp('EOFP')
    Jn = 0;
    for k = 1:Nt-1
        Jn = Jn + (xn(:,k)-xf)'*Q*(xn(:,k)-xf) + un(:,k)'*R*un(:,k);
    end
   Jn = 0.5*(Jn + (xn(:,Nt)-xf)'*Qf*(xn(:,Nt)-xf))
    
    %{
     while isnan(Jn) || Jn > (J - 1e-2*alpha*dJ)
        alpha = 0.5*alpha
        for k = 1:(Nt-1)
            un(:,k) = u(:,k) - alpha*d(:,k) - (K(:,:,k)*(xn(:,k)-x(:,k)));
            xn(:,k+1) = dynamics_midpoint(xn(:,k),un(:,k),dt);
            %xn(:,k+1) = dynamics_rk4(xn(:,k),un(k)
        end
        %Jn = cost(xn,un);
        Jn = 0;
        for k = 1:Nt-1
            Jn = Jn + (xn(:,k) - xf)'*Q*(xn(:,k) - xf) + un(:,k)'*R*un(:,k);
        end
     Jn = 0.5*(Jn + (xn(:,Nt) - xf)'*Qf*(xn(:,Nt) - xf))
    end
 
    %}
   
    J = Jn;
    x = xn;
    u = un;
   %if iter > 5
    %   break
    %end
  end

   
    



















