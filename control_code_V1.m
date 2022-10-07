clear; clc; close all

global A B Nx Nu pert MI L m  nx ny tx ty g ra  

%n = 4; % dimensions of system
m1 = 50.172;  m2=7.4;  m3 = 3.411; m4 = 1.073;
m5 = 7.4;  m6 = 3.411; m7 = 1.073;
L1= 2*0.3382;
L2=0.7743;  L3=0.43788;  L4=0.06942;
L5= 0.7743;  L6=0.43788;  L7=0.06942;

MI1 = 2.21304881;
MI2 = 1.29371245; MI3 = 0.182329; MI4 = 1.59264E-05;
MI5 = 1.293712458; MI6 = 0.182329; MI7 = 1.59264E-05;
r1 = 0.126486;
r2 = 0.33527; r3 = 0.1896; r4 = 0.0195;
r5 = 0.33527; r6 = 0.1896; r7 = 0.0195;
MI = [MI1;MI2;MI3;MI4;MI5;MI6;MI7];
L =  [L1;L2;L3;L4;L5;L6;L7];
m = [m1;m2;m3;m4;m5;m6;m7];
ra = [r1;r2;r3;r4;r5;r6;r7];
g = 9.8; % gravity
Nx = 18;
Nu  = 9;
Tf = 1;
dt = 0.01;
Nt = round(Tf/dt)+1;
A = zeros(Nx,Nx);
B = zeros(Nx,Nu);
pert = 0.001;
nx = 0;
tx = 1;
ny = 1;
ty = 0;


%xdot = @(x,u,dt)[x(10:18)'; fun_qddot(x,u,dt)]


dynamics_midpoint = @(x,u,dt) x + fun_xdot(x + fun_xdot(x,u,dt)*dt/2,u,dt)*dt;


% initial conditions
BV  =   readmatrix("BV.xlsx");
tht1=BV(1,4);tht2=BV(2,4);tht3=BV(3,4);tht4=BV(4,4);tht5=BV(5,4);tht6=BV(6,4);tht7=BV(7,4);hx=BV(8,4);hy=BV(9,4);
omg1 = 0.2;omg2 = 0.2;omg3 = 0.2;omg4 = 0;omg5 = 0.5;omg6 = 0.2;omg7 = 0;vhx =0.2;vhy = 0.2;
x0 = [tht1;tht2;tht3;tht4;tht5;tht6;tht7;hx;hy;omg1;omg2;omg3;omg4;omg5;omg6;omg7;vhx;vhy]

% goal
thtf1=BV(1,5);thtf2=BV(2,5);thtf3=BV(3,5);thtf4=BV(4,5);thtf5=BV(5,5);thtf6=BV(6,5);thtf7=BV(7,5);hfx=BV(8,5);hfy=BV(9,5);
omgf1 = 1;omgf2 = 1;omgf3 = 1;omgf4 = 0;omgf5 = 1;omgf6 = 1;omgf7 = 0;vhfx = 1;vhfy = 1;
xf = [thtf1;thtf2;thtf3;thtf4;thtf5;thtf6;thtf7;hfx;hfy;omgf1;omgf2;omgf3;omgf4;omgf5;omgf6;omgf7;vhfx;vhfy];

% costs
%Q = 1e-5*eye(4);
Q =  1e-5*eye(Nx);
Qf = 15*eye(Nx);
R = 5*1e-4*eye(Nu);
%I = eye(Nu);

e_dJ = 1e-12;

%simulation
%dt = 0.1;
%tf = 1;
%N = floor(tf/dt);
%t = linspace(0,tf,N);
%iterations = 100;

% initialization
%u = rand(Nu,Nt-1)*100;
u = zeros(Nu,Nt-1);
x = zeros(Nx,Nt);
x_prev = zeros(Nx,Nt);
x(:,1) = x0;

%A = fun_amat(x(:,1),u,dt);
%B = fun_bmat(x(:,1),u,dt);



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



%{

%%%%%%%%%%%%%%%% ILQR Algorithm  %%%%%%%%%%%%%%%%%%%%
p = ones(Nx,Nt);
P = zeros(Nx,Nx,Nt);
%d = ones(Nu,Nu,Nt-1);
d = ones(Nu,Nt-1);
K = zeros(Nu,Nx,Nt-1);
%pdim = ones(Nx,Nu,Nt);
dJ = 0.0;  % change in cost

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
   
          %Calculate derivatives of stage cost
           q = Q*( x(:,k)-xf);     % lx
           r = R*u(:,k);       % lu
            
            A = fun_amat(x(:,k),u(:,k),dt);
            B = fun_bmat(x(:,k),u(:,k),dt);

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
          %% 
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
 disp('ITERATIPONSSSSSSSSSSSSSSSSSSSSS')
       k
      % pause()
    end
    disp('ITERATIPONSSSSSSSSSSSSSSSSSSSSS- cOMPLETED')
       k
       pause()
    
  
%%%% End of Backward Pass %%%%%
      
    %Forward rollout with line search
    xn(:,1) = x(:,1);
    alpha = 1.0;
   for k = 1:(Nt-1)
        un(:,k) = u(:,k) - alpha*d(:,k) - (K(:,:,k)*(xn(:,k)-x(:,k)));
        xn(:,k+1) = dynamics_midpoint(xn(:,k),un(:,k),dt);
        disp('ITERATIPONSSSSSSSSSSSSSSSSSSSSS')
       k
      pause()
    end
    disp('EOFP')
     pause() 
    Jn = 0;
    for k = 1:Nt-1
        Jn = Jn + (xn(:,k)-xf)'*Q*(xn(:,k)-xf) + un(:,k)'*R*un(:,k);
    end
   Jn = 0.5*(Jn + (xn(:,Nt)-xf)'*Qf*(xn(:,Nt)-xf))
    
     disp('line search')
     pause() 
 
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
      pause()
    end
        pause() 
 
    J = Jn;
    x = xn;
    u = un;
   %if iter > 5
    %   break
    %end
  end

%}

    




















