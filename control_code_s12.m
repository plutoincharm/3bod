clc; close all

clearvars -except GRFx GRFy GRFxr GRFyr GRFxl GRFyl; 
global A B Nx Nu pert MI L m  nx ny tx ty g r lam vars misc alp alpval indic kc lamall xdata lamx lamy val

%n = 4; % dimensions of system
%m1 = 50.172;  
%m2=7.4;  m3 = 3.411; m4 = 1.073;
%m5 = 7.4;  m6 = 3.411; m7 = 1.073;
%L1= 2*0.3382;
%L1= 0.51264;  
%L2= 0.4418;  L3= 0.4033;  L4=0.06942;
%L5= 0.4418;  L6= 0.4033;  L7=0.06942;

%MI1 = 2.21304881;
%MI1 =  1.95;
%MI2 =0.1038; MI3 = 0.05916;  MI4 = 0.01;
%MI5 = 0.1038; MI6 = 0.05916; MI7 = 0.01;
%r1 = 0.126486;
%{
r1 = 0.19172736;
r2 = 0.33527; r3 = 0.1896;% r4 = 0.0595; 
r4 = 0.0195;
r5 = 0.33527; r6 = 0.1896;%r7 = 0.0595; 
r7 = 0.0195;

%}
%r1 =L1/2;
%r2 = 0.1872; r3 =0.1738;  r4 = 0.044;
%r5 = 0.1872; r6 = 0.1738; r7 = 0.044;



%MI = [MI1;MI2;MI3;MI4;MI5;MI6;MI7 ];
%L = [L1;L2;L3;L4;L5;L6;L7];
%m = [m1;m2;m3;m4;m5;m6;m7];
%ra = [r1;r2;r3;r4;r5;r6;r7];
g = 9.81; % gravity
Nx = 18;
Nu  = 6;
Tf = 0.56;
%Tf = 1;
dt = 0.01;
Nt = round(Tf/dt)+1;
A = zeros(Nx,Nx);
B = zeros(Nx,Nu);
pert = 0.0001;
nx = 0;tx = 1;
ny = 1;ty = 0;
frame = 20;

dynamics_midpoint = @(x,u,dt) x + fun_xdot(x + fun_xdot(x,u,dt)*dt/2,u,dt)*dt;

vars = fun_data();

% initial conditions
%BV     =  readmatrix("BV.xlsx");
%vel    =  readmatrix("vel.xlsx");
BC      =  readmatrix("BC.xlsx");
omg     =  readmatrix("omg.xlsx");
%alpval  =  readmatrix("alp.xlsx");
tht1=BC(1,frame);tht2=BC(2,frame);tht3=BC(3,frame);tht4=BC(4,frame);tht5=BC(5,frame);tht6=BC(6,frame);tht7=BC(7,frame);hx=BC(8,frame) ;hy=BC(9,frame);
%omg1 = 0.1;omg2 = 0;omg3 = 0.2;omg4 = 0.0;omg5 = 0.5;omg6 = 0.5;omg7 = 0.5;vhx =-0.8;vhy = 0.2;
omg1 = omg(1,frame); omg2 =omg(2,frame); omg3 = omg(3,frame);  omg5 = omg(5,frame); omg6 = omg(6,frame); omg7 = omg(7,frame);vhx =  omg(8,frame); vhy = omg(9,frame);
%omg4 = 0;
 omg4 = omg(4,frame);


%omg1 = vel(1,5);omg2 = vel(2,5); omg3 =  vel(3,5); omg4 = vel(4,5);  omg5 =  vel(5,5);
%omg6 =   vel(6,5);  omg7 =   vel(7,5); vhx=  vel(8,5); vhy =  vel(9,5);
 x0 = [tht1;tht2;tht3;tht4;tht5;tht6;tht7;hx;hy;omg1;omg2;omg3;omg4;omg5;omg6;omg7;vhx;vhy];

% goal
thtf1=BC(1,76);thtf2=BC(2,76);thtf3=BC(3,76);thtf5=BC(5,76);thtf6=BC(6,76);thtf7=BC(7,76);hfx=BC(8,76);hfy=BC(9,76);
omgf1 = omg(1,76);omgf2 = omg(2,76); omgf3 =  omg(3,76); omgf5 = omg(5,76);
omgf6 =   omg(6,76);  omgf7 =  omg(7,76); vhfx=  omg(8,76); vhfy =  omg(9,76);
%omgf4 = 0; 
omgf4 = omg(4,76); 
thtf4=BC(4,76);

xf = [thtf1;thtf2;thtf3;thtf4;thtf5;thtf6;thtf7;hfx;hfy;omgf1;omgf2;omgf3;omgf4;omgf5;omgf6;omgf7;vhfx;vhfy];

% costs
%Q = 1e-5*eye(4);
Q =  1e-5*eye(Nx);
Qf = 25*eye(Nx);
R = 5*1e-7*eye(Nu);
I = eye(Nu);
e_dJ = 1e-12;

%simulation
%dt = 0.1;
%tf = 1;
%N = floor(tf/dt);
%t = linspace(0,tf,N);
%iterations = 100;

% initialization
%u = rand(Nu,Nt-1)*15;
%u = ones(Nu,Nt-1)*10;
u = zeros(Nu,Nt-1);
uf  = ones(Nu,Nt-1);
x = zeros(Nx,Nt);
x_prev = zeros(Nx,Nt);
x(:,1) = x0;

%A = fun_amat(x(:,1),u,dt);
%B = fun_bmat(x(:,1),u,dt);


%pause()
xdata = x0;
lamx = [];
lamallx = [lamx];
lamy = [];
lamally = [lamy];
%alp = zeros(Nx/2);%
alp =alpval(:,frame); 
%alp(4) = 0;
% first roll-out
for k = 2:Nt-1
         x(:,k) = dynamics_midpoint(x(:,k-1),u(:,k-1),dt)
        %time(k-1) = (k-1)*dt + 0.54 ;
        %alp = zeros(Nx/2);
        alp =  alpval(:,frame+k-1);
        % alp(4) = 0;
        thetval = BC(:,frame+k-1);
        omgval = omg(:,frame+k-1);
        xdata = [thetval;omgval];
        lamallx = [lamallx,lamx];
        lamally = [lamally,lamy];
        k;
        %lamx(k-1)  =  lamall(1,k-1) + lamall(3,k-1);
        % lamy(k-1)  =  lamall(2,k-1) + lamall(4,k-1);
        pause()           
        % fc() 
end
tdt = 0.54:0.01:0.73;
tt = 0.23:dt:1.58;

%{
figure;
plot(tt,GRFxl,'b-','LineWidth',1);
grid on;
hold on;
plot(tt,GRFx,'r-','LineWidth',1);
grid on;
hold on;
plot(tdt,lamallx,'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('GRF_x (N) \rightarrow');
legend('2D Lagrangian','Experimental');
%
 



figure;
plot(tt,GRFyl,'b-','LineWidth',1);
grid on;
hold on;
plot(tt,GRFy,'r-','LineWidth',1);
grid on;
hold on;
plot(tdt,lamally,'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('GRF_y (N) \rightarrow');
legend('2D Lagrangian','Experimental');
%}






% original cost
J = 0;                                              
for k = 1:Nt-1
    J = J + (x(:,k)-xf)'*Q*( x(:,k)-xf) + (u(:,k))'*R*(u(:,k)) %+(u(:,k)-uf(:,k))'*R*(u(:,k)-uf(:,k))
end 
disp('Original cost:')
J = 0.5*(J + (x(:,Nt)-xf)'*Qf*(x(:,Nt)-xf)) 
val = x;


%{

pause()
disp(' ILQR starts--------------------------------------- ')

%%%%%%%%%%%%%%%% ILQR Algorithm  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   
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
           gu = r + B'*p(:,k+1)   % gu = dg/du
    
          %iLQR (Gauss-Newton) version
          %Hessian
             Gxx = Q + A'*(P(:,:,k+1))*A;
             Guu = R + B'*(P(:,:,k+1)+ mu_reg*eye(Nx))*B
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
       A
       B
      %pause()
    end
    disp('ITERATIPONSSSSSSSSSSSSSSSSSSSSS- cOMPLETED')
       k
       pause()
    
  
%%%% End of Backward Pass %%%%%
     %alp = zeros(Nx/2);%
    alp =alpval(:,frame) 
    %alp(4) = 0; 
    %Forward rollout with line search
    xn(:,1) = x(:,1);
    alpha = 1.0;
    indic = 15;
    xdata = x0;   
       
   for kc = 1:(Nt-1)
        un(:,kc) = u(:,kc) - alpha*d(:,kc) - (K(:,:,kc)*(xn(:,kc)-x(:,kc)));
        xn(:,kc+1) = dynamics_midpoint(xn(:,kc),un(:,kc),dt);
        thetval = BC(:,frame+kc);
        omgval = omg(:,frame+kc);
        xdata = [thetval;omgval];
        %alp = zeros(Nx/2);
        alp =  alpval(:,frame+kc);
       % alp(4) = 0;
        disp('ITERATIPONSSSSSSSSSSSSSSSSSSSSS')
       kc
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
 liter =1;
    while (isnan(Jn) || Jn > (J - 1e-2*alpha*dJ))%&&(liter<=15)
        alpha = 0.5*alpha
         xdata = x0;   
         alp =alpval(:,frame) 
        for k = 1:(Nt-1)
            un(:,k) = u(:,k) - alpha*d(:,k) - (K(:,:,k)*(xn(:,k)-x(:,k)));
            xn(:,k+1) = dynamics_midpoint(xn(:,k),un(:,k),dt);
            thetval = BC(:,frame+k);
            omgval = omg(:,frame+k);
            xdata = [thetval;omgval];
            alp =  alpval(:,frame+k);
            %xn(:,k+1) = dynamics_rk4(xn(:,k),un(k)
            liter = liter + 1;
        end
        %Jn = cost(xn,un);
        Jn = 0;
        for k = 1:Nt-1
            Jn = Jn + (xn(:,k) - xf)'*Q*(xn(:,k) - xf) + un(:,k)'*R*un(:,k);
        end
     Jn = 0.5*(Jn + (xn(:,Nt) - xf)'*Qf*(xn(:,Nt) - xf))
     % pause()
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






