clc; close all

clearvars -except GRFx GRFy GRFxr GRFyr GRFxl GRFyl; 
global A B Nx Nu pert MI L m  nx ny tx ty g r lam vars misc alp alpval indic kc lamall xdata lamx lamy val af acal fx fy Mmat2 invM phi

g = 9.81; % gravity
Nx = 8;
Nu  = 4;
Tf = 5;
%Tf = 1;
dt = 0.002;
Nt = round(Tf/dt)+1;
A = zeros(Nx,Nx);
B = zeros(Nx,Nu);
pert = 0.0001;
nx = 0;tx = 1;
ny = 1;ty = 0;
%s_frame = 20;      
s_frame = 1; 
e_frame = 76; 


dynamics_midpoint = @(x,u,dt) x + fun_xdot((x + fun_xdot(x,u,dt)*dt/2),u,dt)*dt;

vars = fun_data(); 
                                             
% initial conditions

BC      =  readmatrix("BC.xlsx");
omg     =  readmatrix("omg.xlsx");

%%%% extract the heel data

tht1=BC(3,s_frame);tht2=BC(4,s_frame);tht3=BC(5,s_frame);tht4=BC(6,s_frame);tht5=BC(7,s_frame);tht6=BC(8,s_frame);tht7=BC(9,s_frame);hx=BC(1,s_frame) ;hy=BC(2,s_frame);
omg1 = omg(3,s_frame); omg2 =omg(4,s_frame); omg3 = omg(5,s_frame);  omg5 = omg(7,s_frame); omg6 = omg(8,s_frame); omg7 = omg(9,s_frame);vhx =  omg(1,s_frame); vhy = omg(2,s_frame);
omg4 = omg(6,s_frame);
%x0 = [hx;hy;tht1;tht2;tht3;tht4;tht5;tht6;tht7;vhx;vhy;omg1;omg2;omg3;omg4;omg5;omg6;omg7];
x0 = [hx;hy;tht5;tht6;vhx;vhy;omg5;omg6];
%x0 = [0;0;pi/3;pi/4;1;1;2;3];

% goal
thtf1=BC(3,e_frame);thtf2=BC(4,e_frame);thtf3=BC(5,e_frame);thtf5=BC(7,e_frame);thtf6=BC(8,e_frame);thtf7=BC(9,e_frame);hfx=BC(1,e_frame);hfy=BC(2,e_frame);
omgf1 = omg(3,e_frame);omgf2 = omg(4,e_frame); omgf3 =  omg(5,e_frame); omgf5 = omg(7,e_frame);
omgf6 =   omg(8,e_frame);  omgf7 =  omg(9,e_frame); vhfx=  omg(1,e_frame); vhfy =  omg(2,e_frame);
omgf4 = omg(6,e_frame); 
thtf4=BC(6,e_frame);
%xf = [hfx;hfy;thtf1;thtf2;thtf3;thtf4;thtf5;thtf6;thtf7;vhfx;vhfy;omgf1;omgf2;omgf3;omgf4;omgf5;omgf6;omgf7];
xf = [hfx;hfy;thtf5;thtf6;vhfx;vhfy;omgf5;omgf6];



% costs
%Q = 1e-5*eye(4);
Q =  1e-5*eye(Nx);
Qf = 25*eye(Nx);
R = 5*1e-7*eye(Nu);
I = eye(Nu);
e_dJ = 1e-12;


% initialization
u = rand(Nu,Nt-1)*20;
u = ones(Nu,Nt-1)*20;
%u = zeros(Nu,Nt-1);
uf  = ones(Nu,Nt-1);
x = zeros(Nx,Nt);
x_prev = zeros(Nx,Nt);
x(:,1) = x0;



handle_f = @fun_xdot;
% first roll-out
for k = 2:Nt-1
           x(:,k) = fun_explicitEuler(handle_f,x(:,k-1),dt,u(:,k-1));
        k;
        %pause()           
        
end




val = x;
tdt = 0.2:0.01:0.74;
tt = 0.23:dt:1.58;
fr = 1:1:125;
cr = s_frame:1:76;
hr = s_frame:1:75;



%%% theta2
figure;
plot(cr,val(4,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(8, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta2 (N) \rightarrow');
legend('calc','dataset');


%%% theta1
figure;
plot(cr,val(3,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(7, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta1 (N) \rightarrow');
legend('calc','dataset');




%%% hx
figure;
plot(cr,val(1,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(1,s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('hx (N) \rightarrow');
legend('calc','dataset');


%%% hy
figure;
plot(cr,val(2,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(2, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('hy (N) \rightarrow');
legend('calc','dataset');

%


 
% original cost
J = 0;                                              
for k = 1:Nt-1
    J = J +  (x(:,k)-xf)'*Q*( x(:,k)-xf) + (u(:,k))'*R*(u(:,k)) 
end 
disp('Original cost:')
J = 0.5*(J + (x(:,Nt)-xf)'*Qf*(x(:,Nt)-xf)) 
val = x;



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
           r = R*u(:,k);           % lu
            
           
           % select alpha value for finding  A and B for each config at K

           A = fun_amat(x(:,k),u(:,k),dt)
            B = fun_bmat(x(:,k),u(:,k),dt)

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
             dJ = dJ +  gu'*d(:,k)
 disp('ITERATIPONSSSSSSSSSSSSSSSSSSSSS')
       k
       iter
    
     % pause()
    end
    disp('ITERATIPONSSSSSSSSSSSSSSSSSSSSS- cOMPLETED')
         % pause()
    
  
%%%% End of Backward Pass %%%%%
   %Forward rollout with line search
    xn(:,1) = x(:,1);
    alpha = 1.0;
    indic = 15;
    xdata = x0;   
       
   for kc = 1:(Nt-1)
        un(:,kc) = u(:,kc) - alpha*d(:,kc) - (K(:,:,kc)*(xn(:,kc)-x(:,kc)));
        %xn(:,kc+1) = dynamics_midpoint(xn(:,kc),un(:,kc),dt);
          xn(:,kc+1) = fun_explicitEuler(handle_f,xn(:,kc),dt,un(:,kc));
        
       disp('ITERATIPONSSSSSSSSSSSSSSSSSSSSS')
       kc
     %pause()
    end
    disp('EOFP')
     %pause() 
    Jn = 0;
    for k = 1:Nt-1
        Jn = Jn + (xn(:,k)-xf)'*Q*(xn(:,k)-xf) + un(:,k)'*R*un(:,k);
    end
   Jn = 0.5*(Jn + (xn(:,Nt)-xf)'*Qf*(xn(:,Nt)-xf))
    
     disp('line search')
     %pause() 
 liter =1;
    while isnan(Jn) || Jn > (J - 1e-2*alpha*dJ)
         disp('Inside search')
        alpha = 0.5*alpha
          
        for k = 1:(Nt-1)
            un(:,k) = u(:,k) - alpha*d(:,k) - (K(:,:,k)*(xn(:,k)-x(:,k)));
            %xn(:,k+1) = dynamics_midpoint(xn(:,k),un(:,k),dt);
            xn(:,k) = fun_explicitEuler(handle_f,xn(:,k),dt,un(:,k));
        end
        %Jn = cost(xn,un);
        Jn = 0;
        for k = 1:Nt-1
            Jn = Jn + (xn(:,k) - xf)'*Q*(xn(:,k) - xf) + un(:,k)'*R*un(:,k);
        end
     Jn = 0.5*(Jn + (xn(:,Nt) - xf)'*Qf*(xn(:,Nt) - xf))
     % pause()
    end
     %  pause() 
 
   % J = Jn;
    x = xn;
    u = un;
   %if iter > 5
    %   break
    %end
  end

val = x;
%}