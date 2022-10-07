%%%%%%%%%%%%% GAIT SIMULATION %%%%%%%%%%%%%%%%%%%%%
% Red colour - left leg
% Blue  color - right leg 
% Gait cycle time =2.1 seconds
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clearvars;
close all
cord=readmatrix("joint_coordinates.xlsx");
hip=readmatrix("lhip.xlsx");
i=1;
iteration=0;
counter=1;
  txt = '-';
x1 = cord(:,14)
xiphoid_x=cord(:,14);
y1 =cord(:,15);
xiphoid_z=cord(:,15);

time  = cord(:,1);
lhip_x=cord(:,2);
lhip_z=cord(:,3);
lknee_x=cord(:,6);
lknee_z=cord(:,7);
lankle_x=cord(:,10);
lankle_z=cord(:,11);
lmeta_x=cord(:,22);
lmeta_z=cord(:,23);
lpost_x=cord(:,18);
lpost_z=cord(:,19);
lfeetmid_x=(lmeta_x+lpost_x)./2;
lfeetmid_z=(lmeta_z+lpost_z)./2;

rhip_x=cord(:,4);
rhip_z=cord(:,5);
rknee_x=cord(:,8);
rknee_z=cord(:,9);
rankle_x=cord(:,12); 
rankle_z=cord(:,13);
rmeta_x=cord(:,20);
rmeta_z=cord(:,21);
rpost_x=cord(:,16);
rpost_z=cord(:,17);
rfeetmid_x=(rmeta_x+rpost_x)./2;
rfeetmid_z=(rmeta_z+rpost_z)./2

hip_x=(rhip_x+rhip_x)./2;
hip_z=(rhip_z+rhip_z)./2;
x0 = zeros(9,8);
xv  = zeros(9,1);
vel  = zeros(9,8);

%figure; hold on; grid on;
while iteration <counter
    i=1;
    while i<214
        base =line([-1000 2000],[20 20],'LineWidth',1,'Color','black');
        pelvic =line([lhip_x(i) rhip_x(i)],[lhip_z(i) rhip_z(i)],'LineWidth',1,'Color','black');
        T=line([xiphoid_x(i) ((lhip_x(i)+rhip_x(i))/2)],[xiphoid_z(i) ((lhip_z(i)+rhip_z(i))/2)],'LineWidth',1,'Color','black');
     
         u=line([lhip_x(i) lknee_x(i)],[lhip_z(i) lknee_z(i)],'LineWidth',1,'Color','red');
        v=line([lknee_x(i) lankle_x(i)],[lknee_z(i) lankle_z(i)],'LineWidth',1,'Color','red');
        w=line([lankle_x(i)  lmeta_x(i)],[lankle_z(i)  lmeta_z(i)],'LineWidth',1,'Color','red');
        x=line([lankle_x(i) lpost_x(i)],[lankle_z(i) lpost_z(i)],'LineWidth',1,'Color','red');
        y=line([lpost_x(i)  lmeta_x(i)],[lpost_z(i)  lmeta_z(i)],'LineWidth',1,'Color','red');
       % legend('Left leg')

        u1=line([rhip_x(i) rknee_x(i)],[rhip_z(i) rknee_z(i)],'LineWidth',1,'Color','blue');
        v1=line([rknee_x(i) rankle_x(i)],[rknee_z(i) rankle_z(i)],'LineWidth',1,'Color','blue'); 
        w1=line([rankle_x(i)  rmeta_x(i)],[rankle_z(i)  rmeta_z(i)],'LineWidth',1,'Color','blue');
        x1=line([rankle_x(i) rpost_x(i)],[rankle_z(i) rpost_z(i)],'LineWidth',1,'Color','blue');
        y1=line([rpost_x(i)  rmeta_x(i)],[rpost_z(i)  rmeta_z(i)],'LineWidth',1,'Color','blue');
        %legend('','Torso','Left leg','','','','','Right leg','','','','')
        t1= text(1000,1500,txt,'Color','red','FontSize',14);
     





    theta_HAT(i) = ( atan2(xiphoid_z(i) - hip_z(i),xiphoid_x(i)- hip_x(i)));
    theta_l_thigh(i) =(atan2(hip_z(i)-cord(i,7),hip_x(i) - cord(i,6)));
    theta_l_leg(i) = ( atan2(cord(i,7)-cord(i,11),cord(i,6)-cord(i,10)));
    theta_r_thigh(i) = (atan2(hip_z(i)-cord(i,9),hip_x(i)- cord(i,8)));
    theta_r_leg(i) = (atan2(cord(i,9) - cord(i,13),cord(i,8) - cord(i,12)));
    theta_l_ankle(i) = (  atan2(cord(i,11)-lpost_z(i),cord(i,10)-lpost_x(i)));
    theta_r_ankle(i) = (   atan2(cord(i,13)-rpost_z(i),cord(i,12)-rpost_x(i)));


      
        if cord(i,1) == 0
           x0(:,1) = [ theta_HAT(i); theta_l_thigh(i);theta_l_leg(i);theta_l_ankle(i);theta_r_thigh(i);theta_r_leg(i);theta_r_ankle(i);xiphoid_x(i)/1000;xiphoid_z(i)/1000];
           %  xv(:,1) = [ theta_HAT(i-1); theta_l_thigh(i-1);theta_l_leg(i-1);theta_l_ankle(i-1);theta_r_thigh(i-1);theta_r_leg(i-1);theta_r_ankle(i-1);xiphoid_x(i-1)/1000;xiphoid_z(i-1)/1000] ;
            % vel(:,1) = (x0(:,4) -  xv(:,1))/0.01;  

        elseif cord(i,1) == 0.05
           x0(:,2) = [ theta_HAT(i); theta_l_thigh(i);theta_l_leg(i);theta_l_ankle(i);theta_r_thigh(i);theta_r_leg(i);theta_r_ankle(i);xiphoid_x(i)/1000;xiphoid_z(i)/1000] ;  
             xv(:,1) = [ theta_HAT(i-1); theta_l_thigh(i-1);theta_l_leg(i-1);theta_l_ankle(i-1);theta_r_thigh(i-1);theta_r_leg(i-1);theta_r_ankle(i-1);xiphoid_x(i-1)/1000;xiphoid_z(i-1)/1000] ;
             vel(:,2) = (x0(:,2) -  xv(:,1))/0.01;  
        elseif cord(i,1)==0.25
        x0(:,3) = [ theta_HAT(i); theta_l_thigh(i);theta_l_leg(i);theta_l_ankle(i);theta_r_thigh(i);theta_r_leg(i);theta_r_ankle(i);xiphoid_x(i)/1000;xiphoid_z(i)/1000] ;
             xv(:,1) = [ theta_HAT(i-1); theta_l_thigh(i-1);theta_l_leg(i-1);theta_l_ankle(i-1);theta_r_thigh(i-1);theta_r_leg(i-1);theta_r_ankle(i-1);xiphoid_x(i-1)/1000;xiphoid_z(i-1)/1000] ;
             vel(:,3) = (x0(:,3) -  xv(:,1))/0.01;  
        elseif cord(i,1)==0.54
             x0(:,4) = [ theta_HAT(i); theta_l_thigh(i);theta_l_leg(i);theta_l_ankle(i);theta_r_thigh(i);theta_r_leg(i);theta_r_ankle(i);xiphoid_x(i)/1000;xiphoid_z(i)/1000] ;
             xv(:,1) = [ theta_HAT(i-1); theta_l_thigh(i-1);theta_l_leg(i-1);theta_l_ankle(i-1);theta_r_thigh(i-1);theta_r_leg(i-1);theta_r_ankle(i-1);xiphoid_x(i-1)/1000;xiphoid_z(i-1)/1000] ;
             vel(:,4) = (x0(:,4) -  xv(:,1))/0.01;
        
        elseif cord(i,1)==0.75
             x0(:,5) = [ theta_HAT(i); theta_l_thigh(i);theta_l_leg(i);theta_l_ankle(i);theta_r_thigh(i);theta_r_leg(i);theta_r_ankle(i);xiphoid_x(i)/1000;xiphoid_z(i)/1000] ;
             xv(:,1) = [ theta_HAT(i-1); theta_l_thigh(i-1);theta_l_leg(i-1);theta_l_ankle(i-1);theta_r_thigh(i-1);theta_r_leg(i-1);theta_r_ankle(i-1);xiphoid_x(i-1)/1000;xiphoid_z(i-1)/1000] ;
             vel(:,5) = (x0(:,5) -  xv(:,1))/0.01;

        elseif cord(i,1)==1
             x0(:,6) = [ theta_HAT(i); theta_l_thigh(i);theta_l_leg(i);theta_l_ankle(i);theta_r_thigh(i);theta_r_leg(i);theta_r_ankle(i);xiphoid_x(i)/1000;xiphoid_z(i)/1000] ;    
             xv(:,1) = [ theta_HAT(i-1); theta_l_thigh(i-1);theta_l_leg(i-1);theta_l_ankle(i-1);theta_r_thigh(i-1);theta_r_leg(i-1);theta_r_ankle(i-1);xiphoid_x(i-1)/1000;xiphoid_z(i-1)/1000] ;
             vel(:,6) = (x0(:,6) -  xv(:,1))/0.01 
       
        elseif cord(i,1)==1.24
             x0(:,7) = [ theta_HAT(i); theta_l_thigh(i);theta_l_leg(i);theta_l_ankle(i);theta_r_thigh(i);theta_r_leg(i);theta_r_ankle(i);xiphoid_x(i)/1000;xiphoid_z(i)/1000] ;
             xv(:,1) = [ theta_HAT(i-1); theta_l_thigh(i-1);theta_l_leg(i-1);theta_l_ankle(i-1);theta_r_thigh(i-1);theta_r_leg(i-1);theta_r_ankle(i-1);xiphoid_x(i-1)/1000;xiphoid_z(i-1)/1000] ;
             vel(:,7) = (x0(:,7) -  xv(:,1))/0.01;          
        elseif cord(i,1)==1.39
             x0(:,8) = [ theta_HAT(i); theta_l_thigh(i);theta_l_leg(i);theta_l_ankle(i);theta_r_thigh(i);theta_r_leg(i);theta_r_ankle(i);xiphoid_x(i)/1000;xiphoid_z(i)/1000] ;
             xv(:,1) = [ theta_HAT(i-1); theta_l_thigh(i-1);theta_l_leg(i-1);theta_l_ankle(i-1);theta_r_thigh(i-1);theta_r_leg(i-1);theta_r_ankle(i-1);xiphoid_x(i-1)/1000;xiphoid_z(i-1)/1000] ;
             vel(:,8) = (x0(:,8) -  xv(:,1))/0.01;      
        end   
          
        i=i+1;
        pause()
    end
 iteration=iteration+1;
 pause(0.002);
end

time = cord(:,1);
gtc = 1:125; % gait cycle time is from 0 to 1.24 seconds
tfit = time(gtc);
dt = 0.01; % dt should be set to 0.01 only
tt = 0:dt:1.24;


f1 = fit(tfit,theta_HAT(1:125)','fourier8');
tht1 = f1(tt);
f2 = fit(tfit,theta_l_thigh(1:125)','fourier8');
tht2 = f2(tt);
f3 = fit(tfit,theta_l_leg(1:125)','fourier8');
tht3 = f3(tt);
f4 = fit(tfit,theta_l_ankle(1:125)','fourier8');
tht4 = f4(tt);
f5 = fit(tfit,theta_r_thigh(1:125)','fourier8');
tht5 = f5(tt);
f6 = fit(tfit,theta_r_leg(1:125)','fourier8');
tht6 = f6(tt);
f7 = fit(tfit,theta_r_ankle(1:125)','fourier8');
tht7 = f7(tt);


fx1 = fit(tfit,x1,'fourier8');
x1 = fx1(tt);
fy1 = fit(tfit,y1,'fourier8');
y1 = fy1(tt);

%%% accelerations
[vx1,ax1] = differentiate(fx1,tt);
[vy1,ay1] = differentiate(fy1,tt);
[omg1,alp1] = differentiate(f1,tt);
[omg2,alp2] = differentiate(f2,tt);
[omg3,alp3] = differentiate(f3,tt);
[omg4,alp4] = differentiate(f4,tt);
[omg5,alp5] = differentiate(f5,tt);
[omg6,alp6] = differentiate(f6,tt);
[omg7,alp7] = differentiate(f7,tt);

omg = [omg1';omg2';omg3';omg4';omg5';omg6';omg7'];


%{
% obtaining GRF and GRT data for the four stages of the gait cycle
GRdata = readmatrix("GRdata.xlsx");
%t1 = time(24:76);
gf1 = GRdata(24:76,:);
%t2 = time(77:91);
gf2 = GRdata(77:91,:);
%t3 = time(92:144);
gf3 = GRdata(92:144,:);
%t4 = time(145:159);
f41 = GRdata(77:91,2:7);
f42 = GRdata(77:91,8:13);
f4 = [t4 f42 f41];
time_plot = [t1;t2;t3;t4];
GRdata_new = [f1;f2;f3;f4];


%}














%{


lhip_x=time(:,2);
lhip_z=time(:,3);
lknee_x=time(:,6);
lknee_z=time(:,7);
lankle_x=time(:,10);
lankle_z=time(:,11);

rhip_x=time(:,4);
rhip_z=time(:,5);
rknee_x=time(:,8);
rknee_z=time(:,9);
rankle_x=time(:,12);
rankle_z=time(:,13);


i=1;
while i<214
    theta_l_thigh(i)=atan2d(cord(i,7)-cord(i,3),cord(i,6)-cord(i,2) );
    theta_l_leg(i)=atan2d(cord(i,11)-cord(i,7),cord(i,10)-cord(i,6) );
    theta_r_thigh(i)=atan2d(cord(i,9)-cord(i,5),cord(i,8)-cord(i,4) );
    theta_r_leg(i)=atan2d(cord(i,13)-cord(i,9),cord(i,12)-cord(i,8) );
    theta_l_ankle(i)=atan2d(cord(i,11)-lfeetmid_z(i),cord(i,10)-lfeetmid_x(i) );
    theta_r_ankle(i)=atan2d(cord(i,13)-rfeetmid_z(i),cord(i,12)-rfeetmid_x(i) );
    i=i+1;
end

%}


%(lhip_z(i,1)-lknee_z(i,1)),(lhip_x(i,1)-lknee_x(i,1))
% ),(time(2,i)-time(6,i)))

%theta = linspace(0,10,180);
%y = exp(x/10).*sin(4*x);
figure; hold on; grid on;
plot(cord(1:213,1),theta_l_thigh,'-o');
title('Theta left thigh ');
ylabel('Angle in degrees ');
xlabel('time frame(s)');

figure; hold on; grid on;
plot(cord(1:213,1),theta_l_leg,'-o');
title('Theta left leg');
ylabel('Angle in degrees ');
xlabel('time frame(s)');

figure; hold on; grid on;
plot(cord(1:213,1),theta_r_thigh,'-o');
title('Theta right thigh');
ylabel('Angle in degrees ');
xlabel('time frame(s)');
%{
figure; hold on; grid on;
plot(cord(1:213,1),theta_r_leg,'-o');
title('Theta right leg');
ylabel('Angle in degrees ');
xlabel('time frame(s)');
time=cord(1:213,1);
%p=polyfit(cord(1:213,1),lthigh,100)
% pp = spline(lthigh,cord(1:213,1) );
% [~, coeffs] = unmkpp(pp);


%}
writematrix(x0,'BV.xlsx');
writematrix(vel,'vel.xlsx');
writematrix(omg,'omg.xlsx');



