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
xiphoid_x=cord(:,14);
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
lfeetmid_x=(lmeta_x+lpost_x)./2
lfeetmid_z=(lmeta_z+lpost_z)./2

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
rfeetmid_z=(rmeta_z+rpost_z)./2;


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
      %delete(t1)
        
        if cord(i,1)==0.05
            delete(t1)
            txt = 'Left foot heel strike'
            t1= text(1000,1500,txt,'Color','red','FontSize',14);
            pause()
        elseif cord(i,1)==0.25
             delete(t1)
             txt = 'Right foot toe off'
             t1= text(1000,1500,txt,'Color','red','FontSize',14);
               pause()
            
        elseif cord(i,1)==0.75
            delete(t1)
            txt = 'Right foot heel strike'
            t1= text(1000,1500,txt,'Color','red','FontSize',14);
            pause()
            
        elseif cord(i,1)==1
            delete(t1)
            txt = 'Left foot toe off'
            t1= text(1000,1500,txt,'Color','red','FontSize',14);
            pause()
            
        elseif cord(i,1)==1.39
            delete(t1)
            txt = 'Left foot heel strike'
            t1= text(1000,1500,txt,'Color','red','FontSize',14);  
            pause()
        elseif cord(i,1)==1.66
            delete(t1)
            txt = 'Right foot toe off'
            t1= text(1000,1500,txt,'Color','red','FontSize',14);
            pause()
        elseif cord(i,1)==2.1
            delete(t1)
            txt = 'Left foot heel strike'
            t1= text(5000,1500,txt,'Color','red','FontSize',14);
            pause()
        elseif i==213
             pause(2)

        end   

        


         
        %delete(t1)
        pause(0.1);
        delete(t1)
        delete(u);
        delete(w);
        delete(x);
       delete(y);
        delete(v);
        delete(u1);
        delete(v1);
        delete(w1);
        delete(x1);
         delete(y1);
        delete(T);
        delete(pelvic);
        grid on; %if you want the grid to show up.
        axis('equal'); %make the axis equal, to avoid scaling effect
        xlim([-1000,2000]);
        ylim([-10,1800]);
        xlabel('x'); ylabel('y');
        txt1 = 'Left foot Red';
        txt2 = 'Right foot blue';
       % t2= text(1500,1500,txt1);
        i=i+1;

    end
 iteration=iteration+1;
 pause(0.002);
end

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
%}

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

figure; hold on; grid on;
plot(cord(1:213,1),theta_r_leg,'-o');
title('Theta right leg');
ylabel('Angle in degrees ');
xlabel('time frame(s)');
time=cord(1:213,1);
%p=polyfit(cord(1:213,1),lthigh,100)
% pp = spline(lthigh,cord(1:213,1) );
% [~, coeffs] = unmkpp(pp);








