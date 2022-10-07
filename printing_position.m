close all;
clc;
clearvars -except val;

frame  = 20; 
% frame 20 taken for start of SSP1
% frame 55 midswing phase
% frame 75 is taken as end of ssp
%

BC      =  readmatrix("BC.xlsx"); 
omg     =  readmatrix("omg.xlsx");


Lf = 0.1968;
L4  = 0.0841;
L1=  0.4041;
L2 = 0.4418;
L5 = 0.4418;
L3 = 0.4033;
L6 = 0.4033;
r4 = 0.0432;
r7 = 0.0442;
gamma71 =  1.2147;
gamma72 = 3.9843; 
gamma61 = 1.2104;
gamma62 = 3.9711;
r7t = 0.1115;
r7h = 0.0887;
r4t = 0.1114 ;
r4h = 0.0877 ;
r4c =   0.0873;%    0.0992;
gamma43 =  3.9381 %3.6666;
gamma2 = 0.7110;

%i = 20;
i = 1;
%while i < 76
while i <  size(val)
  
  %tht1=val(1,i);tht2=val(2,i);tht3=val(3,i);tht4=val(4,i);tht5=val(5,i);tht6=val(6,i);tht7=val(7,i);hx=val(8,i);hy=val(9,i);
 %tht1=BC(3,i);tht2=BC(4,i);tht3=BC(5,i);tht4=BC(6,i);tht5=BC(7,i);tht6=BC(8,i);tht7=BC(9,i);hx=BC(1,i);hy=BC(2,i); copx=BC(11,i);
 %copy=BC(10,i)
%  x0 = [tht1;tht2;tht3;tht4;tht5;tht6;tht7;hx;hy]
%tht3=val(1,i);tht4=val(2,i);hx=val(3,i);hy=val(4,i);
%hx=BC(10,frame+i);
%hy=BC(11,frame+i); 
tht1=val(3,i);tht2=val(4,i);
hx=val(1,i);hy=val(2,i);  




xk = hx  + L1*cos(tht1+pi); 
yk = hy  + L1*sin(tht1+pi); 
xa = hx + L1*cos(tht1+pi) + L2*cos(tht2+pi);
ya = hy + L1*sin(tht1+pi) + L2*sin(tht2+pi);







%xcop = hx + L1*cos(pi + tht1) + L2*cos(pi + tht2) + L3*cos(pi + tht3) + r4*cos(pi + tht4) + r4c*cos(gamma43 + tht4);
%ycop = hy + L1*sin(pi + tht1) + L2*sin(pi + tht2) + L3*sin(pi + tht3) + r4*sin(pi + tht4) +  r4c*sin(gamma43 + tht4);





axis([-0.5 2 -2 0.5]) 
%axis([-0.5 1.5 -0.5 0.5]) 
base =line([-1 2],[0.02 0.02],'LineWidth',1,'Color','black');
%pelvic =line([lhip_x(i) rhip_x(i)],[lhip_z(i) rhip_z(i)],'LineWidth',1,'Color','black');
  T=line([hx  xk],[hy  yk],'LineWidth',1,'Color','black');
  u=line([xk xa],[yk ya],'LineWidth',1,'Color','red');

  

 i = i +1;
 %pause
end  



