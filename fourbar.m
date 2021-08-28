clear all;
clc;
close all;

slc_data1; % This is an external file that is being included here. 

 
% create an array of theta2. We will be solving for theta3 and r1 for
% each of these 
thetas2 = 0.0:15*pi/180:2*pi;
% % velocity and  acceleration

positions = [];
velocities = [];
accelerations = [];


for I = 1: length(thetas2)
    [pos,vels,accs] = slidercrank_soln(r1,r2,r3,r4,thetas2(I), theta2dot, theta2ddot, init_values);
    positions = [positions pos];
    velocities = [velocities vels];
    accelerations = [accelerations accs];
end

thetas3 = positions(1,:);
thetas4=positions(2,:);
R1 = r1.*(ones(1,25));
R2 = r2*exp(i*thetas2);
R3 = r3*exp(i*thetas3);
R4 = r4*exp(i*thetas4);

    %%the magnitudes of vectors r1, r2, r3, and r4 are being
    %%sorted, and grashof's condition of s+p+q>l is being checked here
Order=[r1 r2 r3 r4];
MS=sort(Order);
L=MS(1,4);              
S=MS(1,1);
P=MS(1,2);
Q=MS(1,3);

if S+P+Q>L
    fprintf('This is a Grashof linkage\n');
end

figure;
cp = R1 + R3*(2/3); %intermediate point on link 3

%%setting up plot for mechanism
   X(1) = 0; Y(1) = 0;
   X(2) = real(R2(1)); Y(2) = imag(R2(1));
   X(3) = real(R2(1)+R3(1)); Y(3) = imag(R2(1)+R3(1));
   X(4) = real(R2(1)+R3(1)-R4(1)); Y(4)=imag(R2(1)+R3(1)-R4(1));
   X(5)=0;Y(5)=0;      
   X=X';Y=Y';
   plot(X,Y,'ro');      
   line(X,Y);               %%plots the first frame of mechanism
   axis([-5 16 -5 16]);     %%to allow ginput from user
   
%%ginput cases
%%user selects point, and this point is projected onto a nearest link
   cog=0;
   R=[];
   [a,b]=ginput(1);  
   coord=[a,b];
   A2=-(Y(2)-0)/(X(2)-0);
   B2=1;
   C2=0;
   A3=-(Y(3)-Y(2))/(X(3)-X(2));
   B3=1;
   C3=-Y(2)+X(2)*((Y(3)-Y(2))/(X(3)-X(2)));
   A4=-(Y(4)-Y(3))/(X(4)-X(3));
   B4=1;
   C4=-Y(3)+X(3)*((Y(4)-Y(3))/(X(4)-X(3)));
   
   if A2>10^3
       d2=abs(x_input-X(2));
   else
       d2=abs((A2*a+B2*b+C2)/sqrt(A2^2+B2^2));
   end
   
   if A3>10^3
       d3=abs(a-X(3));
   else
       d3=abs((A3*a+B3*b+C3)/sqrt(A3^2+B3^2));
   end
   if A4>10^3
       d4=abs(a-X(4));
   else
       d4=abs((A4*a+B4*b+C4)/sqrt(A4^2+B4^2));
   end
   
   xperp2=(B2*(B2*a-A2*b)-A2*C2)/(A2^2+B2^2); yperp2=(A2*(-B2*a+A2*b)-B2*C2)/(A2^2+B2^2);
   xperp3=(B3*(B3*a-A3*b)-A3*C3)/(A3^2+B3^2); yperp3=(A3*(-B3*a+A3*b)-B3*C3)/(A3^2+B3^2);
   xperp4=(B4*(B4*a-A4*b)-A4*C4)/(A4^2+B4^2); yperp4=(A4*(-B4*a+A4*b)-B4*C4)/(A4^2+B4^2);

   
   if  d2<d4 && d2<d3                
       cog=sqrt(xperp2^2+yperp2^2); %%projects onto link 2
       R=cog*exp(i*thetas2);
       elseif d3<d2 && d3<d4          
       cog=sqrt((xperp3-X(2))^2+(yperp3-Y(2))^2); %%projects onto link 3
       R=R2+cog*exp(i*thetas3);
       elseif d4<d2 && d4<d3                      
       cog=sqrt((xperp4-X(4))^2+(yperp4-Y(4))^2); %%projects onto link 4
       R=R1+cog*exp(i*thetas4);
   end
      
%%plotting loop   
for k = 1:length(thetas2)
   X(1) = 0; Y(1) = 0;
   X(2) = real(R2(k)); Y(2) = imag(R2(k));
   X(3) = real(R2(k)+R3(k)); Y(3) = imag(R2(k)+R3(k));
   X(4) = real(R2(k)+R3(k)-R4(k)); Y(4) = imag(R2(k)+R3(k)-R4(k));
   X(5)=0;Y(5)=0;
   X = X';
   Y = Y';
  clf;
  hold on;
  axis([-5 16 -5 16]);
  plot(R(1:k),'o');
  line(X,Y);
  plot(X,Y,'ro');
  title('Four Bar Mechanism')
  pause(0.01);
  
end


%%force analysis
Fstatic = [];
Fdynamic= [];
Rg2s = [];
Rg3s = [];
Rg4s=[];
ag2s = [];
ag3s = [];
ag4s=[];
thetas2dot = theta2dot.*(ones(1,length(thetas2)));
thetas3dot = velocities(1,:);
thetas4dot = velocities(2,:);
thetas2ddot = theta2ddot*(ones(1,length(thetas2)));
thetas3ddot= accelerations(1,:);
thetas4ddot=accelerations(2,:);
% Ans=[thetas3(2) thetas4(2) thetas3dot(2) thetas4dot(2) thetas3ddot(2) thetas4ddot(2)] %%answers for data set: 
V1 = 0+0i;
V2 = r2.*i.*exp(i.*thetas2).*thetas2dot;
V3 = r3.*i.*exp(i.*thetas3).*thetas3dot;
V4 = r4.*i.*exp(i.*thetas4).*thetas4dot;
T12=300; T34=-100;

for K= 1: length(thetas2)
   theta2=thetas2(K); theta3=thetas3(K); theta4=thetas4(K);
   theta2dot=thetas2dot(K); theta2ddot=thetas2ddot(K);
   theta3dot=thetas3dot(K); theta3ddot=thetas3ddot(K);
   theta4dot=thetas4dot(K); theta4ddot=thetas4ddot(K);

   Rg2=b2.*exp(i*(theta2+phi2));
   Rg3=R2(K)+b3*exp(i*(theta3+phi3)); 
   Rg4=R1(K)+b4*exp(i*(theta4+phi4));
   
   R12g2=0-Rg2; R12g2x = real(R12g2); R12g2y = imag(R12g2);
   R32g2 = R2(K)-Rg2; R32g2x = real(R32g2); R32gy = imag(R32g2);
   R23g3 = R2(K)-Rg3; R23g3x = real(R23g3); R23g3y = imag(R23g3);
   R43g3 = R1(K)-Rg3; R43g3x = real(R43g3); R43g3y = imag(R43g3);
   R34g4 = R3(K)-Rg4; R34g4x = real(R34g4); R34g4y = imag(R34g4);
   R14g4 = R4(K)-Rg4; R14g4x = real(R14g4); R14g4y = imag(R14g4);
    
%    a4=0*i;
   ag2= -r2*exp(i*theta2).*(theta2dot.^2)+i*r2*exp(i*theta2).*theta2ddot;
   ag3= -r3*exp(i*theta3).*(theta3dot.^2)+i*r3*exp(i*theta3).*theta3ddot;
   ag4= -r4*exp(i*theta4).*(theta4dot.^2)+i*r4*exp(i*theta4).*theta4ddot;

    Rg2s(K)=Rg2;
    Rg3s(K)=Rg3;
    Rg4s(K)=Rg4;
    ag2s(K)=ag2;
    ag3s(K)=ag3;
    ag4s(K)=ag4;
    
    ag2x=real(ag2);ag2y=imag(ag2);
    ag3x=real(ag3);ag3y=imag(ag3);
    ag4x=real(ag4);ag4y=imag(ag4);
    
    F=zeros(13,13);
    F(1,1)=1;F(1,4)=1;
    F(2,2)=1;F(2,5)=1;
    F(3,1)=-R12g2x;F(3,2)=R12g2y;F(3,5)=-1;F(3,3)=-R32g2x; F(3,5)=-R32g2x;
    F(4,6)=1;F(4,8)=1;
    F(5,7)=1;F(5,9)=1;
    F(6,6)=-R23g3x; F(6,7)=R23g3y; F(6,8)=-R43g3x;F(6,9)=R43g3y;
    F(7,10)=1;F(7,12)=1;
    F(8,11)=1;F(8,13)=1;
    F(9,10)=R14g4x; F(9,11)=-R14g4y; F(9,12)=R14g4x; F(9,13)=-R14g4y;
    F(10,3)=1; F(10,6) = 1;
    F(11,4)=1; F(11,7) = 1;
    F(12,8) = 1; F(12,10)=1;
    F(13,9) = 1; F(13,11)=1;

    RHS=[0 0 0 0 0 0 0 0 T12 0 0 0 0 ]';
    RHS2=[m2*ag2x; m2*ag2y; Ig2*theta2ddot; m3*ag3x; m3*ag3y; Ig3*theta3ddot; m4*ag4x; m4*ag4y; Ig4*theta4ddot-T34; 0; 0;0;0];

fstatic=F\RHS;
fdynamic=F\RHS2;
shaking_f_x = fstatic(2)+fdynamic(2)+fstatic(12)+fdynamic(12);
shaking_f_y = fstatic(3)+fdynamic(3)+fstatic(13)+fdynamic(13);
shaking_f=shaking_f_x+j.*shaking_f_y;

Fstatic=[Fstatic fstatic];
Fdynamic=[Fdynamic fdynamic];
shakingforce(K)=shaking_f;
end
figure;
title('Static Force Analysis');
polarplot(thetas2, Fstatic(3,:));
figure;
title('Dynamic Force Analysis');
polarplot(thetas2, Fdynamic(3,:));


