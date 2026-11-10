function z=Objective2(x)
global teta L E v n Kc M Sigmay z x1 x2 x3 Sai Pi
format long
%%
x1=x(1);
x2=((.29/19)*x(2));
x3=((.9/19)*x(3));
%% Kf Kc
Kf1=2.4*cos(teta)*((x2/x3)^3)+1;
Kf2=1.2*cos(teta)*((x2/x3)^3)+1;
Kf=(Kf1/Kf2)^2;
% Kf=10;
Kc1=2.2+1.2*sec(teta)*((x3/x2)^3);
Kc2=1.6+.6*sec(teta)*((x3/x2)^3);
Kc=1.375*(Kc1/Kc2);
% Kc=10;
%% g1
g11=(E/Sigmay)*(1/x3);
g12=x1-x3;
g13=g11/g12;
g1=((v^2)/E*M)*g13-1;
%% g2
g2=((v^2)/E*M)*(E/Sigmay)*(1/x2)*((E*csc(teta))/n)-1;
%% g3
g31=48*(x1-x3);
g32=Kf*pi^2*n^2*tan(teta)*(x3^3);
g3=((v^2)/E*M)*(g31/g32);
%% g4
g41=12*(x1-x3)^2;
g42=Kc*(pi^2)*(n^3)*((sin(teta))^3);
g4=((v^2)/E*M)*(g41/g42)-1;
%%
if g1<=0 && g2<=0 && g3<=0 && g4<=0
Sai=-min(2*x3+sec(teta)*n*x2);
Pi=max(.0035*E*(L^3)*x1*x2*x3);  
%%
z=[Sai;Pi];
end
end