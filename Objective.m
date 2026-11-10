function z=Objective1(x)
global teta L E v n Kc M Sigmay y Kb Ks z
format long
%%
x1=x(1);
x2=(.29/19)*x1;
x3=(.9/19)*x1;
%% Objective  Function

%% g1
g01=((24*(1-(v^2)))/(Kc*(pi^2)*(n^2)))*(cot(teta))^2;
g02=(v^2)/(E*M);
g03=(x1-x3)/(x3^2);
g1=.5*x3+((n*x2*sec(teta))/12);
g1=(1/g1)*g01*g02*g03-1;
%% g2
g21=(E/Sigmay)*sqrt((n^2)+3*((cot(teta))^2)*(x1-x3));
g22=(n*x3+((x2*n*sec(teta))/6))*(x1-x3);
g2=((v^2)/(E*M))*(g21/g22)-1;
%% g3
g31=(E/Sigmay)*(x3+((n*sec(teta))/6)*x2)*(x1-x3);
g32=sqrt(4*((y/L)/(x1-x3))^2+3*((cot(teta)/n)*(x1-x3)*(x3/x2)+csc(teta)*(x1-x3)*(.25-((y/L)/(x1-x3))^2))^2);
g3=g31/g32;
g3=1/g3;
g3=((v^2)/(E*M))*g3-1;
%% g4
g41=(((n*pi*sin(teta))^2)/(12*(1-v^2)))*(x2^2)*(x3+((n*x2*sec(teta))/6));
g42=(x1-x3)*sqrt(((1/Kb)^2)+((1/Ks)^2)*((x1-x3)^2)*(((cot(teta)*x3)/(n*x2))+.25*(csc(teta))^2));
g4=g41/g42;
g4=1/g4;
g4=(((v^2)/E*M)*g4)-1;
if g1<=0 & g2<=0 & g3<=0 & g4<=0
Sai=2*x3+sec(teta)*n*x2;
Pi=.0035*E*(L^3)*x1*x2*x3;  
%%
z=[Sai;Pi];
end
end