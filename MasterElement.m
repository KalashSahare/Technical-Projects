function [Ke]=MasterElement(P1,D)
%% Physical to Master (cubic shape function)
% cubic polynomial a1+a2x+a3y+a4z+a5xy+a6yz+a7zx+a8xyz
% L=10;
% B1=10;
% H=10;
% P1(1,:)=[0,0,0];
% P1(2,:)=[L,0,0];
% P1(3,:)=[L,B1,0];
% P1(4,:)=[0,B1,0];
% P1(5,:)=[0,0,H];
% P1(6,:)=[L,0,H];
% P1(7,:)=[L,B1,H];
% P1(8,:)=[0,B1,H];

% E=200*10^9;%N/m2
% mu=0.3;
% D=(E/((1+mu)*(1-2*mu)))*[1-mu,mu,mu,0,0,0;mu,1-mu,mu,0,0,0;mu,mu,1-mu,0,0,0;0,0,0,(1-2*mu)/2,0,0;0,0,0,0,(1-2*mu)/2,0;0,0,0,0,0,(1-2*mu)/2];


%% Master Element
P(1,:)=[-1,-1,-1];
P(2,:)=[1,-1,-1];
P(3,:)=[1,1,-1];
P(4,:)=[-1,1,-1];
P(5,:)=[-1,-1,1];
P(6,:)=[1,-1,1];
P(7,:)=[1,1,1];
P(8,:)=[-1,1,1];

syms xi yi zi real
% var=[1,xi,yi,zi,xi^2,yi^2,zi^2,xi*yi,yi*zi,xi*zi,xi^3,yi^3,zi^3,xi^2*yi,xi^2*zi,yi^2*zi,yi^2*xi,zi^2*xi,zi^2*yi,xi*yi*zi];
var=[1,xi,yi,zi,xi*yi,yi*zi,xi*zi,xi*yi*zi];

for shpfunc=1:8
    b=zeros(8,1);
    b(shpfunc)=1;
    A=zeros(8,8);
    for node=1:8
        A(node,1:8)=double(subs(var, [xi, yi, zi], P(node,:)));
    end
    coef=A\b;
    N(shpfunc)=var*coef;
end
syms x y z real
x=N*P1(:,1);%+10^-12*xi+10^-12*yi+10^-12*zi;
y=N*P1(:,2);%+10^-12*xi+10^-12*yi+10^-12*zi;
z=N*P1(:,3);%+10^-12*xi+10^-12*yi+10^-12*zi;
J=double([diff(x,xi),diff(x,yi),diff(x,zi);diff(y,xi),diff(y,yi),diff(y,zi);diff(z,xi),diff(z,yi),diff(z,zi)]);
invJ=inv(J);
for shpfunc=1:8
    DN(shpfunc,:)=[diff(N(shpfunc),xi),diff(N(shpfunc),yi),diff(N(shpfunc),zi)];
end
for shpfunc=1:8
    DNx(shpfunc)=DN(shpfunc,:)*invJ(:,1);
    DNy(shpfunc)=DN(shpfunc,:)*invJ(:,2);
    DNz(shpfunc)=DN(shpfunc,:)*invJ(:,3);
end 
% 
B(1,:)=[DNx(1),sym(0),sym(0),DNx(2),sym(0),sym(0),DNx(3),sym(0),sym(0),DNx(4),sym(0),sym(0),DNx(5),sym(0),sym(0),DNx(6),sym(0),sym(0),DNx(7),sym(0),sym(0),DNx(8),sym(0),sym(0)];
B(2,:)=[sym(0),DNy(1),sym(0),sym(0),DNy(2),sym(0),sym(0),DNy(3),sym(0),sym(0),DNy(4),sym(0),sym(0),DNy(5),sym(0),sym(0),DNy(6),sym(0),sym(0),DNy(7),sym(0),sym(0),DNy(8),sym(0)];
B(3,:)=[sym(0),sym(0),DNz(1),sym(0),sym(0),DNz(2),sym(0),sym(0),DNz(3),sym(0),sym(0),DNz(4),sym(0),sym(0),DNz(5),sym(0),sym(0),DNz(6),sym(0),sym(0),DNz(7),sym(0),sym(0),DNz(8)];
B(4,:)=[DNy(1),DNx(1),sym(0),DNy(2),DNx(2),sym(0),DNy(3),DNx(3),sym(0),DNy(4),DNx(4),sym(0),DNy(5),DNx(5),sym(0),DNy(6),DNx(6),sym(0),DNy(7),DNx(7),sym(0),DNy(8),DNx(8),sym(0)];
B(5,:)=[sym(0),DNz(1),DNy(1),sym(0),DNz(2),DNy(2),sym(0),DNz(3),DNy(3),sym(0),DNz(4),DNy(4),sym(0),DNz(5),DNy(5),sym(0),DNz(6),DNy(6),sym(0),DNz(7),DNy(7),sym(0),DNz(8),DNy(8)];
B(6,:)=[DNz(1),sym(0),DNx(1),DNz(2),sym(0),DNx(2),DNz(3),sym(0),DNx(3),DNz(4),sym(0),DNx(4),DNz(5),sym(0),DNx(5),DNz(6),sym(0),DNx(6),DNz(7),sym(0),DNx(7),DNz(8),sym(0),DNx(8)];
BTDB=B'*D*B;
BTD=B'*D;
BTDB_func=matlabFunction(BTDB+10^-12*xi+10^-12*yi+10^-12*zi, 'Vars', [xi, yi, zi]);
BTD_func=matlabFunction(BTD+10^-12*xi+10^-12*yi+10^-12*zi, 'Vars', [xi, yi, zi]);
N_func=matlabFunction(N+10^-12*xi+10^-12*yi+10^-12*zi, 'Vars', [xi, yi, zi]);
%% Integration
GPw=[ -0.9815606342,0.0471753364;-0.9041172564,0.1069393259;-0.7699026741,0.1600783285;-0.5873179542,0.2031674267;-0.3678314989,0.2334925365;-0.1252334085,0.2491470458;0.1252334085,0.2491470458;0.3678314989,0.2334925365;0.5873179542,0.2031674267;0.7699026741,0.1600783285;0.9041172564,0.1069393259;0.9815606342,0.0471753364];
BTDBv=zeros(24,24);
BTDv=zeros(24,6);
N_funcv=zeros(1,8);
for i=1:length(GPw(:,1))
    for j=1:length(GPw(:,1))
        for k=1:length(GPw(:,1))
            BTDBijk=GPw(i,2)*GPw(k,2)*GPw(k,2)*BTDB_func(GPw(i,1),GPw(j,1),GPw(k,1));
            BTDBv=BTDBv+BTDBijk;
        end
    end
end
Ke=BTDBv*det(J);
clear BTDBv
end