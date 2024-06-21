function [Ke]=MasterSQBLElement(P1,D)
% E=1;
% mu=0.3;
% D=(E/(1-mu^2))*[1,mu,0;mu,1,0;0,0,0.5*(1-mu)];
%% Master Element
% P is point of master element
% P(1,:)=[-1,-1];
% P(2,:)=[1,-1];
% P(3,:)=[1,1];
% P(4,:)=[-1,1];
% % P1 is points of physical element
% % P1(1,:)=[0,0];
% % P1(2,:)=[10,0];
% % P1(3,:)=[10,5];
% % P1(4,:)=[0,5];
syms xi yi real
N(1,1)=0.25*(1-xi)*(1-yi);
N(1,2)=0.25*(1+xi)*(1-yi);
N(1,3)=0.25*(1+xi)*(1+yi);
N(1,4)=0.25*(1-xi)*(1+yi);

syms x y real
x=N*P1(:,1);%+10^-12*xi+10^-12*yi+10^-12*zi;
y=N*P1(:,2);%+10^-12*xi+10^-12*yi+10^-12*zi;
J=double([diff(x,xi),diff(x,yi);diff(y,xi),diff(y,yi)]);
invJ=inv(J);
for shpfunc=1:4
    DN(shpfunc,:)=[diff(N(shpfunc),xi),diff(N(shpfunc),yi)];
end
for shpfunc=1:4
    DNx(shpfunc)=DN(shpfunc,:)*invJ(:,1);
    DNy(shpfunc)=DN(shpfunc,:)*invJ(:,2);
end 
% 
B(1,:)=[DNx(1),sym(0),DNx(2),sym(0),DNx(3),sym(0),DNx(4),sym(0)];
B(2,:)=[sym(0),DNy(1),sym(0),DNy(2),sym(0),DNy(3),sym(0),DNy(4)];
B(3,:)=[DNy(1),DNx(1),DNy(2),DNx(2),DNy(3),DNx(3),DNy(4),DNx(4)];
BTDB=B'*D*B;
BTD=B'*D;
BTDB_func=matlabFunction(BTDB+10^-12*xi+10^-12*yi, 'Vars', [xi, yi]);
BTD_func=matlabFunction(BTD+10^-12*xi+10^-12*yi, 'Vars', [xi, yi]);
N_func=matlabFunction(N+10^-12*xi+10^-12*yi, 'Vars', [xi, yi]);
%% Integration
% a shape function are bilinear so 2nd order guassian integration is used
GPw=[ -1/sqrt(3),1;1/sqrt(3),1];
BTDBv=zeros(8,8);
BTDv=zeros(8,3);
N_funcv=zeros(1,8);
for i=1:length(GPw(:,1))
    for j=1:length(GPw(:,1))
        BTDBijk=GPw(i,2)*GPw(j,2)*BTDB_func(GPw(i,1),GPw(j,1));
        BTDBv=BTDBv+BTDBijk;
    end
end
Ke=BTDBv*det(J);
clear BTDBv
end