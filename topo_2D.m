%% Topology
clear all
close all
% Rectangular Domain
% sq bilinear element
%% initial Data
origin=[0,0];
L=1;
W=1;
N_x=40;
N_y=40;
Eledx=L/N_x;
Eledy=W/N_y;
Tot_ele=N_x*N_y;
Tot_node=(N_x+1)*(N_y+1);
Tot_dof=2*Tot_node;
V=L*W; % domain area/ volume
Ve=Eledx*Eledy; % area/volume of element
Vf=0.2;
penal=3;
rho=Vf*ones(N_x*N_y,1);
rhomin=0.001;
r_fil=sqrt(Eledx^2+Eledy^2);
%% Material Properties
E=200*10^9;%N/m2
mu=0.3;
cond=0; % 0 is for plane stress 1 is for plane strain
if cond==0
    D=(E/(1-mu^2))*[1,mu,0;mu,1,0;0,0,0.5*(1-mu)];
elseif cond==1
    D=(E/((1+mu)*(1-2*mu)))*[1-mu,mu,0;mu,1-mu,0;0,0,0.5*(1-2*mu)];
end
%% Meshing
[Nodepos,MeshData]=Meshing_SqBilinear(origin,N_x,N_y,L,W);
Center = zeros(Tot_ele,2);
for i=1:Tot_ele
    Center(i,:)=0.5*(Nodepos(MeshData(i,1),:)+Nodepos(MeshData(i,3),:));
end
%% Boundary 
% Cantelevier beam with x=0 fixed end
% FEN is fixed end node
count=1;
    for i=1:N_y+1
        FEN(count,1)=1+(i-1)*(N_x+1);
        count=count+1;
    end
clear count
% fEN free end nodes
count=1;
    for i=1:N_y+1
        fEN(count,1)=i*(N_x+1);
        count=count+1;
    end
clear count
%% Forces 
% Concentrated Point load PL specify the node and load value and direction x--1, y--2
% taking the condition of end loaded concentrated beam
load=1*10^6; % N
% Pc=[-load,0.5*N_y*(N_x+1),2];% eg 10N, node no 10, direction x  [10,10,1]
for i=1:length(fEN(:,1))
    Pc(i,:)=[load/length(fEN(:,1)),fEN(i,1),1];
end
%% Local Stiffness matrix
for j=1:length(MeshData(1,:))
    P1(j,:)=Nodepos(MeshData(1,j),:);
end
[Ke]=MasterSQBLElement(P1,D);
%% K local to K global
% KgE=zeros(2*Tot_node,2*Tot_node,Tot_ele);
% for N=1:Tot_ele
%     [KgE(:,:,N)]=local_TO_global(MeshData(N,:),Ke,Tot_node);
% end
% clear Ke
%% Concentrated load vector
CPL=zeros(2*Tot_node,1);
for i=1:length(Pc(:,1))
    CPL(2*Pc(i,2)-2+Pc(i,3),1)=Pc(i,1);
end
%% topology Optimisation
change=1;
itr=0;
U=zeros(2*Tot_node,1);
dCdrho=zeros(Tot_ele,1);
while change>0.01
    itr=itr+1; 
    rho_old=rho;
    %% FEM solver
    KG=zeros(2*Tot_node,2*Tot_node);
    for N=1:Tot_ele
        n1=MeshData(N,1);
        n2=MeshData(N,2);
        n3=MeshData(N,3);
        n4=MeshData(N,4);
%         [KgE]=local_TO_global(MeshData(N,:),Ke,Tot_node);
        e_dof=[2*n1-1;2*n1;2*n2-1;2*n2;2*n3-1;2*n3;2*n4-1;2*n4];
        KG(e_dof,e_dof)=KG(e_dof,e_dof)+rho(N,1)^penal*Ke;
    end
    %Applying boundary condition
    for i=1:length(FEN(:,1))
        node=FEN(i,1);
        dofx=2*node-1;
        dofy=2*node;
        KG(dofx,:)=0;
        KG(:,dofx)=0;
        KG(dofx,dofx)=1;
        CPL(dofx,1)=0;
    
        KG(dofy,:)=0;
        KG(:,dofy)=0;
        KG(dofy,dofy)=1;
        CPL(dofy,1)=0;
    end
    %Solution
    U=KG\CPL;
    
    %% calculating Strain Energy and Sentivity Analysis
    C=0;
    for i=1:Tot_ele
        n1=MeshData(i,1);
        n2=MeshData(i,2);
        n3=MeshData(i,3);
        n4=MeshData(i,4);
        Ue=U([2*n1-1;2*n1;2*n2-1;2*n2;2*n3-1;2*n3;2*n4-1;2*n4],1);
        C=C+rho(i,1)^(penal)*Ue'*Ke*Ue;
        dCdrho(i,1)=-penal*rho(i,1)^(penal-1)*Ue'*Ke*Ue;
    end
    %%  Filtering
    [dCdrho]=FilterSensitivityAna(rho,dCdrho,Center,Tot_ele,r_fil);
    %% Optimality and Updating
    [rho]= UpdateOptimal(rho,dCdrho,Ve,Vf,rhomin,Tot_ele,V);
    change=max(abs(rho-rho_old));
    disp([' It.: ' sprintf('%4i',itr) ' Obj.: ' sprintf('%10.4f',C) ...
       ' Vol.: ' sprintf('%6.3f',sum(rho)*Ve) ...
        ' ch.: ' sprintf('%6.3f',change )])
    %% plotting
%     plot(Nodepos(:,1),Nodepos(:,2),'*')
%     hold on
    for ele=1:Tot_ele
        [cord]=Sqplot1(Eledx,Eledy,Nodepos(MeshData(ele,1),:));
        fill(cord(:,1),cord(:,2),'k','FaceAlpha',rho(ele,1))
        hold on
    end
    axis equal
    axis tight
    hold off
    pause(0.001)
end

%% function
function [rho_n]= UpdateOptimal(rho,dCdrho,Ve,Vf,rhomin,Tot_ele,V)
lam1=0;
lam2=10000;
xi=0.2; % move factor
neta=0.5;
rho_n=zeros(Tot_ele,1);
while (lam2-lam1)>10^-4
    lam=0.5*(lam1+lam2);
    B=-dCdrho./(lam*Ve);
    Bn=B.^neta;
    rhoBn=rho.*Bn;
    for ele=1:Tot_ele
        if rhoBn(ele,1)<=max(rhomin,rho(ele,1)*(1-xi))
            rho_n(ele,1)=max(rhomin,rho(ele,1)*(1-xi));
        elseif rhoBn(ele,1)>=min(1,rho(ele,1)*(1+xi))
            rho_n(ele,1)=min(1,rho(ele,1)*(1+xi));
        else
            rho_n(ele,1)=rhoBn(ele,1);
        end
    end
    if sum(rho_n)*Ve>Vf*V
        lam1=lam;
    else
        lam2=lam;
    end
end
end
function [P1,MD]=Meshing_SqBilinear(origin,Nx,Ny,L,B)
dx=L/Nx;
dy=B/Ny;
count1=0;
    for j=1:Ny+1
        for i=1:Nx+1
            count1=count1+1;
            P1(count1,:)=origin+[(i-1)*dx,(j-1)*dy];
        end
    end             
ele=1;
    for j=1:Ny
        for i=1:Nx
            indy=(j-1)*(Nx+1);
            ind=indy;
            MD(ele,1:4)=[ind+i,ind+i+1,ind+i+Nx+1+1,ind+i+Nx+1];
            ele=ele+1;
        end
    end

end

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
BTDB_func=matlabFunction(BTDB+10^-12*xi+10^-12*yi, 'Vars', [xi, yi]);

%% Integration
% a shape function are bilinear so 2nd order guassian integration is used
GPw=[ -1/sqrt(3),1;1/sqrt(3),1];
BTDBv=zeros(8,8);
for i=1:length(GPw(:,1))
    for j=1:length(GPw(:,1))
        BTDBijk=GPw(i,2)*GPw(j,2)*BTDB_func(GPw(i,1),GPw(j,1));
        BTDBv=BTDBv+BTDBijk;
    end
end
Ke=BTDBv*det(J);
clear BTDBv N x y J invJ B BTDB BTDB_func DN DNx DNy
end

function [dCdrho_new]=FilterSensitivityAna(rho,dCdrho,Centre,Tot_ele,r_filter)
dCdrho_new=zeros(Tot_ele,1);
for ele=1:Tot_ele
    H=zeros(Tot_ele,1);
    for ele_i=1:Tot_ele
        dist=sqrt((Centre(ele,1)-Centre(ele_i,1))^2+(Centre(ele,2)-Centre(ele_i,2))^2);
        if dist<= r_filter
            H(ele_i,1)=(r_filter-dist)/r_filter;
        end
    end
    sumHi=sum(H(:,1));
    sumHiRidci=sum(H.*rho.*dCdrho);
    dCdrho_new(ele,1)=sumHiRidci/(rho(ele,1)*sumHi);    
end
clear dist H sumHi sumHiRidci
end

function [cord]=Sqplot1(L,B,origin)
P1=origin;
P2=P1+[L,0];
P3=P1+[L,B];
P4=P1+[0,B];
cord=[P1;P2;P3;P4];
end

function Sqplot(L,B,origin)
P1=origin;
P2=P1+[L,0];
P3=P1+[L,B];
P4=P1+[0,B];
L1=[P1;P2];
L2=[P2;P3];
L3=[P3;P4];
L4=[P4;P1];
plot(L1(:,1),L1(:,2),'b')
hold on
plot(L2(:,1),L2(:,2),'b')
hold on
plot(L3(:,1),L3(:,2),'b')
hold on
plot(L4(:,1),L4(:,2),'b')
hold on
end

