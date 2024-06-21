clear all
close all
tic
%% FEM Project
% 2D beam/plane stress analysis sig zz==sig xz==sig yz==0
% 2D square Bilinear element 
% 4 nodes
% each node 2 dof u,v
% model Length L, width W
% IP of the type--- bilinear in x and y (a1+a2x)*(a3+a4y)
%% Material Properties
E=200*10^9;%N/m2
mu=0.3;
cond=0; % 0 is for plane stress 1 is for plane strain
if cond==0
    D=(E/(1-mu^2))*[1,mu,0;mu,1,0;0,0,0.5*(1-mu)];
elseif cond==1
    D=(E/((1+mu)*(1-2*mu)))*[1-mu,mu,0;mu,1-mu,0;0,0,0.5*(1-2*mu)];
end
%% Geometry
% cantilever beam with length L width W
% fixed end x=0
% free end x=L
origin=[0,0];
L=10;%m
W=1;%m
figure
Sqplot(L,W,origin)
hold on
%% Discreatisation
N_x=10;
N_y=10;
Eledx=L/N_x;
Eledy=W/N_y;
Tot_ele=N_x*N_y;
Tot_node=(N_x+1)*(N_y+1);
[Nodepos,MeshData]=Meshing_SqBilinear(origin,N_x,N_y,L,W);
plot(Nodepos(:,1),Nodepos(:,2),'*')
hold on
for ele=1:Tot_ele
    Sqplot(Eledx,Eledy,Nodepos(MeshData(ele,1),:))
end
%% Local Stiffness matrix

    for j=1:length(MeshData(1,:))
        P1(j,:)=Nodepos(MeshData(1,j),:);
    end
    [Ke]=MasterSQBLElement(P1,D);
%% K local to K global
KG=zeros(2*Tot_node,2*Tot_node);
KgE=zeros(2*Tot_node,2*Tot_node,Tot_ele);
for N=1:Tot_ele
    [KgE(:,:,N)]=local_TO_global(MeshData(N,:),Ke,Tot_node);
    KG=KG+KgE(:,:,N);
end
clear Ke
pause()
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
% Pc=[load,5,3;load,10,3];% eg 10N, node no 10, direction x  [10,10,1]
for i=1:length(fEN(:,1))
    Pc(i,:)=[load/length(fEN(:,1)),fEN(i,1),2];
end
%% Concentrated load vector
CPL=zeros(2*Tot_node,1);
for i=1:length(Pc(:,1))
    CPL(2*Pc(i,2)-2+Pc(i,3),1)=Pc(i,1);
end

%% Applying boundary condition
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
%% Solution
Q=KG\CPL;
%% Post Processing
tipdefA=load*L^3/(3*E*(W*1^3/12));
tipdef=0;
for i=1:length(fEN(:,1))
    tipdef=tipdef+Q(2*fEN(i,1),1);
end
tipdef=tipdef/length(fEN(:,1));

Error=(tipdefA-tipdef)*100/tipdefA; % percent error
%% Plotting

for i=1:Tot_node
    NewNodePos(i,:)=Nodepos(i,:)+[Q(2*i-1,1),Q(2*i,1)];
end
plot(NewNodePos(:,1),NewNodePos(:,2),'*r')
hold off
toc
%% Function
function [Kg]=local_TO_global(MD,ke,Ntotal)
Kg=zeros(2*Ntotal,2*Ntotal);
for i=1:4
    gdof(2*i-1)=2*MD(1,i)-1;
    gdof(2*i)=2*MD(1,i);
end
for ldofi=1:8
    for ldofj=1:8
        gdofi=gdof(ldofi);
        gdofj=gdof(ldofj);
        Kg(gdofi,gdofj)=ke(ldofi,ldofj);
    end
end
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