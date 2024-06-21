clear all
close all
tic
%% FEM Project
% 3D beam 
% 3D hexahedral element 
% 8 nodes
% each node 3 dof u,v,w
% model Length L, breath B height W
% IP of the type--- cubic polynomial a1+a2x+a3y+a4z+a5xy+a6yz+a7zx+a8xyz
%% Material Properties
E=200*10^9;%N/m2
mu=0.3;
D=(E/((1+mu)*(1-2*mu)))*[1-mu,mu,mu,0,0,0;mu,1-mu,mu,0,0,0;mu,mu,1-mu,0,0,0;0,0,0,(1-2*mu)/2,0,0;0,0,0,0,(1-2*mu)/2,0;0,0,0,0,0,(1-2*mu)/2];

%% Geometery
origin=[0,0,0];
L=10;%m
B=1;%m
H=1;%m
figure
baseplot(L,B,H,origin)
hold on
%% Discretisation
N_x=100;
N_y=4;
N_z=4;
Eledx=L/N_x;
Eledy=B/N_y;
Eledz=H/N_z;
Tot_ele=N_x*N_y*N_z;
[Nodpos,Tot_node,MeshData]=Meshing_1oHex(origin,N_x,N_y,N_z,L,B,H);
plot3(Nodpos(:,1),Nodpos(:,2),Nodpos(:,3),'*')
hold on
for ele=1:Tot_ele
    baseplot(Eledx,Eledy,Eledz,Nodpos(MeshData(ele,1),:))
end
axis equal
hold on

%% Boundary 
% Cantelevier beam with x=0 fixed end
% FEN is fixed end node
count=1;

for j=1:N_z+1
    for i=1:N_y+1
        FEN(count,1)=1+(j-1)*(N_x+1)*(N_y+1)+(i-1)*(N_x+1);
        count=count+1;
    end
end
clear count
% fEN free end nodes
count=1;
for j=1:N_z+1
    for i=1:N_y+1
        fEN(count,1)=(j-1)*(N_x+1)*(N_y+1)+i*(N_x+1);
        count=count+1;
    end
end
clear count
%% Forces
% body forces PHI
PHIx=0;
PHIy=0;
PHIz=0;
% distributed surface forces phi
phix=0;
phiy=0;
phiz=0;
% Concentrated Point load PL specify the node and load value and direction x--1, y--2, z--3
% taking the condition of end loaded concentrated beam
load=1*10^6; % N
% Pc=[load,5,3;load,10,3];% eg 10N, node no 10, direction x  [10,10,1]
for i=1:length(fEN(:,1))
    Pc(i,:)=[load/length(fEN(:,1)),fEN(i,1),3];
end
%% Mapping to master element
% as our basis funxtion are cubic then in x y z so 
% As straight edge element the mapping would be liner in xi yi zi

 for j=1:length(MeshData(1,:))
     P1(j,:)=Nodpos(MeshData(1,j),:);
 end
[Ke]=MasterElement(P1,D);


%% K local to K global
KG=zeros(3*Tot_node,3*Tot_node);
for N=1:Tot_ele
    [KgN]=local_TO_global(MeshData(N,:),Ke(:,:),Tot_node);
    KG=KG+KgN;
    clear KgN
end

%% Concentrated load vector
CPL=zeros(3*Tot_node,1);
for i=1:length(Pc(:,1))
    CPL(3*Pc(i,2)-3+Pc(i,3),1)=Pc(i,1);
end
%% Applying boundary condition
for i=1:length(FEN(:,1))
    node=FEN(i,1);
    dofx=3*node-2;
    dofy=3*node-1;
    dofz=3*node;
    KG(dofx,:)=0;
    KG(:,dofx)=0;
    KG(dofx,dofx)=1;
    CPL(dofx,1)=0;
    
    KG(dofy,:)=0;
    KG(:,dofy)=0;
    KG(dofy,dofy)=1;
    CPL(dofy,1)=0;
    
    KG(dofz,:)=0;
    KG(:,dofz)=0;
    KG(dofz,dofz)=1;
    CPL(dofz,1)=0;
end
%% Solution
Q=KG\CPL;
%% Post Processing
tipdefA=load*L^3/(3*E*(B*H^3/12));
tipdef=0;
for i=1:length(fEN(:,1))
    tipdef=tipdef+Q(3*fEN(i,1),1);
end
tipdef=tipdef/length(fEN(:,1));

Error=(tipdefA-tipdef)*100/tipdefA; % percent error
%% Plotting

for i=1:Tot_node
    NewNodePos(i,:)=Nodpos(i,:)+[Q(3*i-2,1),Q(3*i-1,1),Q(3*i,1)];
end
plot3(NewNodePos(:,1),NewNodePos(:,2),NewNodePos(:,3),'*r')
hold off
toc
%% Functions

function baseplot(L,B,H,origin)
P1=origin;
P2=P1+[L,0,0];
P3=P1+[L,B,0];
P4=P1+[0,B,0];
P5=P1+[0,0,H];
P6=P1+[L,0,H];
P7=P1+[L,B,H];
P8=P1+[0,B,H];
L1=[P1;P2];
L2=[P3;P2];
L3=[P4;P3];
L4=[P1;P4];
L5=[P5;P1];
L6=[P6;P2];
L7=[P7;P3];
L8=[P8;P4];
L9=[P6;P5];
L10=[P7;P6];
L11=[P8;P7];
L12=[P5;P8];
plot3(L1(:,1),L1(:,2),L1(:,3),'b')
hold on
plot3(L2(:,1),L2(:,2),L2(:,3),'b')
hold on
plot3(L3(:,1),L3(:,2),L3(:,3),'b')
hold on
plot3(L4(:,1),L4(:,2),L4(:,3),'b')
hold on
plot3(L5(:,1),L5(:,2),L5(:,3),'b')
hold on
plot3(L6(:,1),L6(:,2),L6(:,3),'b')
hold on
plot3(L7(:,1),L7(:,2),L7(:,3),'b')
hold on
plot3(L8(:,1),L8(:,2),L8(:,3),'b')
hold on
plot3(L9(:,1),L9(:,2),L9(:,3),'b')
hold on
plot3(L10(:,1),L10(:,2),L10(:,3),'b')
hold on
plot3(L11(:,1),L11(:,2),L11(:,3),'b')
hold on
plot3(L12(:,1),L12(:,2),L12(:,3),'b')
hold on
end

function [Kg]=local_TO_global(MD,ke,Ntotal)
Kg=zeros(3*Ntotal,3*Ntotal);
% for lni=1:8
%     gni=MD(1,lni);
%     for lnj=1:8
%         gnj=MD(1,lnj);
%         for i=1:3
%             ldofi=3*lni-(i-1);
%             gdofi=3*gni-(i-1);
%             for j=1:3
%                 ldofj=3*lnj-(j-1);
%                 gdofj=3*gnj-(j-1);
%                 Kg(gdofi,gdofj)=ke(ldofi,ldofj);
%             end
%         end
%     end
% end
for i=1:8
    gdof(3*i-2)=3*MD(1,i)-2;
    gdof(3*i-1)=3*MD(1,i)-1;
    gdof(3*i)=3*MD(1,i);
end
for ldofi=1:24
    for ldofj=1:24
        gdofi=gdof(ldofi);
        gdofj=gdof(ldofj);
        Kg(gdofi,gdofj)=ke(ldofi,ldofj);
    end
end
end