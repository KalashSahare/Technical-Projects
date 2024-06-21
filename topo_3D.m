%% Topology
clear all
close all
% Rectangular Domain
% sq bilinear element
%% initial Data
origin=[0,0,0];
L=1;%m
B=1;%m
H=1;%m
%% Material Properties
E=200*10^9;%N/m2
mu=0.3;
D=(E/((1+mu)*(1-2*mu)))*[1-mu,mu,mu,0,0,0;mu,1-mu,mu,0,0,0;mu,mu,1-mu,0,0,0;0,0,0,(1-2*mu)/2,0,0;0,0,0,0,(1-2*mu)/2,0;0,0,0,0,0,(1-2*mu)/2];
%% Discretisation
N_x=10;
N_y=10;
N_z=10;
Eledx=L/N_x;
Eledy=B/N_y;
Eledz=H/N_z;
Tot_ele=N_x*N_y*N_z;
Tot_node=(N_x+1)*(N_y+1)*(N_z+1);
Tot_dof=3*Tot_node;
V=L*B*H; % domain area/ volume
Ve=Eledx*Eledy*Eledz; % area/volume of element
Vf=0.25;
penal=3;
rho=Vf*ones(Tot_ele,1);
rhomin=0.001;
r_fil=sqrt(Eledx^2+Eledy^2+Eledz^2);

%% Meshing
[Nodepos,Tot_node,MeshData]=Meshing_1oHex(origin,N_x,N_y,N_z,L,B,H);
Center = zeros(Tot_ele,3);
for i=1:Tot_ele
%     Center(i,:)=0.5*(Nodepos(MeshData(i,1),:)+Nodepos(MeshData(i,7),:));
    Center(i,:)=Nodepos(MeshData(i,1),:)+0.5*[Eledx,Eledy,Eledz];
end
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
Pc=[load/4,(N_x+1),1;load/4,(N_x+1)*(N_y+1),1;load/4,(N_x+1)*(N_y+1)*N_z+(N_x+1),1;load/4,(N_x+1)*(N_y+1)*(N_z+1),1];
% for i=1:length(fEN(:,1))
%     Pc(i,:)=[load/length(fEN(:,1)),fEN(i,1),1];
% end
%% Local stifness matrix in local coordinates
% as our basis funxtion are cubic then in x y z so 
% As straight edge element the mapping would be liner in xi yi zi

 for j=1:length(MeshData(1,:))
     P1(j,:)=Nodepos(MeshData(1,j),:);
 end
[Ke]=MasterElement(P1,D);

%% Concentrated load vector
CPL=zeros(3*Tot_node,1);
for i=1:length(Pc(:,1))
    CPL(3*Pc(i,2)-3+Pc(i,3),1)=Pc(i,1);
end
%% topology Optimisation
change=1;
itr=0;
U=zeros(Tot_dof,1);
dCdrho=zeros(Tot_ele,1);
while change>0.01
    itr=itr+1; 
    rho_old=rho;
    %% FEM solver
    KG=zeros(Tot_dof,Tot_dof);
    for N=1:Tot_ele
        n1=MeshData(N,1);
        n2=MeshData(N,2);
        n3=MeshData(N,3);
        n4=MeshData(N,4);
        n5=MeshData(N,5);
        n6=MeshData(N,6);
        n7=MeshData(N,7);
        n8=MeshData(N,8);
        e_dof=[3*n1-2;3*n1-1;3*n1;3*n2-2;3*n2-1;3*n2;3*n3-2;3*n3-1;3*n3;3*n4-2;3*n4-1;3*n4;3*n5-2;3*n5-1;3*n5;3*n6-2;3*n6-1;3*n6;3*n7-2;3*n7-1;3*n7;3*n8-2;3*n8-1;3*n8];
        KG(e_dof,e_dof)=KG(e_dof,e_dof)+rho(N,1)^penal*Ke;
    end
    %Applying boundary condition
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
    %Solution
    U=KG\CPL;
    
    %% calculating Strain Energy and Sentivity Analysis
    C=0;
    for i=1:Tot_ele
        n1=MeshData(i,1);
        n2=MeshData(i,2);
        n3=MeshData(i,3);
        n4=MeshData(i,4);
        n5=MeshData(i,5);
        n6=MeshData(i,6);
        n7=MeshData(i,7);
        n8=MeshData(i,8);
        Ue=U([3*n1-2;3*n1-1;3*n1;3*n2-2;3*n2-1;3*n2;3*n3-2;3*n3-1;3*n3;3*n4-2;3*n4-1;3*n4;3*n5-2;3*n5-1;3*n5;3*n6-2;3*n6-1;3*n6;3*n7-2;3*n7-1;3*n7;3*n8-2;3*n8-1;3*n8],1);
        C=C+(rho(i,1)^(penal)*Ue'*Ke*Ue);
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
    for ele=1:Tot_ele
        [F1,F2,F3,F4,F5,F6]=baseplot(Eledx,Eledy,Eledz,Nodepos(MeshData(ele,1),:));
        fill3(F1(:,1),F1(:,2),F1(:,3),'k','FaceAlpha',rho(ele,1))
        hold on
        fill3(F2(:,1),F2(:,2),F2(:,3),'k','FaceAlpha',rho(ele,1))
        hold on
        fill3(F3(:,1),F3(:,2),F3(:,3),'k','FaceAlpha',rho(ele,1))
        hold on
        fill3(F4(:,1),F4(:,2),F4(:,3),'k','FaceAlpha',rho(ele,1))
        hold on
        fill3(F5(:,1),F5(:,2),F5(:,3),'k','FaceAlpha',rho(ele,1))
        hold on
        fill3(F6(:,1),F6(:,2),F6(:,3),'k','FaceAlpha',rho(ele,1))
        hold on
    end
    axis equal
    hold off
    pause(0.001)
end

%% function
function [rho_n]= UpdateOptimal(rho,dCdrho,Ve,Vf,rhomin,Tot_ele,V)
lam1=0;
lam2=10^9;
xi=0.2; % move factor
neta=0.5;
rho_n=zeros(Tot_ele,1);
q=1;
while ((lam2-lam1)/(lam2+lam1))>10^-3
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
            rho_n(ele,1)=rhoBn(ele,1)^q;
        end
    end
    if sum(rho_n)*Ve>Vf*V
        lam1=lam;
    else
        lam2=lam;
    end
end
end

function [dCdrho_new]= FilterSensitivityAna(rho,dCdrho,Centre,Tot_ele,r_filter)
dCdrho_new=zeros(Tot_ele,1);
for ele=1:Tot_ele
    H=zeros(Tot_ele,1);
    for ele_i=1:Tot_ele
        dist=sqrt((Centre(ele,1)-Centre(ele_i,1))^2+(Centre(ele,2)-Centre(ele_i,2))^2+(Centre(ele,3)-Centre(ele_i,3))^2);
        if dist<=r_filter 
            H(ele_i,1)=r_filter-dist;
        end
    end
    sumHi=sum(H(:,1));
    sumHiRidci=sum(H.*rho.*dCdrho);
    dCdrho_new(ele,1)=sumHiRidci/(rho(ele,1)*sumHi);    
end
clear dist H sumHi sumHiRidci
end

function [F1,F2,F3,F4,F5,F6]= baseplot(L,B,H,origin)
P1=origin;
P2=P1+[L,0,0];
P3=P1+[L,B,0];
P4=P1+[0,B,0];
P5=P1+[0,0,H];
P6=P1+[L,0,H];
P7=P1+[L,B,H];
P8=P1+[0,B,H];
F1=[P1;P2;P3;P4];
F2=[P4;P3;P7;P8];
F3=[P5;P6;P7;P8];
F4=[P1;P2;P6;P5];
F5=[P1;P4;P8;P5];
F6=[P2;P3;P7;P6];
% cord=[P1;P2;P3;P4;P5;P6;P7;P8];
end