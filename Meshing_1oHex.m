function [P1,Totnode,MD]=Meshing_1oHex(origin,Nx,Ny,Nz,L,B,H)
dx=L/Nx;
dy=B/Ny;
dz=H/Nz;
count1=0;
for k=1:Nz+1
    for j=1:Ny+1
        for i=1:Nx+1
            count1=count1+1;
            P1(count1,:)=origin+[(i-1)*dx,(j-1)*dy,(k-1)*dz];
        end
    end
end
Totnode=count1;                
ele=1;
% corner node
for k=1:Nz
    for j=1:Ny
        for i=1:Nx
            indz=(k-1)*(Nx+1)*(Ny+1);
            indz1=k*(Nx+1)*(Ny+1);
            indy=(j-1)*(Nx+1);
            ind=indy+indz;
            ind1=indy+indz1;
            MD(ele,1:8)=[ind+i,ind+i+1,ind+i+Nx+1+1,ind+i+Nx+1,ind1+i,ind1+i+1,ind1+i+Nx+1+1,ind1+i+Nx+1];
            ele=ele+1;
        end
    end
end
end