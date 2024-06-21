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