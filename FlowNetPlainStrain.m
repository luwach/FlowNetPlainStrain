
clc;
clear all;

%Parameters of the model

Lx = input('Insert width of the model Lx [cm]: ');
Ly = input('Insert height of the model Ly (soil) [cm]: ');
D = input('Insert diameter of the caisson D [cm]: ');
H = input('Insert height of the caisson H [cm]: ');
w = input('Insert parameter w [cm]: ');
pw = input('Insert inner pressure pw [m]: ');
pz = input('Insert outer pressure pz [m]: ');
d = input('Insert nodes spacing [cm]: ');

%Number of nodes

nw=w/d;
nx=Lx/d+1;
nu=(5/d+1);
ny=Ly/d+1+nw+nu;
nk=D/d+1;
nx1=(nx-nk)/2;
hk=H/d+1;
    
h=zeros(ny,nx);

for i = 1:ny
    for j = 1:nx
        
        h(i,j)= -Inf;
        
    end    
end

%Potential lines

%BOUNDARY CONDITIONS

%Left nodes of the model

for i = ny-nw-nu+1:ny
    for j = 1:nx1
        
        h(i,j)=NaN;
        
    end    
end

%Right nodes of the model

for i = ny-nw-nu+1:ny
    for j = nx1+nk+1:nx
        
        h(i,j)=NaN;
        
    end    
end

%Center nodes of the model

for i = ny-nu+1:ny
    for j = nx1+1:nx1+nk
        
        h(i,j)=NaN;
        
    end    
end

%Pressure inside of the model

h(ny-nu,nx1+1:nx1+nk)=linspace(pw,pw,nk);

%Pressure outside of the model

h(ny-nw-nu,1:nx1)=linspace(pz,pz,nx1);
h(ny-nw-nu,nx1+nk+1:nx)=linspace(pz,pz,nx1);

%Pressure on the inner boundary of the caisson

h(ny-hk-nu+2:ny-nu,nx1+1)=linspace(0,pw,hk-1);
h(ny-hk-nu+2:ny-nu,nx1+nk)=linspace(0,pw,hk-1);

%Pressure on the outer boundary of the caisson

h(ny-nu-hk+2:ny-nu-nw,nx1)=linspace(pz+1,pz,hk-nw-1);
h(ny-nu-hk+2:ny-nu-nw,nx1+nk+1)=linspace(pz+1,pz,hk-nw-1);

%Pressure on the inner and outer boundary of the model

h(1:ny-nu-nw,1)=linspace(pz+1,pz,ny-nu-nw);
h(1:ny-nu-nw,nx)=linspace(pz+1,pz,ny-nu-nw);

%Pressure at the bottom boundary of the model

h(1,1:nx)=linspace(pz+1,pz+1,nx);

h