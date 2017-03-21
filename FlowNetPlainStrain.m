
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

u_to_w = zeros(ny, nx);
w_to_u = zeros(2, nx*ny);
m = 0;

for i = 1:ny
    for j = 1:nx
        if h(i,j) == -Inf
            m=m+1;
            u_to_w(i,j) = m;
            w_to_u(:, m) = [i, j]';
        end
    end
end

%Matrix M building

M0=zeros(m,m);
M1=zeros(m,m);

for m1=1:m
    
    c = w_to_u(:,m1);
    
    if h(c(1)-1,c(2)) < 10000 && h(c(1)-1,c(2)) > -10000
    
        M1(m1,m1)=-1;
        
    elseif h(c(1)-1,c(2)) == -Inf
        
        M1(m1,m1)=-1;
        M0(m1, u_to_w(c(1)-1, c(2)))=1;
        
    end
    
    if h(c(1)+1,c(2)) < 10000 && h(c(1)+1,c(2)) > -10000
    
        M1(m1,m1)=M1(m1,m1)-1;
      
    elseif h(c(1)+1,c(2)) == -Inf
        
        M1(m1,m1)=M1(m1,m1)-1;
        M0(m1, u_to_w(c(1)+1, c(2)))=1;
        
    end
    
    if h(c(1),c(2)-1) < 10000 && h(c(1),c(2)-1) > -10000
    
        M1(m1,m1)=M1(m1,m1)-1;
       
    elseif h(c(1),c(2)-1) == -Inf
        
        M1(m1,m1)=M1(m1,m1)-1;
        M0(m1, u_to_w(c(1), c(2)-1))=1;
        
    end
    
    if h(c(1),c(2)+1) < 10000 && h(c(1),c(2)+1) > -10000
    
        M1(m1,m1)=M1(m1,m1)-1;
        
    elseif h(c(1),c(2)+1) == -Inf
        
        M1(m1,m1)=M1(m1,m1)-1;
        M0(m1, u_to_w(c(1), c(2)+1))=1;
        
    end
    
end

M=M0+M1;

%Calculation of the h values

%NaN changed to 0

%Left nodes

for i = ny-nw-nu+1:ny
    for j = 1:nx1
        
        h(i,j)=0;
        
    end    
end

%Right nodes

for i = ny-nw-nu+1:ny
    for j = nx1+nk+1:nx
        
        h(i,j)=0;
        
    end    
end

%Center nodes

for i = ny-nu+1:ny
    for j = nx1+1:nx1+nk
        
        h(i,j)=0;
        
    end    
end

tol=1d-6;
err=1;
t = waitbar(0,'h is calculating...');

%Calculation of the matrix h

while err > tol
    
    waitbar(tol/err)
    
    for m1=1:m
    
        c = w_to_u(:,m1);
   
        h(c(1), c(2))= -Inf;
    
    end

    b=zeros(m,1);
   
    for m1=1:m
    
        c = w_to_u(:,m1);
    
        if h(c(1)-1,c(2)) < 10000 && h(c(1)-1,c(2)) > -10000
    
            b(m1,1)=-h(c(1)-1,c(2));
        
        end
    
        if h(c(1)+1,c(2)) < 10000 && h(c(1)+1,c(2)) > -10000
    
            b(m1,1)=b(m1,1)-h(c(1)+1,c(2));
        
        end
    
        if h(c(1),c(2)-1) < 10000 && h(c(1),c(2)-1) > -10000
    
            b(m1,1)=b(m1,1)-h(c(1),c(2)-1);
        
        end
    
        if h(c(1),c(2)+1) < 10000 && h(c(1),c(2)+1) > -10000

            b(m1,1)=b(m1,1)-h(c(1),c(2)+1);
        
        end
    
    end
   
    u=M\b;

    for m1=1:m
    
        c = w_to_u(:,m1);
   
        h(c(1), c(2))=u(m1,1);
    
    end

    hkp1=h;
   
    %Boundary conditions

    %Left boundary of the model
    
    for i = 2:ny-nu-nw-1
        
        hkp1(i,1) = 0.25*(h(i+1,1)+2*h(i,2)+h(i-1,1));
        
    end
    
    %Left corner
    
    hkp1(1,1) = 0.5*(h(2,1)+h(1,2));
    
    %Bottom boundary
  
    for j = 2:nx-1
        
        hkp1(1,j) = 0.25*(h(1,j-1)+h(1,j+1)+2*h(2,j));
        
    end
    
    %Right corner

    hkp1(1,nx) = 0.5*(h(1,nx-1)+h(2,nx));
    
    %Right boundary of the model

    for i = 2:ny-nu-nw-1
        
        hkp1(i,nx) = 0.25*(h(i+1,nx)+2*h(i,nx-1)+h(i-1,nx));
        
    end
    
    %Left inner wall of the caisson

    for i = ny-hk-nu+2:ny-nu-1
        
        hkp1(i,nx1+1) = 0.25*(h(i+1,nx1+1)+h(i-1,nx1+1)+2*h(i,nx1+2));
        
    end
    
    %Right inner wall of the caisson
    
    for i = ny-hk-nu+2:ny-nu-1
        
        hkp1(i,nx1+nk) = 0.25*(h(i+1,nx1+nk)+h(i-1,nx1+nk)+2*h(i,nx1+nk-1));
        
    end
    
     %Left outer wall of the caisson
    
    for i = ny-nu-hk+2:ny-nu-nw-1
        
        hkp1(i,nx1) = 0.25*(h(i+1,nx1)+h(i-1,nx1)+2*h(i,nx1-1));
    
    end
        
    %Right outer wall of the caisson
    
    for i = ny-nu-hk+2:ny-nu-nw-1
        
        hkp1(i,nx1+nk+1) = 0.25*(h(i+1,nx1+nk+1)+h(i-1,nx1+nk+1)+2*h(i,nx1+nk+2));
    
    end
    
    err = sqrt(sum(sum((hkp1-h).^2)));
    
    h = hkp1;
    
end

close(t)

for i = ny-nw-nu+1:ny
    for j = 1:nx1
        
        h(i,j)=NaN;
        
    end    
end

%Right nodes

for i = ny-nw-nu+1:ny
    for j = nx1+nk+1:nx
        
        h(i,j)=NaN;
        
    end    
end

%Center nodes

for i = ny-nu+1:ny
    for j = nx1+1:nx1+nk
        
        h(i,j)=NaN;
        
    end    
end

q1=0;

for j=2:nx1-1
    
    q1 = q1 + abs(h(ny-nu-nw-1,j))-abs(h(ny-nu-nw-3,j));
    
end

q = 0.5*(abs(h(ny-nu-nw-1,1))-abs(h(ny-nu-nw-3,1)))+q1+0.5*(abs(h(ny-nu-nw-1,nx1))-abs(h(ny-nu-nw-3,nx1)))

h

%Stream lines

s=ones(ny,nx);

for i = 1:ny
    for j = 1:nx
        
        s(i,j)= -Inf;
        
    end    
end

%Left nodes

for i = ny-nw-nu+1:ny
    for j = 1:nx1
        
        s(i,j)=NaN;
        
    end    
end

%Right nodes

for i = ny-nw-nu+1:ny
    for j = nx1+nk+1:nx
        
        s(i,j)=NaN;
        
    end    
end

%Center nodes

for i = ny-nu+1:ny
    for j = nx1+1:nx1+nk
        
        s(i,j)=NaN;
        
    end    
end

%Pressure inside the caisson

s(ny-nu,nx1+1:nx1+nk)=linspace(0,0,nk);

s(1:ny-nu,nx1+(nk-1)/2+1)=linspace(q,q,ny-nu);

%Pressure outside the caisson

s(ny-nw-nu,1:nx1)=linspace(q,0,nx1);
s(ny-nw-nu,nx1+nk+1:nx)=linspace(0,q,nx1);

%Pressure on the inner boundary of the caisson

s(ny-hk-nu+1:ny-nu,nx1+1)=linspace(0,0,hk);
s(ny-hk-nu+1:ny-nu,nx1+nk)=linspace(0,0,hk);

%Pressure on the outer boundary of the caisson

s(ny-nu-hk+1:ny-nu-nw,nx1)=linspace(0,0,hk-nw);
s(ny-nu-hk+1:ny-nu-nw,nx1+nk+1)=linspace(0,0,hk-nw);

%Pressure on the outer boundary of the soil model

s(1:ny-nu-nw,1)=linspace(q,q,ny-nu-nw);
s(1:ny-nu-nw,nx)=linspace(q,q,ny-nu-nw);

%Pressure on the bottom boundary

s(1,1:nx)=linspace(q,q,nx);

s

s(1,1:nx)=linspace(q,q,nx);

v_to_z = zeros(ny, nx);
z_to_v = zeros(2, nx*ny);
n = 0;

for i = 1:ny
    for j = 1:nx
        if s(i,j) == -Inf
            n=n+1;
            v_to_z(i,j) = n;
            z_to_v(:, n) = [i, j]';
        end
    end
end

%Matrix N created

N0=zeros(n,n);
N1=zeros(n,n);

for n1=1:n
    
    d = z_to_v(:,n1);
    
    if s(d(1)-1,d(2)) < 10000 && s(d(1)-1,d(2)) > -10000
    
        N1(n1,n1)=-1;
        
    elseif s(d(1)-1,d(2)) == -Inf
        
        N1(n1,n1)=-1;
        N0(n1, v_to_z(d(1)-1, d(2)))=1;
        
    end
    
    if s(d(1)+1,d(2)) < 10000 && s(d(1)+1,d(2)) > -10000
    
        N1(n1,n1)=N1(n1,n1)-1;
      
    elseif s(d(1)+1,d(2)) == -Inf
        
        N1(n1,n1)=N1(n1,n1)-1;
        N0(n1, v_to_z(d(1)+1, d(2)))=1;
        
    end
    
    if s(d(1),d(2)-1) < 10000 && s(d(1),d(2)-1) > -10000
    
        N1(n1,n1)=N1(n1,n1)-1;
       
    elseif s(d(1),d(2)-1) == -Inf
        
        N1(n1,n1)=N1(n1,n1)-1;
        N0(n1, v_to_z(d(1), d(2)-1))=1;
        
    end
    
    if s(d(1),d(2)+1) < 10000 && s(d(1),d(2)+1) > -10000
    
        N1(n1,n1)=N1(n1,n1)-1;
        
    elseif s(d(1),d(2)+1) == -Inf
        
        N1(n1,n1)=N1(n1,n1)-1;
        N0(n1, v_to_z(d(1), d(2)+1))=1;
        
    end
    
end

N=N0+N1;

%Calculation of the s values

%NaN changed to 0

%Left nodes

for i = ny-nw-nu+1:ny
    for j = 1:nx1
        
        s(i,j)=0;
        
    end    
end

%Center nodes

for i = ny-nw-nu+1:ny
    for j = nx1+nk+1:nx
        
        s(i,j)=0;
        
    end    
end

%Right nodes

for i = ny-nu+1:ny
    for j = nx1+1:nx1+nk
        
        s(i,j)=0;
        
    end    
end

tol=1d-6;
erro=1;
l = waitbar(0,'s is calculating...');

while erro > tol
    
    waitbar(tol/erro)
    
    for n1=1:n
    
        d = z_to_v(:,n1);
   
        s(d(1), d(2))= -Inf;
    
    end

    f=zeros(n,1);
   
    for n1=1:n
    
        d = z_to_v(:,n1);
    
        if s(d(1)-1,d(2)) < 10000 && s(d(1)-1,d(2)) > -10000
    
            f(n1,1)=-s(d(1)-1,d(2));
        
        end
    
        if s(d(1)+1,d(2)) < 10000 && s(d(1)+1,d(2)) > -10000
    
            f(n1,1)=f(n1,1)-s(d(1)+1,d(2));
        
        end
    
        if s(d(1),d(2)-1) < 10000 && s(d(1),d(2)-1) > -10000
    
            f(n1,1)=f(n1,1)-s(d(1),d(2)-1);
        
        end
    
        if s(d(1),d(2)+1) < 10000 && s(d(1),d(2)+1) > -10000

            f(n1,1)=f(n1,1)-s(d(1),d(2)+1);
        
        end
    
    end
   
    z=N\f;

    for n1=1:n
    
        d = z_to_v(:,n1);
   
        s(d(1), d(2))=z(n1,1);
    
    end

    skp1=s;
   
    %Boundary conditions

    %Upper left boundary of the model
    
    for j = 2:nx1-1
        
        skp1(ny-nu-nw,j) = 0.25*(s(ny-nu-nw,j+1)+2*s(ny-7-nw,j)+s(ny-nu-nw,j-1));
        
    end
    
    %Upper right boundary of the model
    
    for j = nx1+nk+2:nx-1
        
        skp1(ny-nu-nw,j) = 0.25*(s(ny-nu-nw,j+1)+2*s(ny-nu-1-nw,j)+s(ny-nu-nw,j-1));
        
    end
    
    %Inner left side of the model
    
    for j = nx1+2:nx1+(nk-1)/2
        
        skp1(ny-nu,j) = 0.25*(s(ny-nu,j+1)+2*s(ny-nu-1,j)+s(ny-nu,j-1));
        
    end
    
    %Inner right side of the model
    
    for j = nx1+(nk-1)/2+2:nx1+nk-1
        
        skp1(ny-nu,j) = 0.25*(s(ny-nu,j+1)+2*s(ny-nu-1,j)+s(ny-nu,j-1));
        
    end
    
    erro = sqrt(sum(sum((skp1-s).^2)));
    
    s = skp1;
    
end

close(l)

s