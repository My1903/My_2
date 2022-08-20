clear all
close all
clc
tic
% 
n =20;

Nx = n;
Ny = n;
Nt = (8*n);

Lx = 1; %end point
Ly = 1;
T = 10;

% Mesh
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
t = linspace(0, T, Nt);

%add mesh 
x_ = linspace(0, Lx, 2*Nx-1);
y_ = linspace(0, Ly, 2*Ny-1);

% Mesh size
hx = Lx/(Nx-1);  %% Mesh size (x)
hy = Ly/(Ny-1);  %% Mesh size (y)
dt = T/(Nt-1);
alpha_x= dt/(2*(hx)^2);
alpha_y= dt/(2*(hy)^2);
%vector of index
% T = zeros(n);
% T(1,1:n) = 10; %TOP
% T(n,1:n) = 1;  %BOTTOM
% T(1:n,1) = 1;  %LEFT
% T(1:n,n) = 1;  %RIGHT
% dt = dx^2/4;

% number of unknown at one step
N = Nx * Ny;

M = zeros(N, N);   % N rows, N columns
B = zeros(N, 1);   % N rows, 1 columns
U_0 = zeros(N, 1); % N rows, 1 column
% D_x=zeros(Nx, Ny);
% D_y=zeros(Nx, Ny);
 D_=zeros(2*Nx-1, 2*Ny-1);
disp(D(x_(1),y_(1)));

%meshgrid
[X,Y] = meshgrid(x,y);



% 
% for i=1:Nx-1
%     for j=1:Ny
%         D_x(i,j) = D((x(i)+x(i+1))/2,y(j));
%     end
% end
% 
% for i=1:Nx
%     for j=1:Ny-1
%         D_y(i,j) = D(x(i),(y(j)+y(j+1))/2);
%     end
% end


for i=1:2*Nx-1
    for j=1:2*Ny-1
        D_(i,j) = D(x_(i),y_(j));
    end
end


U_0 = reshape(u_0new(X, Y), [], 1);

%% Set up cofficients of matrix M and B and solve M*x = B
% % loop over t-direction
for k = 1:Nt
   % interior points
      for i = 2:(Nx-1)
        for j = 2:(Ny-1)
            % convert the cell (i,j) to the nth grid point in order to establish M
            n = i + (j-1)*Nx;
            M(n, n)    = alpha_x*(D_(2*i-1-1,j) +D_(2*i-1+1,j))+alpha_y*(D_(i,2*j-1-1)+ D_(i,2*j-1+1)) + 1;
            M(n, n-1)  = -alpha_x*D_(2*i-1-1,j);
            M(n, n+1)  = -alpha_x*D_(2*i-1+1,j);
            M(n, n-Nx) = -alpha_y*D_(i,2*j-1-1);
            M(n, n+Nx) = -alpha_y*D_(i,2*j-1+1);
        end
      end
    % B for interior points
    for j = 1:Ny
        for i = 1:Nx
            n=i+(j-1)*Nx;
            B(n) = -dt*f(U_0(n)) + U_0(n);
        end
    end

   %Boundary conditions
    % BC:{u^n+1}_0j={u^n+1}_2j
    i=1;
    for j=2:(Ny-1)
        n=i+(j-1)*Nx;
        M(n, n)    = alpha_x*(D_(2*i-1,j) +D_(2*i-1+1,j))+alpha_y*(D_(i,2*j-1-1)+ D_(i,2*i-1+1)) + 1; %alpha_x*(D_x(i-1,j)
        M(n, n+1)  = -alpha_x*D_(2*i-1,j) -alpha_x*D_(2*i-1+1,j);  %-alpha_x*D(i-1,j) 
        M(n, n-Nx) = -alpha_y*D_(i,2*j-1-1);
        M(n, n+Nx) = -alpha_y*D_(i,2*j-1+1); 

    end
    % BC:{u^n+1}_Nx+1j={u^n+1}_(Nx-1)j
    i=Nx;
    for j=2:(Ny-1)
        n=i+(j-1)*Nx;
         M(n, n)    = alpha_x*(D_(2*i-1-1,j) +D_(2*i-1,j))+alpha_y*(D_(i,2*j-1-1)+ D_(i,2*j-1+1)) + 1;%D_x(i+1,j))
         M(n, n-1)  = -alpha_x*D_(2*i-1-1,j) -alpha_x*D_(2*i-1,j);  %-alpha_x*D(i+1,j) 
         M(n, n-Nx) = -alpha_y*D_(i,2*j-1-1);
         M(n, n+Nx) = -alpha_y*D_(i,2*j-1+1); 

        
    end
    % BC:{u^n+1}_i0={u^n+1}_i2
    j=1;
    for i=2:Nx-1
        n=i+(j-1)*Nx;
        M(n, n)    = alpha_x*(D_(2*i-1-1,j) +D_(2*i-1+1,j))+alpha_y*(D_(i,2*j-1)+ D_(i,2*j-1+1)) + 1;%alpha_y*(D_y(i,j-1)
        M(n, n+1)  =  -alpha_x*D_(2*i-1+1,j);
        M(n, n-1)  = -alpha_x*D_(2*i-1-1,j) ;
        M(n, n+Nx) = -alpha_y*D_(i,2*j-1)-alpha_y*D_(i,2*j-1+1); %-alpha_y*D(i,j-1)
    end
    % BC:{u^n+1}_iNy+1={u^n+1}_iNy-1
    j=Ny;
    for i=2:Nx-1
        n=i+(j-1)*Nx;
        M(n, n)    = alpha_x*(D_(2*i-1-1,j) +D_(2*i-1+1,j))+alpha_y*(D_(i,2*j-1-1)+ D(i,2*j-1)) + 1;% D_y(i,j+1))
        M(n, n+1)  =  -alpha_x*D_(2*i-1+1,j);
        M(n, n-1)  = -alpha_x*D_(2*i-1-1,j) ;
        M(n, n-Nx) = -alpha_y*D_(i,2*j-1-1)-alpha_y*D_(i,2*j-1); %-alpha_y*D(i,j+1)
     
    end 
     
       j=Ny; i=1;
        n=i+(j-1)*Nx;
        M(n, n)    = alpha_x*(D_(2*i-1,j) +D_(2*i-1+1,j))+alpha_y*(D_(i,2*j-1-1)+ D_(i,2*j-1)) + 1;% D_y(i,j+1)) ,D_x(i-1,j)
        M(n, n+1)  =  -alpha_x*D_(2*i-1+1,j)-alpha_x*D_(2*i-1,j);%-alpha_x*D(i-1,j)
        M(n, n-Nx) = -alpha_y*D_(i,2*j-1-1)-alpha_y*D_(i,2*j-1); %-alpha_y*D(i,j+1)

        j=1; i=Nx;
         n=i+(j-1)*Nx;
         M(n, n)    = alpha_x*(D_(2*i-1-1,j) +D_(2*i-1,j))+alpha_y*(D_(i,2*j-1)+ D_(i,2*j-1+1)) + 1;%(D_y(i,j-1),D_x(i+1,j)
         M(n, n-1)  = -alpha_x*D_(2*i-1-1,j) -alpha_x*D_(2*i-1,j);  %-alpha_x*D(i+1,j) 
         M(n, n+Nx) = -alpha_y*D_(i,2*j-1+1)-alpha_y*D_(i,2*j-1); %-alpha_y*D(i,j-1)


       i=1;j=1;
       n=i+(j-1)*Nx;
        M(n, n)    = alpha_x*(D_(2*i-1,j) +D_(2*i-1+1,j))+alpha_y*(D_(i,2*j-1)+ D(i,2*j-1+1)) + 1;%D_x(i-1,j),D_y(i,j-1)
        M(n, n+1)  =  -alpha_x*D_(2*i-1+1,j)-alpha_x*D_(2*i-1,j);%-alpha_x*D(i-1,j)
        M(n, n+Nx) = -alpha_y*D_(i,2*j-1)-alpha_y*D_(i,2*j-1+1); %-alpha_y*D(i,j-1)
       i=Nx; j=Ny;
         n=i+(j-1)*Nx;
         M(n, n)    = alpha_x*(D_(2*i-1-1,j) +D_(2*i-1,j))+alpha_y*(D_(i,2*j-1-1)+ D_(i,2*j-1)) + 1;%D_x(i+1,j)), D_y(i,j+1)
         M(n, n-1)  = -alpha_x*D_(2*i-1-1,j) -alpha_x*D_(2*i-1,j);  %-alpha_x*D(i+1,j) 
         M(n, n-Nx) = -alpha_y*D_(i,2*j-1-1)-alpha_y*D_(i,2*j-1); %-alpha_y*D(i,j+1)
    

    %% Mx= B ===> x = B\M
    %U_1_vec= inv(M)*B;
    U_1_vec= B\M;
    U_1 = reshape(U_1_vec, Nx, Ny);
   



    % reset the value of U_0   
    U_0 = U_1_vec;
    %%
    
    
% print step of calculation to follow time   
fprintf('Step %d \n', k);
 end
 toc


surf(x,y,U_1);
 shading interp

xlabel('x');
ylabel('y');
