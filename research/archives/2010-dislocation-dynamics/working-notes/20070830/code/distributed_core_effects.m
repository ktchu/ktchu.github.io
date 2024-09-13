%
% This script demonstrates that the use of a distributed core when computing
% the stress field around dislocation lines produces errors on the order 
% of the "size" of the distributed core.
% 
% The dislocation configuration is taken to be an array of edge dislocations
% with Burgers vector in the (1,0,0) direction and line direction in the 
% (0,0,1) direction.  The array ordered such that the normal to the plane
% of the dislocation array is in the (1,0,0) direction.  The dislocation
% array is periodically repeated in the x-direction.
%
% The stress fields for this configuration are computed using the Fourier
% transform formulation in Xiang et. al.  The analytical stress fields
% are computed using the formulae from Section 19-5 of Hirth & Lothe.
%
% Kevin T. Chu
% MAE, Princeton University
% 08/2007
%

clear;
format long;
addpath ..


% set grid parameters
Nx = 2400;  % third-order delta function
Ny = 1920;  % third-order delta function
Nx = 100;
Ny = 80;
Nx = 200;
Ny = 160;
Nx = 400;
Ny = 320;
Nx = 800;
Ny = 640;
x_lo = -50;
x_lo = -50;
x_hi =  50;
y_lo = -40;
y_hi =  40;
dx = (x_hi-x_lo)/Nx;
dy = (y_hi-y_lo)/Ny;

% core parameters 
physical_core_radius  = 4;
numerical_core_radius = 2;
delta_function_type = 1;

% dislocation parameters
burgers_vectors = [1, 0, 0];
positions = [0.5, 0];
positions = [0.35, 0.25];
positions = [0, 0];
num_dislocation_arrays = size(burgers_vectors,1);

% elastic constants
G = 1;
poisson_ratio = 1/3;

% generate grid
x = x_lo:dx:x_hi;
if (x(end) == x_hi)
  x = x(1:end-1);
end
y = y_lo:dy:y_hi;
if (y(end) == y_hi)
  y = y(1:end-1);
end
[Y,X] = meshgrid(y,x);


[ sigma_xx_numerical, sigma_yy_numerical, sigma_xy_numerical, d_z] = ...
  compute_stress_fields_numerical(X, Y, ...
                                  num_dislocation_arrays, ...
                                  burgers_vectors, positions, ...
                                  numerical_core_radius, ...
                                  delta_function_type, ...
                                  G, poisson_ratio);

% calculate analytical stress due to the center dislocation array
% and one image on each side
[ sigma_xx_analytical, sigma_yy_analytical, sigma_xy_analytical] = ...
  compute_stress_fields_sharp(X, Y, ...
                              num_dislocation_arrays, ...
                              burgers_vectors, positions, ...
                              G, poisson_ratio);
[ sigma_xx_analytical_plus, ...
  sigma_yy_analytical_plus, ...
  sigma_xy_analytical_plus] = ...
  compute_stress_fields_sharp(X+(x_hi-x_lo), Y, ...
                              num_dislocation_arrays, ...
                              burgers_vectors, positions, ...
                              G, poisson_ratio);
[ sigma_xx_analytical_minus, ...
  sigma_yy_analytical_minus, ...
  sigma_xy_analytical_minus] = ...
  compute_stress_fields_sharp(X-(x_hi-x_lo), Y, ...
                              num_dislocation_arrays, ...
                              burgers_vectors, positions, ...
                              G, poisson_ratio);
sigma_xx_analytical = sigma_xx_analytical ...
                    + sigma_xx_analytical_plus ...
                    + sigma_xx_analytical_minus;
sigma_yy_analytical = sigma_yy_analytical ...
                    + sigma_yy_analytical_plus ...
                    + sigma_yy_analytical_minus;
sigma_xy_analytical = sigma_xy_analytical ...
                    + sigma_xy_analytical_plus ...
                    + sigma_xy_analytical_minus;


% plot last dislocation line field
figure(1), clf;
pcolor(X-0.5*dx,Y-0.5*dy,d_z);
axis([x_lo+0.5*dx x_hi-0.5*dx y_lo+0.5*dy y_hi-0.5*dy]);
shading flat 
colorbar;


% plot stress fields
figure(2), clf;
inside_core = find( sqrt((X-positions(1,1)).^2+(Y-positions(1,2)).^2) ...
                    <= physical_core_radius);
sigma_xx_numerical(inside_core) = 0;
pcolor(X-0.5*dx,Y-0.5*dy,sigma_xx_numerical);
axis([x_lo+0.5*dx x_hi-0.5*dx y_lo+0.5*dy y_hi-0.5*dy]);
shading flat 
title('\sigma_{xx}');
colorbar;

figure(3), clf;
inside_core = find( sqrt((X-positions(1,1)).^2+(Y-positions(1,2)).^2) ...
                    <= physical_core_radius);
sigma_yy_numerical(inside_core) = 0;
pcolor(X-0.5*dx,Y-0.5*dy,sigma_yy_numerical);
axis([x_lo+0.5*dx x_hi-0.5*dx y_lo+0.5*dy y_hi-0.5*dy]);
shading flat 
title('\sigma_{yy}');
colorbar;

figure(4), clf;
inside_core = find( sqrt((X-positions(1,1)).^2+(Y-positions(1,2)).^2) ...
                    <= physical_core_radius);
sigma_xy_numerical(inside_core) = 0;
pcolor(X-0.5*dx,Y-0.5*dy,sigma_xy_numerical);
axis([x_lo+0.5*dx x_hi-0.5*dx y_lo+0.5*dy y_hi-0.5*dy]);
shading flat 
title('\sigma_{xy}');
colorbar;


% compute errors for stress field outside of dislocation cores
outside_cores = find( sqrt((X-positions(1,1)).^2+(Y-positions(1,2)).^2) ...
                    > physical_core_radius);
for line_num = 2:num_dislocation_arrays 
  outside_core_cur = ...
    find( sqrt((X-positions(line_num,1)).^2+(Y-positions(line_num,2)).^2) ...
        > physical_core_radius);
  outside_cores = intersect(outside_cores,outside_core_cur);
end

err_sigma_xx = zeros(size(X));
err_sigma_yy = zeros(size(X));
err_sigma_xy = zeros(size(X));
err_sigma_xx(outside_cores) = abs(sigma_xx_numerical(outside_cores) ...
                                 -sigma_xx_analytical(outside_cores) );
max_err_sigma_xx = max(max(err_sigma_xx))
err_sigma_yy(outside_cores) = abs(sigma_yy_numerical(outside_cores) ...
                                 -sigma_yy_analytical(outside_cores) );
max_err_sigma_yy = max(max(err_sigma_yy))
err_sigma_xy(outside_cores) = abs(sigma_xy_numerical(outside_cores) ...
                                 -sigma_xy_analytical(outside_cores) );
max_err_sigma_xy = max(max(err_sigma_xy))
sum_err_sigma = zeros(size(X));
sum_err_sigma = err_sigma_xx + err_sigma_yy + err_sigma_xy;
max_sum_err_sigma = max(max(sum_err_sigma))


% plot errors
figure(5), clf;
pcolor(X-0.5*dx,Y-0.5*dy,err_sigma_xx);
axis([x_lo+0.5*dx x_hi-0.5*dx y_lo+0.5*dy y_hi-0.5*dy]);
shading flat 
title('Error \sigma_{xx}');
colorbar;

figure(6), clf;
pcolor(X-0.5*dx,Y-0.5*dy,err_sigma_yy);
axis([x_lo+0.5*dx x_hi-0.5*dx y_lo+0.5*dy y_hi-0.5*dy]);
shading flat 
title('Error \sigma_{yy}');
colorbar;

figure(7), clf;
pcolor(X-0.5*dx,Y-0.5*dy,err_sigma_xy);
axis([x_lo+0.5*dx x_hi-0.5*dx y_lo+0.5*dy y_hi-0.5*dy]);
shading flat 
title('Error \sigma_{xy}');
colorbar;


% check error in known relationships between stress components

% sigma_xx = sigma_yy when x = 0
err_1 = zeros(size(y));
idx = find(x==0);
err_1 = sigma_xx_numerical(idx,:)-sigma_yy_numerical(idx,:);
inside_core = find(abs(y) < physical_core_radius);
err_1(inside_core) = 0;
max_err_1 = norm(err_1,inf)
figure(8), clf;
plot(y,err_1);

% sigma_xx = 0 = sigma_yy when y = 0
err_2 = zeros(size(x));
err_3 = zeros(size(x));
idx = find(y==0);
err_2 = sigma_xx_numerical(:,idx);
err_3 = sigma_yy_numerical(:,idx);
inside_core = find(abs(x) < physical_core_radius);
err_2(inside_core) = 0;
err_3(inside_core) = 0;
max_err_2 = norm(err_2,inf)
max_err_3 = norm(err_3,inf)
figure(9), clf;
plot(x,err_2);
hold on;
plot(x,err_3,'r');
