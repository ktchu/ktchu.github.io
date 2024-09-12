%
% This script demonstrates that the use of a smeared core when computing
% the stress field around dislocation lines produces errors on the order 
% of the "size" of the smearing.
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

% set grid parameters
Nx = 2400;  % third-order delta function
Ny = 1920;  % third-order delta function
Nx = 400;
Ny = 320;
x_lo = -50;
x_lo = -50;
x_hi =  50;
y_lo = -40;
y_hi =  40;
dx = (x_hi-x_lo)/Nx;
dy = (y_hi-y_lo)/Ny;

% smearing factor
physical_core_radius  = 3;
numerical_core_radius = 1;
delta_function_order  = 1;

% dislocation parameters
b = [1, 0, 0];
pos = [0, 0];

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
[X,Y] = meshgrid(x,y);


% initialize numerical and analytical solution for stress
sigma_xx_numerical = zeros(size(X));
sigma_yy_numerical = zeros(size(X));
sigma_xy_numerical = zeros(size(X));
sigma_xx_analytical = zeros(size(X));
sigma_yy_analytical = zeros(size(X));
sigma_xy_analytical = zeros(size(X));


% compute stress field by summing up the stress field from 
% individual dislocations
num_dislocations = size(b,1);
for line_num = 1:num_dislocations

  switch (delta_function_order) 

    case 1 
      % compute delta functions of pos using first-order delta function
      delta_X = zeros(size(X));
      x_in_core = find(abs(X - pos(line_num,1)) <= numerical_core_radius); 
      delta_X(x_in_core) = 0.5/numerical_core_radius ...
        * (1+cos(pi*(X(x_in_core)-pos(line_num,1))/numerical_core_radius));

      delta_Y = zeros(size(Y));
      y_in_core = find(abs(Y - pos(line_num,2)) <= numerical_core_radius); 
      delta_Y(y_in_core) = 0.5/numerical_core_radius ...
        * (1+cos(pi*(Y(y_in_core)-pos(line_num,2))/numerical_core_radius));

    case 2
      % compute delta functions of pos using second-order delta function
 
     delta_X = zeros(size(X));
     x_in_core = find(abs(X - pos(line_num,1)) <= numerical_core_radius); 
     delta_X(x_in_core) = -1/6/numerical_core_radius ...
        * (1+cos(pi*(X(x_in_core)-pos(line_num,1))/numerical_core_radius));
      x_in_core = find(abs(X - pos(line_num,1)) <= 0.5*numerical_core_radius); 
      delta_X(x_in_core) = 4/3/numerical_core_radius ...
        * (1+cos(2*pi*(X(x_in_core)-pos(line_num,1))/numerical_core_radius))...
                          -1/6/numerical_core_radius ...
        * (1+cos(pi*(X(x_in_core)-pos(line_num,1))/numerical_core_radius));

      delta_Y = zeros(size(Y));
      y_in_core = find(abs(Y - pos(line_num,1)) <= numerical_core_radius); 
      delta_Y(y_in_core) = -1/6/numerical_core_radius ...
        * (1+cos(pi*(Y(y_in_core)-pos(line_num,2))/numerical_core_radius));
      y_in_core = find(abs(Y - pos(line_num,2)) <= 0.5*numerical_core_radius); 
      delta_Y(y_in_core) = 4/3/numerical_core_radius ...
        * (1+cos(2*pi*(Y(y_in_core)-pos(line_num,2))/numerical_core_radius))...
                          -1/6/numerical_core_radius ...
        * (1+cos(pi*(Y(y_in_core)-pos(line_num,2))/numerical_core_radius));

    case 3
      % compute delta functions of pos using third-order delta function
 
     delta_X = zeros(size(X));
     x_in_core = find(abs(X - pos(line_num,1)) <= numerical_core_radius); 
     delta_X(x_in_core) = 1/48/numerical_core_radius ...
        * (1+cos(pi*(X(x_in_core)-pos(line_num,1))/numerical_core_radius));
      x_in_core = find(abs(X - pos(line_num,1)) <= 0.5*numerical_core_radius); 
      delta_X(x_in_core) = -16/15/numerical_core_radius ...
        * (1+cos(2*pi*(X(x_in_core)-pos(line_num,1))/numerical_core_radius))...
                          +1/48/numerical_core_radius ...
        * (1+cos(pi*(X(x_in_core)-pos(line_num,1))/numerical_core_radius));
      x_in_core = find(abs(X - pos(line_num,1)) <= numerical_core_radius/3); 
      delta_X(x_in_core) = 243/80/numerical_core_radius ...
        * (1+cos(3*pi*(X(x_in_core)-pos(line_num,1))/numerical_core_radius))...
                           -16/15/numerical_core_radius ...
        * (1+cos(2*pi*(X(x_in_core)-pos(line_num,1))/numerical_core_radius))...
                          +1/48/numerical_core_radius ...
        * (1+cos(pi*(X(x_in_core)-pos(line_num,1))/numerical_core_radius));

     delta_Y = zeros(size(Y));
     y_in_core = find(abs(Y - pos(line_num,1)) <= numerical_core_radius); 
     delta_Y(y_in_core) = 1/48/numerical_core_radius ...
        * (1+cos(pi*(Y(y_in_core)-pos(line_num,1))/numerical_core_radius));
      y_in_core = find(abs(Y - pos(line_num,1)) <= 0.5*numerical_core_radius); 
      delta_Y(y_in_core) = -16/15/numerical_core_radius ...
        * (1+cos(2*pi*(Y(y_in_core)-pos(line_num,1))/numerical_core_radius))...
                          +1/48/numerical_core_radius ...
        * (1+cos(pi*(Y(y_in_core)-pos(line_num,1))/numerical_core_radius));
      y_in_core = find(abs(Y - pos(line_num,1)) <= numerical_core_radius/3); 
      delta_Y(y_in_core) = 243/80/numerical_core_radius ...
        * (1+cos(3*pi*(Y(y_in_core)-pos(line_num,1))/numerical_core_radius))...
                           -16/15/numerical_core_radius ...
        * (1+cos(2*pi*(Y(y_in_core)-pos(line_num,1))/numerical_core_radius))...
                          +1/48/numerical_core_radius ...
        * (1+cos(pi*(Y(y_in_core)-pos(line_num,1))/numerical_core_radius));


    otherwise
      error('Unsupported order for delta function');
  end

  % compute dislocation line field
  d_z = delta_X.*delta_Y;

  % compute FFT of d_z
  d_z_fft = fft2(d_z);

  % compute stress field from current dislocation in Fourier space
  k_x = 0:length(x)-1;
  k_x(ceil((length(x)+1)/2)+1:end) = k_x(ceil((length(x)+1)/2)+1:end) - length(x); 
  k_x = 2*pi*k_x/(x_hi-x_lo);
  k_y = 0:length(y)-1;
  k_y(ceil((length(y)+1)/2):end) = k_y(ceil((length(y)+1)/2):end) - length(y); 
  k_y = 2*pi*k_y/(y_hi-y_lo);
  [K_x,K_y] = meshgrid(k_x, k_y);
  norm_K_sq = K_x.^2 + K_y.^2;
  sigma_xx_fft = -2*G*i/(1-poisson_ratio)*K_y.^2./norm_K_sq.^2 ...
               .*(K_x*b(line_num,2) - K_y*b(line_num,1)).*d_z_fft;
  sigma_yy_fft = -2*G*i/(1-poisson_ratio)*K_x.^2./norm_K_sq.^2 ...
               .*(K_x*b(line_num,2) - K_y*b(line_num,1)).*d_z_fft;
  sigma_xy_fft = 2*G*i/(1-poisson_ratio)*K_x.*K_y./norm_K_sq.^2 ...
               .*(K_x*b(line_num,2) - K_y*b(line_num,1)).*d_z_fft;
  
  % correct DC component of stress fields in Fourier space
  sigma_xx_fft(1,1) = 0;
  sigma_yy_fft(1,1) = 0;
  sigma_xy_fft(1,1) = 0;

  % compute stress contribution in real space
  sigma_xx_cur = ifft2(sigma_xx_fft);
  sigma_yy_cur = ifft2(sigma_yy_fft);
  sigma_xy_cur = ifft2(sigma_xy_fft);

  % check that imaginary part is small
  if ( norm(abs(imag(sigma_xx_cur))) > 1e-10 ...
     | norm(abs(imag(sigma_yy_cur))) > 1e-10 ...
     | norm(abs(imag(sigma_xy_cur))) > 1e-10 )
     norm(abs(imag(sigma_xx_cur)))
     norm(abs(imag(sigma_yy_cur)))
     norm(abs(imag(sigma_xy_cur)))
    error('Imaginary component of stress nonzero!!');
  end

  sigma_xx_numerical = sigma_xx_numerical + real(sigma_xx_cur);
  sigma_yy_numerical = sigma_yy_numerical + real(sigma_yy_cur);
  sigma_xy_numerical = sigma_xy_numerical + real(sigma_xy_cur);

  % compute contribution to analytical solution for stress 
  % including all periodic images in the y-direction and
  % first neighbors in the x-direction
  [sigma_xx_cur, sigma_xy_cur, sigma_yy_cur] = ...
    sigma_D( b(line_num,:), pos(line_num,:), ...
             (y_hi-y_lo), x, y );
  sigma_xx_cur = 0.5*G/(1-poisson_ratio)*sigma_xx_cur; 
  sigma_yy_cur = 0.5*G/(1-poisson_ratio)*sigma_yy_cur; 
  sigma_xy_cur = 0.5*G/(1-poisson_ratio)*sigma_xy_cur; 
  sigma_xx_analytical = sigma_xx_analytical + sigma_xx_cur;
  sigma_yy_analytical = sigma_yy_analytical + sigma_yy_cur;
  sigma_xy_analytical = sigma_xy_analytical + sigma_xy_cur;

  pos_neighbor = pos(line_num,:);
  pos_neighbor(1) = pos_neighbor(1) + (x_hi-x_lo);
  [sigma_xx_cur, sigma_xy_cur, sigma_yy_cur] = ...
    sigma_D( b(line_num,:), pos_neighbor, ...
             (y_hi-y_lo), x, y );
  sigma_xx_cur = 0.5*G/(1-poisson_ratio)*sigma_xx_cur; 
  sigma_yy_cur = 0.5*G/(1-poisson_ratio)*sigma_yy_cur; 
  sigma_xy_cur = 0.5*G/(1-poisson_ratio)*sigma_xy_cur; 
  sigma_xx_analytical = sigma_xx_analytical + sigma_xx_cur;
  sigma_yy_analytical = sigma_yy_analytical + sigma_yy_cur;
  sigma_xy_analytical = sigma_xy_analytical + sigma_xy_cur;

  pos_neighbor = pos(line_num,:);
  pos_neighbor(1) = pos_neighbor(1) - (x_hi-x_lo);
  [sigma_xx_cur, sigma_xy_cur, sigma_yy_cur] = ...
    sigma_D( b(line_num,:), pos_neighbor, ...
             (y_hi-y_lo), x, y );
  sigma_xx_cur = 0.5*G/(1-poisson_ratio)*sigma_xx_cur; 
  sigma_yy_cur = 0.5*G/(1-poisson_ratio)*sigma_yy_cur; 
  sigma_xy_cur = 0.5*G/(1-poisson_ratio)*sigma_xy_cur; 
  sigma_xx_analytical = sigma_xx_analytical + sigma_xx_cur;
  sigma_yy_analytical = sigma_yy_analytical + sigma_yy_cur;
  sigma_xy_analytical = sigma_xy_analytical + sigma_xy_cur;

end


% plot last dislocation line field
figure(1), clf;
pcolor(X-0.5*dx,Y-0.5*dy,d_z);
axis([x_lo+0.5*dx x_hi-0.5*dx y_lo+0.5*dy y_hi-0.5*dy]);
shading flat 
colorbar;


% plot stress fields
figure(2), clf;
pcolor(X-0.5*dx,Y-0.5*dy,sigma_xx_numerical);
axis([x_lo+0.5*dx x_hi-0.5*dx y_lo+0.5*dy y_hi-0.5*dy]);
shading flat 
title('\sigma_{xx}');
colorbar;

figure(3), clf;
pcolor(X-0.5*dx,Y-0.5*dy,sigma_yy_numerical);
axis([x_lo+0.5*dx x_hi-0.5*dx y_lo+0.5*dy y_hi-0.5*dy]);
shading flat 
title('\sigma_{yy}');
colorbar;

figure(4), clf;
pcolor(X-0.5*dx,Y-0.5*dy,sigma_xy_numerical);
axis([x_lo+0.5*dx x_hi-0.5*dx y_lo+0.5*dy y_hi-0.5*dy]);
shading flat 
title('\sigma_{xy}');
colorbar;


% compute errors for stress field outside of dislocation cores
outside_cores = find( sqrt((X-pos(1,1)).^2+(Y-pos(1,2)).^2) ...
                    > physical_core_radius);
for line_num = 2:num_dislocations
  outside_core_cur = ...
    find( sqrt((X-pos(line_num,1)).^2+(Y-pos(line_num,2)).^2) ...
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
err_1 = sigma_xx_numerical(:,idx)-sigma_yy_numerical(:,idx);
inside_core = find(abs(y) < physical_core_radius);
err_1(inside_core) = 0;
max_err_1 = norm(err_1,inf)
figure(8), clf;
plot(y,err_1);

