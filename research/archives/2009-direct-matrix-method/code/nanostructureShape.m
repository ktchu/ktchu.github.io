%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Solve the following boundary-value problem for nanostructure 
% evolution using Newton's method with a pseudospectral discretization 
% based on the rational Chebyshev functions for a semi-infinite domain:
%
%   (d/dx)^3 f - sqrt(f) + 1 = 0
%
% with boundary conditions
%
%  (1) f(0) = 0
%  (2) df/dx(0) = alpha
%  (3) f --> 1 as x --> infinity
%
% The Neumann boundary condition is imposed by replacing the discretized
% differential equation at the second largest grid point. 
%
% NOTES:  
% - Relies on cheb() by N. Trefethen available at
%   http://web.comlab.ox.ac.uk/oucl/work/nick.trefethen/spectral.html
% 
% REFERENCES:
% - D. Margetis, P.-W. Fok, M.J. Aziz and H.A. Stone, "Continuum theory of 
%   nanostructure decay via a microscale condition." Physical Review Letters 
%   97, 096102 (2006).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
format long
set(0,'DefaultAxesFontSize',18,'DefaultAxesFontName','Helvetica')
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultTextFontSize',24,'DefaultTextFontName','Helvetica')

% set print format
print_format = 'eps';
fig_dir = 'figures';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha     = 1;
f_inf     = 1;    % f at x = infinity
N         = 250;
L         = 5;
res_tol   = 1e-8;
max_iters = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute differentiation matrix and grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute Chebyshev differentiation matrix on interval [-1,1]
[D_y,y] = cheb(N);
one_minus_y = diag(1-y);
warning off MATLAB:divideByZero  % disable warning message when computing x
D = 0.5/L*(one_minus_y^2)*D_y; x = L*(1+y)./(1-y);
warning on MATLAB:divideByZero   % re-enable warning message
x_int = x(2:end-1);  % extract interior grid points 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct discrete operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D3 = D^3;
% contributions to D^3 from interior grid points and grid
% point at infinity
D3_int = D3(2:end-1,2:end-1); D3_inf = D3(2:end-1,1);           
% contribution to Neuman BC from interior grid points and
% grid point at infinity
BC_f_prime_int = D(end,2:end-1);  BC_f_prime_inf = D(end,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize Newton iteration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_int      = f_inf*atan(x_int);  % initial iterate
count      = 0; 
res        = D3_int*f_int + D3_inf*f_inf - sqrt(f_int) + 1;
res(1)     = BC_f_prime_int*f_int + BC_f_prime_inf*f_inf - alpha;  
res_norm   = norm(res,inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ( res_norm > res_tol & count < max_iters) 

  % compute Jacobian
  J = D3_int - 0.5*diag(1./sqrt(f_int));
  J(1,:) = BC_f_prime_int;  % Neumann boundary condition

  % compute delta_f and update solution
  delta_f = -J\res; f_int = f_int + delta_f;

  % update residual
  res    = D3_int*f_int + D3_inf*f_inf - sqrt(f_int) + 1;
  res(1) = BC_f_prime_int*f_int + BC_f_prime_inf*f_inf - alpha; 
 
  % update loop variables
  res_norm = norm(res,inf)
  count = count + 1

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
plot(x,[f_inf; f_int; 0],'k-');
axis([0 25 0 1.5]);
xlabel('x'); ylabel('f','Rotation',0);
filename = sprintf('nanostructure_soln.%s', print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spectral coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = [f_inf; f_int; 0];
coefs = abs(fft([f; flipud(f(2:end-1))]));

figure(2); clf;
loglog(coefs(1:N),'ko');
axis([1 N 1e-12 1000]);
xlabel('n'); ylabel('|a_n|','Rotation',0);
filename = sprintf('nanostructure_coefs_loglog.%s', print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(3); clf;
semilogy(coefs(1:N),'ko');
axis([0 N 1e-12 1000]);
xlabel('n'); ylabel('|a_n|','Rotation',0);
filename = sprintf('nanostructure_coefs_semilogy.%s', print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);
