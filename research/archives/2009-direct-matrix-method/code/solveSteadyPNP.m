%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solve the following nonlinear integro-differential equation for 
% steady-state electrochemical thin-films using Newton's method with 
% a pseudospectral discretization based on the Chebyshev grid on [-1,1]:
%
%   epsilon^2 (E'' - 0.5*E^3) - 0.25*(c0(E) + j*(x+1)) E = 0.25*j
%
% where epsilon is the ratio of the Debye screening length to the
% cell size (i.e. = lambda_D/L) and c0(E) is given by
%
%   c0 = (1 - j) + epsilon^2 ( 2*E(1) - 2*E(-1) - int_{-1}^1{E^2 dx} ).
%
% The boundary conditions for this equation are the Butler-Volmer
% rate equations:
%
%  k_c (c(-1) + rho(-1)) - j_r = j
%
%  -k_c (c(1) + rho(1)) + j_r = j
%
% where c(x) = c0(E) + j*(x+1) + 2*epsilon^2*E(x)^2 and
%       rho(x) = 4*epsilon^2*E'(x)
%
% NOTES:  
%  - Relies on cheb() and clencurt() by N. Trefethen available at
%    http://web.comlab.ox.ac.uk/oucl/work/nick.trefethen/spectral.html
%
% REFERENCES:  
% - M. Z. Bazant, K. T. Chu, and B. J. Bayly. "Current-voltage relations for 
%   electrochemical thin films." SIAM J. Appl. Math., 65, 1463–1484 (2005).
% - K. T. Chu and M. Z. Bazant. "Electrochemical thin films at and above the 
%   classical limiting current." SIAM J. Appl. Math., 65, 1485–1505 (2005).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j         = 1.5;  
epsilon   = 0.01; 
k_c       = 10;    
j_r       = 10;    
N         = 200; 
res_tol   = 1e-8;
max_iters = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute grid, differentiation matrix, and quadrature weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D,x] = cheb(N-1);         % Chebyshev differentiation matrix 
L     = D*D;               % Laplacian operator
[x,w] = clencurt(N-1);     % Clenshaw-Curtis quadrature weights 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up continuation in j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_count = 0;  % total number of iterations required 
j_start     = 0.5;
dj          = 0.1;
j_cur       = j_start;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate initial iterate for Newton iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1-j_cur;
c = c0 + j_cur*(x+1);
E = -2*j_cur./(j_cur*(x+1)+c0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton iteration with simple continuation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ( j_cur <= j & dj > 0 )  

  % display j_cur and total count
  j_cur       = j_cur   
  total_count = total_count

  % initialize Newton iteration
  count     = 0;  % iteration count for single Newton iteration
  c0        = 1-j_cur + epsilon^2*(2*E(1)-2*E(N)-w*(E.^2)); 
  res       = epsilon^2*(L*E-0.5*E.^3) ...
            - 0.25*(c0+j_cur*(x+1)).*E - 0.25*j_cur;
  res(1)    = -k_c*(c0+2*j_cur+epsilon^2*(2*E(1)^2+4*D(1,:)*E)) ...
            + j_r - j_cur;
  res(N)    = k_c*(c0+epsilon^2*(2*E(N)^2+4*D(N,:)*E)) ...
            - j_r - j_cur;
  res_norm  = norm(res,inf);

  while ( (res_norm > res_tol) & (count < max_iters) )
  
    % construct Jacobian for interior grid points
    dc0_dE = -2*epsilon^2*w.*E'; 
    dc0_dE(1) = dc0_dE(1) + 2*epsilon^2;
    dc0_dE(N) = dc0_dE(N) - 2*epsilon^2;
    J = epsilon^2*L - diag(1.5*epsilon^2*(E.*E) + 0.25*(c0+j_cur*(x+1))) ...
      - 0.25*kron(E,dc0_dE);
    
    % construct Jacobian for boundary conditions
    J(1,:) = -k_c*(dc0_dE + 4*epsilon^2*D(1,:)); 
    J(1,1) = J(1,1) - 4*k_c*epsilon^2*E(1);
    J(N,:) = k_c*(dc0_dE + 4*epsilon^2*D(N,:)); 
    J(N,N) = J(N,N) + 4*k_c*epsilon^2*E(N);
 
    % compute delta_E and update solution
    delta_E = -J\res;  E = E + delta_E;

    % update residual
    c0     = 1-j_cur + epsilon^2*(2*E(1)-2*E(N)-w*(E.^2)); 
    res    = epsilon^2*(L*E-0.5*E.^3) ...
           - 0.25*(c0+j_cur*(x+1)).*E - 0.25*j_cur;
    res(1) = -k_c*(c0+2*j_cur+epsilon^2*(2*E(1)^2+4*D(1,:)*E)) ...
           + j_r - j_cur;
    res(N) = k_c*(c0+epsilon^2*(2*E(N)^2+4*D(N,:)*E)) ...
           - j_r - j_cur;

    % update loop variables
    res_norm  = norm(res,inf)
    count = count + 1
  
  end  % Newton iteration loop

  % update continuation variables
  if (j - j_cur < dj)
    dj = j - j_cur;
  end
  j_cur = j_cur + dj;
  total_count = total_count + count;

end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
plot(x,E,'k-');
axis([-1 1 -100 0]);
xlabel('x'); ylabel('E','Rotation',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spectral coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coefs = abs(fft([E; flipud(E(2:end-1))]));
figure(2); clf;
semilogy(coefs(1:N),'ko');
axis([0 N 1e-15 1e4]);
xlabel('n'); ylabel('|a_n|','Rotation',0,'Position',[-32 5e-6]);
