%
% This script summarizes the results from the computational experiments
% on the effects of the numerical core size on the accuracy of the stress
% fields.
%
% For all of the results, the dislocation configuration is taken to be 
% an array of edge dislocations with Burgers vector in the (1,0,0) 
% direction and line direction in the (0,0,1) direction.  The array 
% ordered such that the normal to the plane of the dislocation array is 
% in the (1,0,0) direction.  The dislocation array is periodically 
% repeated in the x-direction.  The simulation cells is [-50,50]x[-40,40],
% which implies that the inter-dislocation spacing is 80 and the 
% distance between periodic images of the dislocation array is 100.
%
% NOTES on numerics:
% - All calculations are on a uniform grid with dx = dy.
% - The errors for the stress fields are only computed OUTSIDE of the 
%   physical core.
%
% Kevin T. Chu
% MAE, Princeton University
% 08/2007
%

clear;
format long;

% max error in sigma as a function of grid resolution
% - err = abs(err_sigma_xx) + abs(err_sigma_yy) + abs(err_sigma_xy)
% - fixed physical core size = 4
% - fixed numerical core size = 1
% - first-order delta function
dx  = [1, 0.5, 0.25, 0.125, 0.0625];
err_vary_dx = [0.047007178457961, ...
               0.002209366222708, ...
               0.002136971107087, ...
               0.002155556102371, ...
               0.002166998068433];
figure(1), clf;
plot(dx,err_vary_dx,'bo-');
axis([0 0.5 0 0.005]);
xlabel('dx = dy');
ylabel('Error');


% max error in sigma outside of physical core with fixed size as a 
% function of the numerical core size
% - err = abs(err_sigma_xx) + abs(err_sigma_yy) + abs(err_sigma_xy)
% - fixed dx = dy = 0.25
% - fixed physical core size = 4
% - first-order delta function
num_core_radius = [0.5, 1, 2, 4, 8];
err_vary_num_core = [5.381040695322851e-04, ...
                     0.002136971107087, ...
                     0.008677219835230, ...
                     0.034509006204891, ...
                     0.082416413759933];
P = polyfit(log(num_core_radius),log(err_vary_num_core),1);
figure(2), clf;
loglog(num_core_radius,err_vary_num_core,'bo');
hold on;
plot(1e-1:1:11, exp(P(1)*log(1e-1:1:11)+P(2)), 'r');
axis([1e-1 10 1e-5 0.1]);
text_string = sprintf('Slope =%5.2f', P(1));
text(0.5,0.01,text_string);
xlabel('Numerical Core Size');
ylabel('Error');


% max error in sigma as a function of the numerical core size with
% physical core equal to numerical core size
% - err = abs(err_sigma_xx) + abs(err_sigma_yy) + abs(err_sigma_xy)
% - fixed dx = dy = 0.25
% - physical core size = numerical core size
% - first-order delta function
num_core_radius = [0.5, 1, 2, 4, 8];
err_vary_core = [0.190322357065602, ...
                 0.123373903604078, ...
                 0.066057254865894, ...
                 0.034509006204891, ...
                 0.017253306456731];
P = polyfit(log(num_core_radius),log(err_vary_core),1);
figure(3), clf;
loglog(num_core_radius,err_vary_core,'bo');
hold on;
%plot(num_core_radius, exp(P(1)*log(num_core_radius)+P(2)), 'r');
%axis([0 8 1e-2 1]);
plot(1e-1:1:11, exp(P(1)*log(1e-1:1:11)+P(2)), 'r');
axis([1e-1 10 1e-2 1]);
text_string = sprintf('Slope =%5.2f', P(1));
text(1,0.5,text_string);
xlabel('Numerical Core Size = Physical Core Size');
ylabel('Error');


% max error in sigma outside of physical core with fixed size as a 
% function of the numerical core size for first-, second-, and 
% third-order smoothed delta functions
% - err = abs(err_sigma_xx) + abs(err_sigma_yy) + abs(err_sigma_xy)
% - fixed dx = dy = 0.125 (first- and second-order)
% - fixed dx = dy = 0.041667 (third-order)
% - fixed physical core size = 4
num_core_radius = [0.5, 1, 2, 4, 8];
err_delta_1 = [5.413797046938378e-04, ...
               0.002155556102371, ...
               0.008677100593490, ...
               0.034508921062433, ...
               0.082416414104439];
err_delta_2 = [2.288983050390715e-05, ...
               2.825462018584898e-05, ...
               4.131105989759552e-04, ...
               0.003452042777642, ...
               0.029121904172413];
num_core_radius_3 = [2, 4, 8];
err_delta_3 = [7.515807169598994e-06, ...
               5.129040433898726e-04, ...
               0.008110110648768];
figure(4), clf;
semilogy(num_core_radius,err_delta_1,'bo-');
hold on;
plot(num_core_radius,err_delta_2,'r*-');
plot(num_core_radius_3,err_delta_3,'m^-');
axis([0 8 1e-6 0.1]);
xlabel('Numerical Core Size');
ylabel('Error');
legend('first-order','second-order','third-order','location','southeast');

