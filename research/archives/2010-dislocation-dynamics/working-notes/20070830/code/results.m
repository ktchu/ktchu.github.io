%
% This script summarizes the results from a comparison of the 
% nonsingular dislocation theory of Cai with numerical solutions
% and the singular theory of dislocation lines.
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
addpath ..

% max error in stress field outside of physical core as a function of the
% core size parameter when solving for the stress fields numerically or
% using Cai's nonsingular expressions for the stress
% 
% - err = abs(err_sigma_xx) + abs(err_sigma_yy) + abs(err_sigma_xy)
% - fixed physical core size = 4
% - fixed dx = dy = 0.25 
% - first-order delta function for numerical solutions
% 
num_core_radius = [0.5, 1, 2, 4, 8];
err_vary_num_core = [5.381040695322851e-04, ...
                     0.002136971107087, ...
                     0.008677219835230, ...
                     0.034509006204891, ...
                     0.082416413759933];
P_num = polyfit(log(num_core_radius),log(err_vary_num_core),1);
cai_core_radius = [0.25, 0.5, 1, 2, 4, 8];
err_vary_cai_core = [5.039268387475697e-04, ...
                     0.001990053048359, ...
                     0.007574320656547, ...
                     0.025366212459843, ...
                     0.061323945282872, ...
                     0.094685517698294];
P_cai = polyfit(log(cai_core_radius),log(err_vary_cai_core),1);
figure(1), clf;
loglog(num_core_radius,err_vary_num_core,'bo');
hold on;
plot(1e-1:1:11, exp(P_num(1)*log(1e-1:1:11)+P_num(2)), 'b');
loglog(cai_core_radius,err_vary_cai_core,'ro');
plot(1e-1:1:11, exp(P_cai(1)*log(1e-1:1:11)+P_cai(2)), 'r');
axis([1e-1 10 1e-5 0.5]);
text_string = sprintf('Slope (numerical) =%5.2f', P_num(1));
text(1,1e-4,text_string);
text_string = sprintf('Slope (Cai) =%5.2f', P_cai(1));
text(0.15,0.02,text_string);
xlabel('Core Radius');
ylabel('Error');


% max error in sigma as a function of the core size with physical core 
% equal to core size 
% - err = abs(err_sigma_xx) + abs(err_sigma_yy) + abs(err_sigma_xy)
% - fixed dx = dy = 0.25
% - physical core size = core size parameter
% - first-order delta function for numerical solution
num_core_radius = [0.5, 1, 2, 4, 8];
err_vary_core_num = [0.190322357065602, ...
                 0.123373903604078, ...
                 0.066057254865894, ...
                 0.034509006204891, ...
                 0.017253306456731];
cai_core_radius = [0.5, 1, 2, 4, 8];
P_num = polyfit(log(num_core_radius),log(err_vary_core_num),1);
err_vary_core_cai = [0.392346407857646, ...
                     0.216812501122843, ...
                     0.116075478562721, ...
                     0.061323945282872, ...
                     0.030657457566354];
P_cai = polyfit(log(cai_core_radius),log(err_vary_core_cai),1);
figure(2), clf;
loglog(num_core_radius,err_vary_core_num,'bo');
hold on;
plot(1e-1:1:11, exp(P_num(1)*log(1e-1:1:11)+P_num(2)), 'b');
loglog(cai_core_radius,err_vary_core_cai,'ro');
plot(1e-1:1:11, exp(P_cai(1)*log(1e-1:1:11)+P_cai(2)), 'r');
axis([1e-1 10 1e-2 1]);
text_string = sprintf('Slope (numerical) =%5.2f', P_num(1));
text(0.2,0.05,text_string);
text_string = sprintf('Slope (Cai) =%5.2f', P_cai(1));
text(1,0.5,text_string);
xlabel('Core Radius = Physical Core Radius');
ylabel('Error');


% comparison of numerical solution with analytical solution computed
% using Cai's nonsingular expressions for the stress field for various
% ratios of numerical core radius to Cai's core radius
%
% - err = abs(err_sigma_xx) + abs(err_sigma_yy) + abs(err_sigma_xy)
% - fixed dx = dy = 0.25
% - first-order delta function for numerical solution
% - fixed physical core size = 4
core_radius_equal = [0.5, 1, 2, 4];
err_num_vs_cai_equal = [0.001483569447328, ...
                        0.005459014031824, ...
                        0.016831146233217, ...
                        0.031238710306055];
cai_core_radius         = [0.5, 1, 2, 4];
err_num_vs_cai_cai_half = [1.712158677345597e-04, ...
                           0.001480658586741, ...
                           0.011937052095883, ...
                           0.021333026014862];
num_core_radius          = [0.5, 1, 2, 4];
err_num_vs_cai_num_half = [0.007067837055516, ...
                           0.023250905835120, ...
                           0.052788879056246, ...
                           0.075739616982004];
figure(3), clf;
semilogy(core_radius_equal,err_num_vs_cai_equal,'bo-');
hold on;
plot(cai_core_radius,err_num_vs_cai_cai_half,'ro-');
plot(num_core_radius,err_num_vs_cai_num_half,'go-');
axis([0 4 1e-4 0.1]);
xlabel('Core Radius');
ylabel('Max Difference');
legend('Numerical Core = Cai Core', ...
       'Numerical Core = 2*(Cai Core)', ...
       'Numerical Core = 0.5*(Cai Core)', ...
       'location','southeast');
