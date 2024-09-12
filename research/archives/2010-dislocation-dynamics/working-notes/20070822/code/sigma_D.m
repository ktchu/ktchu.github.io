%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
%
% Function sigma_D accepts as inputs:
% (a) the burgers vector for parallel edge dislocations and its position
% in the periodic simulation box
% (b) the spacing between the dislocation lines with its periodic image D
% (c) points [x] and [y].
% The outputs are:
% (a) sigma_xx
% (b) sigma_xy
% (c) sigma_yy
% for the points (x, y) in the domain.
%
% Note: mu/(2*(1-nu)) has been factored out of the expressions for the
% stress.
%
% Adele T. Lim
% 08/2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma_xx, sigma_xy, sigma_yy] = sigma_D (b, r, D, x, y)

x = ( x - r(1) )/D;
y = ( y - r(2) )/D;

[X, Y] = meshgrid(x,y);

        sigma_0_D = 1/D./( cosh(2*pi*X) - cos(2*pi*Y) ).^2; %prefactor for sigma of D

        sigma_xx = -sigma_0_D.*(...
                            b(1)*sin(2*pi*Y).*( cosh(2*pi*X) - cos(2*pi*Y) + 2*pi*X.*sinh(2*pi*X) )...
                            -b(2)*2*pi*X.*( cosh(2*pi*X).*cos(2*pi*Y) - 1 )...
                         );

        sigma_xy =  sigma_0_D.*(...
                            b(1)*2*pi*X.*( cosh(2*pi*X).*cos(2*pi*Y) - 1 )...
                            -b(2)*sin(2*pi*Y).*( cosh(2*pi*X) - cos(2*pi*Y) - 2*pi*X.*sinh(2*pi*X) )...
                         );

        sigma_yy = -sigma_0_D.*(...
                            b(1)*sin(2*pi*Y).*( cosh(2*pi*X) - cos(2*pi*Y) - 2*pi*X.*sinh(2*pi*X) )...
                            -b(2)*( 2*sinh(2*pi*X).*( cosh(2*pi*X) - cos(2*pi*Y) ) - 2*pi*X.*( cosh(2*pi*X).*cos(2*pi*Y) - 1 ) )...
                         );

