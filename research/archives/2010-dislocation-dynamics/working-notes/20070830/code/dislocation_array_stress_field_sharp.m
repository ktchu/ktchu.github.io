%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% dislocation_array_stress_field_sharp
%
% Description: MATLAB function for calculating the stress field of edge
% dislocation array when each dislocation core is infinitely thin.
%
% Input arguments:
%   b              = a row vector for the burgers vector b of the array
%                    of parallel edge dislocations
%   r              = a position vector for the dislocation in the simulation
%                    box
%   D              = the spacing between dislocation lines in the infinite
%                    dislocation array which along the y-axis
%   X, Y           = computational grid (as generated by meshgrid()) on
%                    which to calculate the stress
%   G              = shear modulus
%   poisson_ratio  = poisson ratio
%   a              = core radius
%
% Returns values:
%   sigma_xx = xx-component of stress
%   sigma_yy = yy-component of stress
%   sigma_xy = xy-component of stress
%
% Function dislocation_array_stress_field_Cai assumes that the Burgers
% vector lies in the x-y plane with the tangent vector in the z-direction.
% The dislocation array is periodic and inifinite along the y-axis.
% The stress field is the infinite sum of the stress field due to a single
% singular edge dislocation with the specified core radius (Hirth & Lothe 
% 2006).
%
% Notes:
% (a) Reference:  J. Hirth & J. Lothe, Theory of Dislocations, 1982.
%
% File: dislocation_array_stress_field_sharp.m
% Written by: Adele T. Lim, MAE, Princeton University
% Date modified: Aug 2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHANGE LOG
% - 2007/08/30: Kevin T. Chu
%   - Modified to take physical parameters are arguments.
%   - Cleaned up code and documentation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma_xx, sigma_yy, sigma_xy] = ...
  dislocation_array_stress_field_sharp (b, r, D, X, Y, G, poisson_ratio)

X = ( X - r(1) )/D;
Y = ( Y - r(2) )/D;


sigma_0 = 0.5*G/(1-poisson_ratio)/D./( cosh(2*pi*X) - cos(2*pi*Y) ).^2; 

sigma_xx = -sigma_0.*(...
   b(1)*sin(2*pi*Y).*( cosh(2*pi*X) - cos(2*pi*Y) + 2*pi*X.*sinh(2*pi*X) )...
  -b(2)*2*pi*X.*( cosh(2*pi*X).*cos(2*pi*Y) - 1 ) );

sigma_xy =  sigma_0.*(...
   b(1)*2*pi*X.*( cosh(2*pi*X).*cos(2*pi*Y) - 1 )...
  -b(2)*sin(2*pi*Y).*( cosh(2*pi*X) - cos(2*pi*Y) - 2*pi*X.*sinh(2*pi*X) ));

sigma_yy = -sigma_0.*(...
   b(1)*sin(2*pi*Y).*( cosh(2*pi*X) - cos(2*pi*Y) - 2*pi*X.*sinh(2*pi*X) )...
  -b(2)*( 2*sinh(2*pi*X).*( cosh(2*pi*X) - cos(2*pi*Y) ) ...
        - 2*pi*X.*( cosh(2*pi*X).*cos(2*pi*Y) - 1 ) ) );

