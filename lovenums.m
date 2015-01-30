function [lkli,lkl]=lovenums(mod,els)
% [lkli,lkl]=LOVENUMS(mod,els)
%
% Returns a set of Love numbers. Intermediate numbers can be linearly
% interpolated with  errors of less than 0.05% for all l<200, as opposed
% to direct calculation. See Wahr et al. (1998) and Han & Wahr (1995).
%
% INPUT:
%
% mod   'Wahr' Wahr et al. (1998) after Han & Wahr (1995) after PREM.
% els   The spherical harmonic degrees where you want them
%
% OUTPUT:
%
% lkli   The spherical harmonic degrees l under consideration and k_l, 
%        the elastic Love numbers in PREM, linearly interpolated from:
% lkl    The original table from the model specified as the mod input.
%
% Last modified by charig-at-princeton.edu, 01/25/2010
% Last modified by fjsimons-at-alum.mit.edu, 02/21/2012

% Only one option, really
defval('mod','Wahr')

switch mod
  case 'Wahr'
    lkl=[0    0.000;
	 1    0.027;
	 2   -0.303;
	 3   -0.194;
	 4   -0.132;
	 5   -0.104;
	 6   -0.089;
	 7   -0.081;
	 8   -0.076;
	 9   -0.072;
	 10  -0.069;
	 12  -0.064;
	 15  -0.058;
	 20  -0.051;
	 30  -0.040;
	 40  -0.033;
	 50  -0.027;
	 70  -0.020;
	 100 -0.014;
	 150 -0.010;
	 200 -0.007];
 otherwise
  error('Specify valid Love number table')
end

% Now interpolate to where you want them
lkli=[els(:) interp1(lkl(:,1),lkl(:,2),els(:))];



