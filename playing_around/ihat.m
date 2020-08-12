%[h,I]=ihat(H,k,r,n;I)
%---------------------
%
%Inverse Hankel transform of order n.
%
%Input:
% H      Spectrum K(k)
% k      Spatial frequencies [rad/m]   {pi/numel(H)*(0:numel(H)-1)}
% r      Radial positions [m]          {0:numel(H)-1}
% n      Transform order               {0}
%   or
% I      Integration kernel �          {default}
%
%Output:
% h      Signal h(r)
% I      Integration kernel
%
%
% �)  If the integration kernel is missing, it is
%     recomputed from the Bessel functions (slow).
%

%     Marcel Leutenegger � June 2006
%
function [h,I]=ihat(H,k,r,n)
if sum(size(H) > 1) > 1
   error('Spectrum must be a vector.');
end
if nargin < 2 | isempty(k)
   disp('test1')
   k=pi/numel(H)*(0:numel(H)-1).';
else
   disp('test2')
   [k,w]=sort(k(:));
   H=H(w);
end
if nargin < 3 | isempty(r)
   disp('test3')
   r=0:numel(H)-1;
end
if nargin < 4 | isempty(n)
   disp('test4')
   n=0;
end
if numel(n) > 1
   if exist('w','var')
      disp('test5')
      I=n(w,:);
   else
      disp('test6')
      I=n;
   end
else
   I=besselj(n,k*r(:).');
end
h=reshape(H*I,size(r));
