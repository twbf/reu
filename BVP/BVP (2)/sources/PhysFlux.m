%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Distributed under GNU General Public License       %%%
%%% -------------------------------------------------- %%%

function flux = PhysFlux (w)
    
    global g2
    
    flux = 0.0*w;
    
    flux(:,1) = w(:,2);
    flux(:,2) = w(:,2).^2./(w(:,1)+eps) + g2*w(:,1).^2;
    
end % PhysFlux ()
