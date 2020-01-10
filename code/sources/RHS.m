%{

    Copyright (C) 2020 Denys DUTYKH

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

%}

%%% -------------------------------------------------- %%%
%%% Right-hand side of the Dysthe-Lo-Mei equation      %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, USMB           %%%
%%% E-mail: Denys.Dutykh@univ-smb.fr                   %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Distributed under GNU General Public License       %%%
%%% -------------------------------------------------- %%%

%%% Right-hand side of the ODE system
function rhs = RHS (w, dt)

    global e g k op

    v_hat = turn(w, -dt);
    v     = ifft(v_hat);
    vx    = ifft(1i*k.*v_hat);
    va2   = conj(v).*v;
    Hv    = ifft(op.*fft(va2));
    rh    = fft(-1i*va2.*v - 8.0*e*g*va2.*vx - 4i*e*g*v.*Hv);
    rhs   = turn(rh, dt);

end % RHS ()