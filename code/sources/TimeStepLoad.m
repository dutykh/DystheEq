%{
    Coefficients of a very high order Runge-Kutta scheme

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
%%% Coefficients of the Verner 9(8) Runge-Kutta method %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, USMB           %%%
%%% E-mail: Denys.Dutykh@univ-smb.fr                   %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Distributed under GNU General Public License       %%%
%%% -------------------------------------------------- %%%

% Coefficients of the Verner 16-stages, 9(8) method
s6 = sqrt(6);

c1 = 103/1680;
c8 = -27/140;
c9 = 76/105;
c10 = -201/280;
c11 = 1024/1365;
c12 = 3/7280;
c13 = 12/35;
c14 = 9/280;

a2 = 1/12;
a3 = 1/9;
a4 = 1/6;
a5 = 2*(1+s6)/15;
a6 = (6+s6)/15;
a7 = (6-s6)/15;
a8 = 2/3;
a9 = 1/2;
a10 = 1/3;
a11 = 1/4;
a12 = 4/3;
a13 = 5/6;
a15 = 1/6;

b21 = 1/12;

b31 = 1/27;
b32 = 2/27;

b41 = 1/24;
b43 = 1/8;

b51 = (4+94*s6)/375;
b53 = -(94+84*s6)/125;
b54 = (328+208*s6)/375;

b61 = (9-s6)/150;
b64 = (312+32*s6)/1425;
b65 = (69+29*s6)/570;

b71 = (927-347*s6)/1250;
b74 = (-16248+7328*s6)/9375;
b75 = (-489+179*s6)/3750;
b76 = (14268-5798*s6)/9375;

b81 = 2/27;
b86 = (16-s6)/54;
b87 = (16+s6)/54;

b91 = 19/256;
b96 = (118-23*s6)/512;
b97 = (118+23*s6)/512;
b98 = -9/256;

b10_1 = 11/144;
b10_6 = (266-s6)/864;
b10_7 = (266+s6)/864;
b10_8 = -1/16;
b10_9 = -8/27;

b11_1 = (5034-271*s6)/61440;
b11_7 = (7859-1626*s6)/10240;
b11_8 = (-2232+813*s6)/20480;
b11_9 = (-594+271*s6)/960;
b11_10 = (657-813*s6)/5120;

b12_1 = (5996-3794*s6)/405;
b12_6 = -(4342+338*s6)/9;
b12_7 = (154922-40458*s6)/135;
b12_8 = (-4176+3794*s6)/45;
b12_9 = (-340864+242816*s6)/405;
b12_10 = (26304-15176*s6)/45;
b12_11 = -26624/81;

b13_1 = (3793+2168*s6)/103680;
b13_6 = (4042+2263*s6)/13824;
b13_7 = (-231278+40717*s6)/69120;
b13_8 = (7947-2168*s6)/11520;
b13_9 = (1048-542*s6)/405;
b13_10 = (-1383+542*s6)/720;
b13_11 = 2624/1053;
b13_12 = 3/1664;

b14_1 = -137/1296;
b14_6 = (5642-337*s6)/864;
b14_7 = (5642+337*s6)/864;
b14_8 = -299/48;
b14_9 = 184/81;
b14_10 = -44/9;
b14_11 = -5120/1053;
b14_12 = -11/468;
b14_13 = 16/9;

b15_1 = (33617-2168*s6)/518400;
b15_6 = (-3846+31*s6)/13824;
b15_7 = (155338-52807*s6)/345600;
b15_8 = (-12537+2168*s6)/57600;
b15_9 = (92+542*s6)/2025;
b15_10 = -(1797+542*s6)/3600;
b15_11 = 320/567;
b15_12 = -1/1920;
b15_13 = 4/105;

b16_1 = -(36487+30352*s6)/279600;
b16_6 = -(29666+4499*s6)/7456;
b16_7 = (2779182-615973*s6)/186400;
b16_8 = (-94329+91056*s6)/93200;
b16_9 = (-232192+121408*s6)/17475;
b16_10 = (101226-22764*s6)/5825;
b16_11 = -169984/9087;
b16_12 = -87/30290;
b16_13 = 492/1165;
b16_15 = 1260/233;

e1 = -7/400;
e8 = 63/200;
e9 = -14/25;
e10 = 21/20;
e11 = -1024/975;
e12 = -21/36400;
e13 = -3/25;
e14 = -9/280;
e15 = 9/25;
e16 = 233/4200;

% a constant for adaptive time stepping
C36 = 1/36;