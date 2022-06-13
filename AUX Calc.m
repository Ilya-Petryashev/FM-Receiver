%% Входная цепь
clear
clc 

w0 = 2 * pi * 250.7e6;
Q = 20;
L = 1e-9;
C = 1 / ((w0)^2 * L);
R = Q / sqrt(C / L);

%% Эмиттерный повторитель(вц + урч, упч + ао)
clc
clear

E = 3.5;
Ie = 5e-3;
Ue = E/2; 
Re = Ue/Ie       % сопротивление эмиттера
Ikl = Ie; 
beta = 147.8;
Ib = Ikl/beta      % ток на базе
Idel = Ib*10        % ток делителя напряжения, берём его в 10 раз больше 
Ub_rt = 780.359e-3;         % напряжение из РТ транзистора
Ub = Ue + Ub_rt      % напряжение на базе
Rbaz_e = Ub/Idel      % сопротивление база-эмиттер
Rbaz = (Rbaz_e * (E - Ub)) / Ub     %сопротивление базы

%% УРЧ (П) - новая добротность
clear
clc

E = 3.5;
Polosa = 2 * 40e3 * (1 + 8 + sqrt(8));
Polosa = 1e6; % задаем чуть больше расчетной
f0 = 250.7e6;
Q = f0/Polosa;
Q = 250;    

K = 10;
Ic_Rt = 5e-3;
Ube_Rt = 774.659e-3;
S = 0.112;
Ic_m = Ic_Rt / 3;
Recv = K / S;
Rn = 5e3; % входное сопротивление смесителя 2-5кОм
Rc = (Recv * Rn) / (Recv + Rn)
Uce_m = Rc * Ic_m;
Uce = Uce_m * 3;
Re = (E - Uce - (Rc * Ic_Rt)) / Ic_Rt;
Ue = Re * Ic_Rt;
Ub = Ube_Rt + Ue;
R1 = 10e3;
R2 = (R1 * (E - Ub)) / Ub
Ck = Q / (2 * pi * f0 * Rc)
Lk = 1 / ((2 * pi * f0) ^ 2 * Ck)
S_proverka = K / Rc;
%% Гетеродин
fg = 240e6; Ls = 10e-9; QL = 15; fr = 2e9;
G_sub = 4e4; C_sub = 1.6e-6; p_al = 0.028e-6;m = 4e-7 * pi * 1.000023;
sigma_al = 3.571e7; 
t = 1e-6; tm = 0.9e-6; w = 12e-6; t_ox = 35e-9; delta = sqrt(p_al / (pi * m * fg));
% расчет индуктивности
SL = QL / fg
Rs = 6.06 * (Ls / SL)

delta = sqrt(p_al / (pi * m * fg));
l = Rs * w * sigma_al * delta * (1 - exp(- t / delta))

e_ox = 35e-12;
C_ox = (w * l * e_ox) / t_ox

C_tot =  1 / ((2 * pi * fr) ^ 2 * Ls)

R_si = 2 / (w * l * G_sub)

C_si = (w * l * C_sub) / 2

Rp = 1 / (( 2* pi * fg) ^ 2 * C_ox ^ 2 * R_si) + (R_si * (C_ox + C_si) ^ 2) / (C_ox ^ 2)

Cp = C_ox * (1 + ((2 * pi * fg) ^ 2 * (C_ox + C_si) * C_si * R_si ^ 2)) / (1 + ((2 * pi * fg) ^ 2 * (C_ox + C_si) ^ 2 * R_si ^ 2))

C_s = C_tot - Cp 

CL = C_s + Cp;
R = Rs + Rp;
Y1 = (1i * (2 * pi * fg) * CL) + ((R + 1i * (2 * pi * fg) * Ls)/(Rp * (Rs + 1i * (2 * pi * fg) * Ls)));
B1 = imag(Y1);
B2 = imag(Y1);
G11 = 1e-6;
C = (-B1 / 2) / (2 * pi * fg);
gm = G11 
W = (gm ^ 2 * 0.4e-6) / (2 * C_ox * 4.384e-2 * 1.5e-3)
WM3 = W * 2
%% Смеситель
clc
clear 

Ugs_U0 = 250e-3; % общее для всех транзисторов
I7 = 1e-3;
L7 = 0.4e-6; % общее для всех транзисторов
C_ox = 35.1e-12 / 35e-9; % общее для всех транзисторов
U0 = 0.465; % общее для всех транзисторов
E = 3.5;

K7 =  I7 / Ugs_U0 ^ 2 ;
W7 = (K7 * L7) / (0.5 * C_ox * 4.384e-2);

I3 = I7 / 2; I6 = I3;

K6 =  I6 / Ugs_U0 ^ 2 ; K3 = K6;
W6 = (K6 * L7) / (0.5 * C_ox * 4.384e-2); W3 = W6;

Ug_min = sqrt(2) * (Ugs_U0);

I1 = I7 / 4; I2 = I1; I4 = I1; I5 = I1;

K1 =  I1 / Ugs_U0 ^ 2 ; K2 = K1 ; K4 = K1; K5 = K1;
W1 = (K1 * L7) / (0.5 * C_ox * 4.384e-2); W2 = W1 ; W4 = W1; W5 = W1;

Ugate1245_min = Ugs_U0 + Ugs_U0 + (Ugs_U0 + U0);
Rn_max = (E - (3 + sqrt(2))*(Ugs_U0)) / (I1 + I4);
Rn = 4000;
Ugate1245_max = E - (I1 + I4) * Rn + U0;
Ugetm = Ugate1245_max - Ugate1245_min;

Ugate36_min = Ugs_U0 + (Ugs_U0 + U0);
Ugate36_max = E - (I1 + I4) * Rn - (Ugs_U0 + U0) + 2*U0;

Uc = Ugate36_max - Ugate36_min;

I8 = I7;
Rsource = (E - (Ugs_U0 + U0)) / I8;

gm3 = sqrt(2 * I3 * C_ox * 4.384e-2 * ( W3 / L7 ));
Kcm = (2 / pi) * gm3 * Rn;
%% Эмиттерный повторитель(пч + упч)
clc
clear

E = 3.5;
Ie = 0.875e-3;
Ue = E/2; 
Re = Ue/Ie       % сопротивление эмиттера
Ikl = Ie; 
beta = 147.8;
Ib = Ikl/beta      % ток на базе
Idel = Ib * 10        % ток делителя напряжения, берём его в 10 раз больше 
Ub_rt = 731.185e-3;         % напряжение из РТ транзистора
Ub = Ue + Ub_rt      % напряжение на базе
Rbaz_e = Ub/Idel      % сопротивление база-эмиттер
Rbaz = (Rbaz_e * (E - Ub)) / Ub     %сопротивление базы
%% УПЧ (П)
clear
clc

E = 3.5;
Q = 10;
f0 = 10.7e6;
K = 40;
Ic_Rt = 20e-3;
Ube_Rt = 831.621e-3;
S = 0.429;
Ic_m = Ic_Rt / 3;
Recv = K / S;
Rn = 200;
Rc = (Recv * Rn) / (Recv + Rn)
Uce_m = Rc * Ic_m;
Uce = Uce_m * 3;
Re = (E - Uce - (Rc * Ic_Rt)) / Ic_Rt;
Ue = Re * Ic_Rt;
Ub = Ube_Rt + Ue;
R1 = 10e3;
R2 = (R1 * (E - Ub)) / Ub
Ck = Q / (2 * pi * f0 * Rc)
Lk = 1 / ((2 * pi * f0) ^ 2 * Ck)
S_proverka = K / Rc;
%% Эмиттерный повторитель(ао + дет)
clc
clear

E = 3.5;
Ie = 20e-3;
Ue = E/2; 
Re = Ue/Ie       % сопротивление эмиттера
Ikl = Ie; 
beta = 147.8;
Ib = Ikl/beta      % ток на базе
Idel = Ib*10        % ток делителя напряжения, берём его в 10 раз больше 
Ub_rt = 831.621e-3;         % напряжение из РТ транзистора
Ub = Ue + Ub_rt      % напряжение на базе
Rbaz_e = Ub/Idel      % сопротивление база-эмиттер
Rbaz = (Rbaz_e * (E - Ub)) / Ub     %сопротивление базы
%% Детектор (Я)
clear
clc
f0 = 10.7e6;
F = 40e3;
L3 = 5e-9;
C3 = 1 / ((2 * pi * f0) ^ 2 * L3);
C1 = C3 * 100;
C2 = C1;
L21 = L3 / 200; L22 = L21;
L1 = L3 / 100;
RC_min = 1 / f0; % 93n
RC_max = 1 / F; % 25u
RC = 10e-6;
Rn = 500;
Cn = RC / Rn;
Cl3 = 10e-6;



