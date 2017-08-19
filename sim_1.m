z = 0;

m_sys = 0;
c_m = 0;
rho_inf = 0;
V = 0;
g = 0;
C_D = 0;
W_z = 0;
A = 0;
m_tot = 0;
L_dot_l = 0;
L_dot_g = 0;

V_g = 0;
V_l = 0;

m_g = 0;
R = 0;
T_g = 0;
M_g = 0;
P_inf = 0;
m_l = 0;
rho_l = 0;
U_z = 0;

m_pg = 0;
c_pg = 0;
T_pg = 0;
q_dot_pg = 0;
V_pg = 0;

m_pf = 0;
c_pf = 0;
T_pf = 0;
q_dot_pf = 0;

c_g = 0;
T_g = 0;
q_dot_g = 0;

c_l = 0;
T_l = 0;
q_dot_l = 0;

m_f = 0;
c_f = 0;
T_f = 0;
q_dot_f = 0;

G = 0;
alpha_pgeff = 0;
r_v = 0;
epsilon_pint = 0;
sigma = 0;
T_gp = 0;


g = 9.8;
m_sys = 3.0054;
% y(1) -> dz/dt
% y(2) -> z
% y(3) -> T_pg
% y(4) -> T_pf
% y(5)-> T_g
% y(6)-> T_l
% y(7)-> T_f
F=@(t,y) [
;
;
;
;
];


F=@(t,y) [t.*y(1)./(y(1).^2+y(2).^2+1);
(y(2)-y(1)).^2./(y(2).^2+y(3).^2+1);
t.^2.*y(3).^2./(y(1).^2+y(3).^2+1)];

[t y]=ode45(F,[0 7200], [1 0 -1]);