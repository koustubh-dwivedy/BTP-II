global state;
state = [0, 0, 1]; %because T_boil =276.6K and initial conditions are at 300K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t y]=ode45(@F,[0 7200], [0 500 300 300 300 300 300 0]);
% y(1)-> dz/dt
% y(2)-> z
% y(3)-> T_pg
% y(4)-> T_pf
% y(5)-> T_sg
% y(6)-> T_sl
% y(7)-> T_sf
% y(8)-> m_sl