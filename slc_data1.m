    %%
    %% Data set 1: four bar linkage
    %%
r1 = 8;
r2 = 3.5;
r3 = 10.0;
r4=9.0;
thetas2=15*pi/180;
theta2dot=36.7;
theta2ddot=0;
guess1=45*pi/180;
guess2=60*pi/180;
init_values=[guess1;guess2];
rho=7700; %density of stainless steel
area=0.001; %cross sectional area
b2 = r2/2; m2 = rho*area*r2; phi2 = 0; Ig2 = (1/12)*m2*r2.^2;
b3 = r3/2; m3 = rho*area*r3; phi3 = 0; Ig3 = (1/12)*m3*r3.^2;
b4 = r4/2; m4 = rho*area*r4; phi4 = 0; Ig4 = (1/12)*m4*r4.^2;
mu = 30;
P = 3;

    %%
    %% Data set 2: four bar linkage
    %%
% r1 = 5;
% r2 = 1;
% r3 = 4;
% r4=3;
% thetas2=30*pi/180;
% theta2dot=30;
% theta2ddot=0;
% guess1=45*pi/180;
% guess2=60*pi/180;
% init_values=[guess1;guess2];
% b2 = 1.25; m2 = 2.46; phi2 = 0; Ig2 = 0;
% b3 = 6.67; m3 = 9.84; phi3 = 0; Ig3 = 82.0;
% b4 = 8.89; m4 = 9.95; phi4 = 0; Ig4 = 50.0;
% m4 = 10.0;
% mu = 30;
% P = 3;

    %%
    %% Data set 3: four bar linkage
    %%
% r1 = 8;
% r2 = 2;
% r3 = 5;
% r4 = 4;
% thetas2=60*pi/180; %deg -> rad
% theta2dot=45; %rad/s
% theta2ddot=0;
% guess1=45*pi/180;
% guess2=60*pi/180;
% init_values=[guess1;guess2];
% b2 = 1.25; m2 = 2.46; phi2 = 0; Ig2 = 0;
% b3 = 6.67; m3 = 9.84; phi3 = 0; Ig3 = 82.0;
% b4 = 8.89; m4 = 9.95; phi4 = 0; Ig4 = 50.0;
% m4 = 10.0;
% mu = 30;
% P = 3;

    %%
    %% Data set 4: four bar linkage
    %%
% r1 = 10;
% r2 = 32;
% r3 = 20;
% r4 = 20;
% thetas2=120*pi/180; %deg -> rad
% theta2dot=-50; %rad/s
% theta2ddot=40; %rad/s^2
% guess1=45*pi/180;
% guess2=60*pi/180;
% init_values=[guess1;guess2];
% b2 = 1.25; m2 = 2.46; phi2 = 0; Ig2 = 0;
% b3 = 6.67; m3 = 9.84; phi3 = 0; Ig3 = 82.0;
% b4 = 8.89; m4 = 9.95; phi4 = 0; Ig4 = 50.0;
% m4 = 10.0;
% mu = 30;
% P = 3;

    %%
    %% Data set 5: four bar linkage
    %%
% r1 = 5;
% r2 = 1;
% r3 = 6;
% r4 = 4;
% thetas2=210*pi/180; %deg -> rad
% theta2dot=-100; %rad/s
% theta2ddot=-40; %rad/s^2
% guess1=45*pi/180;
% guess2=60*pi/180;
% init_values=[guess1;guess2];
% b2 = 1.25; m2 = 2.46; phi2 = 0; Ig2 = 0;
% b3 = 6.67; m3 = 9.84; phi3 = 0; Ig3 = 82.0;
% b4 = 8.89; m4 = 9.95; phi4 = 0; Ig4 = 50.0;
% m4 = 10.0;
% mu = 30;
% P = 3;



