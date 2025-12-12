%% RL Circuit 
clearvars;
%Vars
h = 1/192000; %Sample Time
R = 100; %  Resistance (Ω)
L = 0.1; % Inductance (H)
k_simulation = 0.01; %Simulation Time Span
i_0 = 0; %Initial Current

%System Matrices
A = [1-h*R/L];
B = [h/L];
C = [-R];
D = [1];

%Initialize Vars
k = 0:h:k_simulation; %Time Vector
V_in = ones(size(k));

%Create System
rl_circuit = ss(A, B, C, D, h);

%Run System
[sysout] = lsim(rl_circuit, V_in, k, i_0);

%Save Results to Their Respective Vars
V_R = sysout;
