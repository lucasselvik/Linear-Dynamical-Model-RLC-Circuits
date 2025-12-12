%% RC Circuit 
clearvars;

%Vars
h = 0.0001; %Sample Time
R = 1000; %  Resistance (Ω)
C = 1e-6; % Capacatance (F)
k_simulation = 0.005; %Simulation Time Span
V_C0 = 0; %Initial Capacitor Voltage

%System Matrices
A = 1 - h/(R*C);
B = h/(R*C);
C = [-1;
      1];
D = [1;...
     0];

%Initialize Vars
k = 0:h:k_simulation; %Time Vector
V_in = ones(size(k));

%Create System
rc_circuit = ss(A, B, C, D, h);

%Run System
[sysout] = lsim(rc_circuit, V_in, k, V_C0);

%Save Results to Their Respective Vars
V_C = sysout(:, 2);
V_R = sysout(:, 1);