%% Case study 3: Circuits as Resonators, Sensors, and Filters
% *ESE 105* 
%
% *Name: Lucas Selvik and Rex Paster*
%
% function myResonatorCircuit(Vin,h) receives a time-series voltage sequence
% sampled with interval h, and returns the output voltage sequence produced
% by a circuit
%
% inputs:
% Vin - time-series vector representing the voltage input to a circuit
% h - scalar representing the sampling interval of the time series in
% seconds
%
% outputs:
% Vout - time-series vector representing the output voltage of a circuit

function Vout = myResonatorCircuit(Vin,h)
    %Set Vars 
    DR = 0.00715; % Dampening Ratio
    freq = 440;  % Note To Play (hz)
    L = 0.1; % Inductance (H)
    
    %Calculate C & R
    C = 1/(((2*pi*freq)^2)*L); %Uses the Resonant LC Freq Equation
    R = 2*DR*sqrt(L/C);    % Uses the Dampening Ratio Equation

    %Call RLC Sim
    Vout = simulateRLC(Vin, h, C, L, R);
end

%% simulateRLC Local Function Definition
function V_Out = simulateRLC(V_in, h, C, L, R, V_C0, i_0)
% simulateRLC Simulates an RLC circuit using discrete-time state-space modeling
%
% Inputs:
%   V_in        - Input voltage Discrete Time Matrix (V)
%   h           - Sample rate (s)
%   C           - Capacitance (F)
%   L           - Inductance (H)
%   R           - Resistance (Ohms)
%   V_C0        - Initial capacitor voltage (V) (Optional: Default = 0)
%   i_0         - Initial inductor current (A)  (Optional: Default = 0)
%
% Outputs:
%   V_Out       - Voltage Output Discrete Time Matrix (V)


%% --- Define Unspecified Arguments & Simulation Variables ---
timeVec = (0:length(V_in)-1)*h;

% Check for optional arguments
if (nargin < 6)          % V_C0 not provided
   V_C0 = 0;          % default capacitor voltage
end
if (nargin < 7)          % i_0 not provided
   i_0 = 0;           % default inductor current
end

%% --- Define Discrete-Time State-Space Matrices ---
% State variables: x = [V_C; i_L]
% V_C: voltage across capacitor
% i_L: current through inductor

A = [ 1       , h/C;       % State transition matrix
     -h/L     , 1 - h*R/L ];

B = [0;           % Input matrix (v_in affects inductor current)
     h/L];

C = [0, R];       % Output is voltage across resistor

D = 0;            % V_in has no istantaneous affect on V_out

%% --- Create State-Space System ---
rlc_circuit = ss(A, B, C, D, h);

%% --- Run Simulation ---
% lsim computes the time response of the discrete-time system
% Initial state: [V_C0; i_0]
sysout = lsim(rlc_circuit, V_in, timeVec, [V_C0; i_0]);

%% --- Extract Outputs ---
V_R = sysout;    % Voltage across resistor
V_Out = V_R; 

end

% Note: Some comments in this function were generated with the assistance of ChatGPT and
% have been reviewed for accuracy and clarity. All actionable code is
% authored by Lucas Selvik and/or Rex Paster.

