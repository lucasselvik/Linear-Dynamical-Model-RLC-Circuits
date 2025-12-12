% The Following Code Was Used For The Development of Figures used In our Report

%% RLC Based Figures
h = 1/192000;            % Sampling interval (sample time)
k_simulation = 0.015;    % Total simulation time (seconds)
k = 0:h:k_simulation;    % Time vector
V_in = ones(size(k));    % Step input voltage (constant 1V)
V_C0 = 0;                % Initial capacitor voltage (Volts)
i_0 = 0;                 % Initial current through inductor (Amps)

% Base component values for R, L, and C
RMaster = 5.2;           % Resistance in Ohms
CMaster = 1e-6;          % Capacitance in Farads
LMaster = 0.1;           % Inductance in Henries

% Plot RLC circuit step response varying R, C, and L parameters
figure(Name="Step Input RLC Variations");
sgtitle('RLC Circuit Response for Various Component Values with Step Input');

% Vary resistance R: low, medium, and high values with fixed L and C
subplot(3, 3, 1);
R = RMaster/5; C = CMaster; L = LMaster;
plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5)
title("R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=16)
ylabel("Output Voltage (V)", FontSize=16)
xlabel("Time (s)", FontSize=16);

subplot(3, 3, 2);
R = RMaster; C = CMaster; L = LMaster;
plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5)
title("R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=16)
ylabel("Output Voltage (V)", FontSize=16)
xlabel("Time (s)", FontSize=16);

subplot(3, 3, 3);
R = RMaster*5; C = CMaster; L = LMaster;
plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5)
title("R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=16)
ylabel("Output Voltage (V)", FontSize=16)
xlabel("Time (s)", FontSize=16);

% Vary capacitance C: low, medium, and high values with fixed R and L 
subplot(3, 3, 4);
R = RMaster; C = CMaster/5; L = LMaster;
plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5)
title("R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=16)
ylabel("Output Voltage (V)", FontSize=16)
xlabel("Time (s)", FontSize=16);

subplot(3, 3, 5);
R = RMaster; C = CMaster; L = LMaster;
plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5)
title("R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=16)
ylabel("Output Voltage (V)", FontSize=16)
xlabel("Time (s)", FontSize=16);

subplot(3, 3, 6);
R = RMaster; C = CMaster*5; L = LMaster;
plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5)
title("R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=16)
ylabel("Output Voltage (V)", FontSize=16)
xlabel("Time (s)", FontSize=16);

% Vary inductance L: low, medium, and high values with fixed R and C
subplot(3, 3, 7);
R = RMaster; C = CMaster; L = LMaster/5;
plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5)
title("R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=16)
ylabel("Output Voltage (V)", FontSize=16)
xlabel("Time (s)", FontSize=16);

subplot(3, 3, 8);
R = RMaster; C = CMaster; L = LMaster;
plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5)
title("R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=16)
ylabel("Output Voltage (V)", FontSize=16)
xlabel("Time (s)", FontSize=16);

subplot(3, 3, 9);
R = RMaster; C = CMaster; L = LMaster*5;
plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5)
title("R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=16)
ylabel("Output Voltage (V)", FontSize=16)
xlabel("Time (s)", FontSize=16);

% Plot step responses of three tuned RLC configurations on one figure for comparison
figure(Name="Three Tuned Responses");
sgtitle('Tuning RLC Components to Create Diverse Responses with Step Input');
hold on;

% First tuned configuration with increased R and reduced L
R1 = RMaster*6; C1 = CMaster; L1 = LMaster*0.2;
V_out1 = simulateRLC(V_in, h, C1, L1, R1, V_C0, i_0);
plot(k, V_out1, 'LineWidth', 1.5, 'DisplayName', "R: " + R1 + " Ω, L: " + L1 + " H, C: " + C1 + " F");

% Second tuned configuration: base values
R2 = RMaster; C2 = CMaster; L2 = LMaster;
V_out2 = simulateRLC(V_in, h, C2, L2, R2, V_C0, i_0);
plot(k, V_out2, 'LineWidth', 1.5, 'DisplayName', "R: " + R2 + " Ω, L: " + L2 + " H, C: " + C2 + " F");

% Third tuned configuration: reduced capacitance
R3 = RMaster; C3 = CMaster/10; L3 = LMaster;
V_out3 = simulateRLC(V_in, h, C3, L3, R3, V_C0, i_0);
plot(k, V_out3, 'LineWidth', 1.5, 'DisplayName', "R: " + R3 + " Ω, L: " + L3 + " H, C: " + C3 + " F");

% Plot the input voltage in red for reference
plot(k, V_in, 'LineWidth', 1.5, 'DisplayName', 'Input Voltage', 'Color', 'red');

hold off;
ylabel("Output Voltage (V)", FontSize=12);
xlabel("Time (s)", FontSize=12);
legend('show');

% Analysis of RLC circuit responses for sinusoidal inputs at logarithmically spaced frequencies
frequencies = logspace(log10(100), log10(10000), 9); % 100 Hz to 10 kHz
R = 100; L = 100e-3; C = 0.1e-6; % Base component values
Fs = 192000; % Sampling frequency (Hz)
k_simulation = 0.005; % Simulation time span
k = 0:h:k_simulation; % Time vector
figure(Name="Base RLC Responses");
sgtitle('RLC Circuit Response at Various Frequencies with Oscillating Input');

for i = 1:length(frequencies)
    f = frequencies(i);
    V_in = sin(2 * pi * f * k); % Sinusoidal input at frequency f
    % Adjust simulation time for high frequencies to capture waveform accurately
    if i == 7; k_simulation = 0.0015; end
    if i == 8; k_simulation = 0.001; end
    if i == 9; k_simulation = 0.0006; end
    k = 0:h:k_simulation;
    V_in = sin(2 * pi * f * k);
    subplot(3, 3, i);
    hold on
    plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5)
    plot(k, V_in, LineWidth=1.5)
    hold off
    title("f: " + f + " Hz, R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=14)
    ylabel("Voltage (V)", FontSize=12)
    xlabel("Time (s)", FontSize=12)
    legend("Output","Input", FontSize=10)
end

% Calculate gain for base parameters for plotting gain vs frequency
gain_vector = calculateGain(V_in, h, C, L, R, V_C0, i_0, k);

% Calculate resonant frequency for RLC circuit (Hz)
f_resonant = 1 / (2 * pi * sqrt(L * C));

% Plot response at resonant frequency with sinusoidal input
figure(Name="Resonant Frequency Response");
f = f_resonant; % Resonant frequency
k_simulation = 0.01; % Simulation time for resonant frequency
k = 0:h:k_simulation;
V_in = sin(2 * pi * f * k); % Sinusoidal input at resonant frequency
hold on
plot(k, V_in, LineWidth=1.5)
plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5)
hold off
title("Voltages at Resonant Frequency");
subtitle("f: " + f + " Hz, R: " + R + " Ω, L: " + L + " H, C: " + C + " F");
ylabel("Voltage (V)");
xlabel("Time (s)");
legend("Input","Output")

% Increase Resistance, vary frequencies to observe effect on circuit response
frequencies = logspace(log10(100), log10(10000), 9);
R = 1000; % Increased resistance
figure(Name="Increased R RLC Responses");
sgtitle('Circuit Responses to Oscillating Input with Increased Resistance');

for i = 1:length(frequencies)
    f = frequencies(i);
    V_in = sin(2 * pi * f * k);
    % Adjust simulation time span for high frequencies
    if i == 5; k_simulation = 0.006; end
    if i == 6; k_simulation = 0.004; end
    if i == 7; k_simulation = 0.002; end
    if i == 8; k_simulation = 0.001; end
    if i == 9; k_simulation = 0.0006; end
    k = 0:h:k_simulation;
    V_in = sin(2 * pi * f * k);
    subplot(3, 3, i);
    hold on
    plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5);
    plot(k, V_in, LineWidth=1.5);
    hold off
    title("f: " + f + " Hz, R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=12);
    ylabel("Voltage (V)", FontSize=12);
    xlabel("Time (s)", FontSize=12);
    legend("Output","Input", FontSize=10)
end

% Similar frequency response plotting with increased inductance L
frequencies = logspace(log10(100), log10(10000), 9);
R = 100;
L = 100e-2; % Increased inductance
figure(Name="Increased L RLC Responses");
sgtitle('Circuit Responses to Oscillating Input with Increased Inductance');

for i = 1:length(frequencies)
    f = frequencies(i);
    V_in = sin(2 * pi * f * k);
    % Adjust simulation time span for high frequencies
    if i == 5; k_simulation = 0.006; end
    if i == 6; k_simulation = 0.004; end
    if i == 7; k_simulation = 0.002; end
    if i == 8; k_simulation = 0.001; end
    if i == 9; k_simulation = 0.0006; end
    k = 0:h:k_simulation;
    V_in = sin(2 * pi * f * k);
    subplot(3, 3, i);
    hold on
    plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5);
    plot(k, V_in, LineWidth=1.5);
    hold off
    title("f: " + f + " Hz, R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=12);
    ylabel("Voltage (V)", FontSize=12);
    xlabel("Time (s)", FontSize=12);
    legend("Output","Input", FontSize=10)
end

% Similar frequency response plotting with increased capacitance C
frequencies = logspace(log10(100), log10(10000), 9);
R = 100; L = 100e-3; C = 0.1e-5; % Increased capacitance
figure(Name="Increased C RLC Responses");
sgtitle('Circuit Responses to Oscillating Input with Increased Capacitance');

for i = 1:length(frequencies)
    f = frequencies(i);
    V_in = sin(2 * pi * f * k);
    % Adjust simulation time span for high frequencies
    if i == 5; k_simulation = 0.006; end
    if i == 6; k_simulation = 0.004; end
    if i == 7; k_simulation = 0.002; end
    if i == 8; k_simulation = 0.001; end
    if i == 9; k_simulation = 0.0006; end
    k = 0:h:k_simulation;
    V_in = sin(2 * pi * f * k);
    subplot(3, 3, i);
    hold on
    plot(k, simulateRLC(V_in, h, C, L, R, V_C0, i_0), LineWidth=1.5);
    plot(k, V_in, LineWidth=1.5);
    hold off
    title("f: " + f + " Hz, R: " + R + " Ω, L: " + L + " H, C: " + C + " F", FontSize=12);
    ylabel("Voltage (V)", FontSize=12);
    xlabel("Time (s)", FontSize=12);
    legend("Output","Input", FontSize=10)
end

% Plot output gain for different RLC circuit tunings - summary figure
figure(Name="Output Gain");
sgtitle("Output Gain for Various RLC Circuit Tunings", FontSize=20);

subplot(2,2,1);
hold on
title(['5.1: Gain for R = ' num2str(R) ' Ω, L = ' num2str(L) ' H and C = ' num2str(C) ' F'], FontSize=16);
plot(frequencies, gain_vector, 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Gain (%)', 'FontSize', 20);
grid on;
hold off

% Calculate gains for different resistor values (keeping L and C fixed)
R = 1000; L = 100e-3; C = 0.1e-6;
Fs = 192000;
h=1/Fs;
k_simulation = 0.01;
k = 0:h:k_simulation;
frequencies = logspace(log10(100), log10(10000), 9);
GR1 = calculateGain(V_in, h, C, L, 1, V_C0, i_0, k); % Gain for R=1Ω
GR2 = calculateGain(V_in, h, C, L, 10, V_C0, i_0, k);
GR3 = calculateGain(V_in, h, C, L, 100, V_C0, i_0, k);
GR4 = calculateGain(V_in, h, C, L, 1000, V_C0, i_0, k);
GR5 = calculateGain(V_in, h, C, L, 5000, V_C0, i_0, k);

subplot(2,2,2);
hold on
title(['5.2: Gain for L = ' num2str(L) ' H and C = ' num2str(C) ' F'], FontSize=16);
plot(frequencies, GR1, 'LineWidth', 1.5);
plot(frequencies, GR2, 'LineWidth', 1.5);
plot(frequencies, GR3, 'LineWidth', 1.5);
plot(frequencies, GR4, 'LineWidth', 1.5);
plot(frequencies, GR5, 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Gain (%)', 'FontSize', 20);
grid on;
legend("R = 1 Ω", "R = 10 Ω", "R = 100 Ω", "R = 1000 Ω", "R = 5000 Ω");
hold off

% Calculate gains for different inductor values (R and C fixed)
GL1 = calculateGain(V_in, h, C, L/2, R, V_C0, i_0, k);
GL2 = calculateGain(V_in, h, C, L, R, V_C0, i_0, k);
GL3 = calculateGain(V_in, h, C, L*10, R, V_C0, i_0, k);
GL4 = calculateGain(V_in, h, C, L*100, R, V_C0, i_0, k);
GL5 = calculateGain(V_in, h, C, L*1000, R, V_C0, i_0, k);

subplot(2,2,3);
hold on
title(['5.3: Gain for R = ' num2str(R) ' Ω and C = ' num2str(C) ' F'], FontSize=16);
plot(frequencies, GL1, 'LineWidth', 1.5);
plot(frequencies, GL2, 'LineWidth', 1.5);
plot(frequencies, GL3, 'LineWidth', 1.5);
plot(frequencies, GL4, 'LineWidth', 1.5);
plot(frequencies, GL5, 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('Gain (%)', 'FontSize', 16);
grid on;
legend("L = 50 mH", "L = 100 mH", "L = 1 H", "L = 10 H", "L = 100 H");
hold off

% Calculate gains for different capacitor values (R and L fixed)
GC1 = calculateGain(V_in, h, C/10, L, R, V_C0, i_0, k);
GC2 = calculateGain(V_in, h, C/5, L, R, V_C0, i_0, k);
GC3 = calculateGain(V_in, h, C, L, R, V_C0, i_0, k);
GC4 = calculateGain(V_in, h, 5*C, L, R, V_C0, i_0, k);
GC5 = calculateGain(V_in, h, 50*C, L, R, V_C0, i_0, k);

subplot(2,2,4);
hold on
title(['5.4: Gain for R = ' num2str(R) ' Ω and L = ' num2str(L) ' H'], FontSize=16);
plot(frequencies, GC1, 'LineWidth', 1.5);
plot(frequencies, GC2, 'LineWidth', 1.5);
plot(frequencies, GC3, 'LineWidth', 1.5);
plot(frequencies, GC4, 'LineWidth', 1.5);
plot(frequencies, GC5, 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('Gain (%)', 'FontSize', 16);
grid on;
legend("C = 0.000001 F", "C = 0.000002 F", "C = 0.0001 F", "C = 0.0005 F", "C = 0.005 F");
hold off

clear; clc;

% --- Constants for RLC gain vs frequency plotting ---
R = 1000; L = 100e-3; C = 0.1e-6; Fs = 192000; h = 1/Fs; k_simulation = 0.01;
k = 0:h:k_simulation; V_in = ones(size(k)); V_C0 = 0; i_0 = 0;
frequencies = logspace(log10(100), log10(10000), 9);

% Create figure for gain variation with R, L, C
figure('Name','RLC Gain vs Frequency','Units','inches','Position',[0 0 3.5 6]);
sgtitle('Effect of R, L, and C on RLC Gain','FontSize',12);

% --- Vary resistor R ---
subplot(3,1,1); hold on;
p1 = plot(frequencies, calculateGain2(V_in,h,C,L,1,V_C0,i_0,k),'LineWidth',1.2);
p2 = plot(frequencies, calculateGain2(V_in,h,C,L,10,V_C0,i_0,k),'LineWidth',1.2);
p3 = plot(frequencies, calculateGain2(V_in,h,C,L,100,V_C0,i_0,k),'LineWidth',1.2);
p4 = plot(frequencies, calculateGain2(V_in,h,C,L,1000,V_C0,i_0,k),'LineWidth',1.2);
p5 = plot(frequencies, calculateGain2(V_in,h,C,L,5000,V_C0,i_0,k),'LineWidth',1.2);

set(gca,'XScale','log','FontSize',9);
ylabel('Gain (ratio)','FontSize',10);
title('Varying R','FontSize',10);
grid on;

% Add legend formatting
hLeg = legend([p1 p2 p3 p4 p5], '1 \Omega','10 \Omega','100 \Omega','1 k\Omega','5 k\Omega', ...
   'FontSize',8,'Orientation','horizontal','Location','southoutside','Interpreter','tex');
hLeg.ItemTokenSize = [15,5];

xlabel('Frequency (Hz)','FontSize',10);

% --- Vary inductor L ---
subplot(3,1,2); hold on;
p1 = plot(frequencies, calculateGain2(V_in,h,C,L/2,R,V_C0,i_0,k),'LineWidth',1.2);
p2 = plot(frequencies, calculateGain2(V_in,h,C,L,R,V_C0,i_0,k),'LineWidth',1.2);
p3 = plot(frequencies, calculateGain2(V_in,h,C,10*L,R,V_C0,i_0,k),'LineWidth',1.2);
p4 = plot(frequencies, calculateGain2(V_in,h,C,100*L,R,V_C0,i_0,k),'LineWidth',1.2);
p5 = plot(frequencies, calculateGain2(V_in,h,C,1000*L,R,V_C0,i_0,k),'LineWidth',1.2);

set(gca,'XScale','log','FontSize',9);
ylabel('Gain (ratio)','FontSize',10);
title('Varying L','FontSize',10);
grid on;

hLeg = legend([p1 p2 p3 p4 p5], "50mH","100mH","1H","10H","100H", ...
   'FontSize',8,'Orientation','horizontal','Location','southoutside');
hLeg.ItemTokenSize = [15,5];

xlabel('Frequency (Hz)','FontSize',10);

% --- Vary capacitor C ---
subplot(3,1,3); hold on;
p1 = plot(frequencies, calculateGain2(V_in,h,C/10,L,R,V_C0,i_0,k),'LineWidth',1.2);
p2 = plot(frequencies, calculateGain2(V_in,h,C/5,L,R,V_C0,i_0,k),'LineWidth',1.2);
p3 = plot(frequencies, calculateGain2(V_in,h,C,L,R,V_C0,i_0,k),'LineWidth',1.2);
p4 = plot(frequencies, calculateGain2(V_in,h,5*C,L,R,V_C0,i_0,k),'LineWidth',1.2);
p5 = plot(frequencies, calculateGain2(V_in,h,50*C,L,R,V_C0,i_0,k),'LineWidth',1.2);

set(gca,'XScale','log','FontSize',9);
ylabel('Gain (ratio)','FontSize',10);
title('Varying C','FontSize',10);
grid on;

hLeg = legend([p1 p2 p3 p4 p5], "0.01µF","0.02µF","0.1µF","0.5µF","5µF", ...
   'FontSize',8,'Orientation','horizontal','Location','southoutside');
hLeg.ItemTokenSize = [15,5];

xlabel('Frequency (Hz)','FontSize',10);

%% --- Function Definitions ---

% Calculate gain function variant for plotting gain vs frequency curves
% Inputs:
%   V_in - input voltage vector
%   h - sampling interval
%   C, L, R - capacitor, inductor, and resistor values
%   V_C0, i_0 - initial capacitor voltage and initial current
%   k - time vector
% Output:
%   G - gain vector corresponding to each frequency
function G = calculateGain2(V_in,h,C,L,R,V_C0,i_0,k)
    frequencies = logspace(log10(100), log10(10000), 9);
    G = zeros(1,length(frequencies));
    for i = 1:length(frequencies)
        f = frequencies(i);
        V_in = sin(2*pi*f*k); % sinusoidal input at frequency f
        V_out = simulateRLC(V_in,h,C,L,R,V_C0,i_0);
        if max(abs(V_in))>0
            G(i) = max(abs(V_out))/max(abs(V_in));
        else
            G(i) = 0;
        end
    end
end


% Calculate gain considering output voltage across resistor for bandpass
function G = calculateGain(V_in, h, C, L, R, V_C0, i_0, k)
    frequencies = logspace(log10(100), log10(10000), 9); % Frequencies 100 Hz to 10 kHz
    G = zeros(1, length(frequencies));
    for i = 1:length(frequencies)
        f = frequencies(i);
        V_in = sin(2 * pi * f * k);
        V_out = simulateRLC(V_in, h, C, L, R, V_C0, i_0);
        if max(abs(V_in)) > 0
            G(i) = max(abs(V_out)) / max(abs(V_in));
        else
            G(i) = 0;
        end
    end
end


%% RL Based Figures
clearvars;
% Define sampling time and component values
h = 1/192000;          % Sampling time interval
R = 100;               % Resistance in ohms
L = 0.1;               % Inductance in henries
k_simulation = 0.01;   % Simulation duration
i_0 = 0;               % Initial current

% State-space matrices for discrete-time RL circuit
A = [1 - h*R/L];
B = [h/L];
C = [-R];
D = [1];

% Time vector and input voltage (step)
k = 0:h:k_simulation;
V_in = ones(size(k));

% Create discrete-time state-space system for RL circuit
rl_circuit = ss(A, B, C, D, h);

% Simulate system response to step input voltage
[sysout] = lsim(rl_circuit, V_in, k, i_0);

V_R = sysout; % Voltage across inductor

% Plot input vs voltage across inductor
figure(Name = "RL Circuit: Voltage")
hold on
plot(k, V_in-0.003, LineWidth=1.5); % Input voltage offset slightly for clarity
plot(k, V_R, LineWidth=1.5);
grid on
xlabel("Time (s)")
ylabel("Voltage (V)")
title("RL Circuit Voltage vs Time")
legend("Input Voltage","Voltage across Inductor")
hold off


%% RC Circuit Figures
clearvars;

figure(Name = "RC Circuit")
hold on

% Multiple simulations with varying sampling times to compare capacitor voltage response
for h_val = [0.0009, 0.0007, 0.0005, 0.0003, 0.0001]
    h = h_val;             % Sampling time
    R = 1000;              % Resistance (Ohms)
    C = 1e-6;              % Capacitance (Farads)
    k_simulation = 0.005;  % Total simulation time (seconds)
    V_C0 = 0;              % Initial capacitor voltage
    
    % State-space system matrices for RC circuit
    A = 1 - h/(R*C);
    B = h/(R*C);
    Cmat = [-1; 1];
    D = [1; 0];
    
    k = 0:h:k_simulation;  % Time vector
    V_in = ones(size(k));  % Step input voltage
    
    % Create discrete-time state-space system
    rc_circuit = ss(A, B, Cmat, D, h);
    
    % Simulate response with initial capacitor voltage condition
    [sysout] = lsim(rc_circuit, V_in, k, V_C0);
    
    V_C = sysout(:, 2);    % Capacitor voltage
    plot(k, V_C, LineWidth=1.5);
end

xlabel("Time (s)")
ylabel("Voltage (V)")
title("RC Circuit Capacitor Voltage vs Time")
grid
legend("h = 9.0e-04 s", "h = 7.0e-04 s", "h = 5.0e-04 s", "h = 3.0e-04 s", "h = 1.0e-04 s")

% Cleanup variables keeping only relevant RC variables
clearvars -except rc_circuit V_C V_R;

%% Simulate RLC Function Local Definition
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

% Note: Some comments, figures, and code in this document were generated
% with the assistance of AI. All work was thoroughly reviewed and verified
% before inclusion in our case study. However, this file itself is not
% intended for publication, as it was used solely to generate our figures.