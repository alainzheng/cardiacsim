% Cardiac simulator based on:
% windkessel model for systemic circulation and 
% the work of Defares (Defares, Osborn, and Hara 1963) for the heart model
% Code Author: Mathieu BIAVA, Laura COUTURE,...
%              Emilien MESSIAEN, Alain Kai Rui ZHENG


%% sources of some parameter value or expression
% * from Benoit HAUT

% ** from "A Dynamical State Space Representation and Performance Analysis of a
% Feedback-Controlled Rotary Left Ventricular Assist Device"

% *** from MIT course

clear all; close all;



%% %%%%%%%%%%%%%
%%%PARAMETERS%%%
%%%%%%%%%%%%%%%%

heartRate = 80; %bpm
cardiacCycle = 60/heartRate; %[s]

%%heart
Emax = 2; % [mmHg/ml] systolic max elastance  (=1/C = 1/0.4) ***
Emin = 0.06; % [mmHg/ml] diastolic min elastance (=1/C = 1/15)  ***
UleftAtria = 8; % [mmHg] left atria pressure (cnst) *
UrightAtria = 3; % [mmHg] right atria pressure (cnst) *
V0 = 15; % dead volume [ml] ***

%%systemic circulation (windkessel) W4series
Rperi = 1; % [mmHg/(ml/s)] total peripheral resistance *
Cart = 2; % [ml/mmHg] total arterial capacitance *
Raorta = 0.0398; % [mmHg/(ml/s)] characteristic resistance of aorta **
Laorta = 0.0005; % [mmHg/(ml/s2)] inertance of blood in aorta **

%%valves
RmitralValve = 0.005; % [mmHg/(ml/s)] mitral valve resistance * 0.01
RaorticValve = 0.1; % [mmHg/(ml/s)] aortic valve resistance * 0.05



%% just helps visualise the variable elastance/compliance of the left
%%% ventricle

% t = 0:0.001:6*cardiacCycle;
% x = elastance(t,cardiacCycle,Emax,Emin);
% figure(1)
% plot(t, x)
% xlabel("time [s]");
% ylabel("elastance [mmHg/ml] and compliance [ml/mmHg]");
% hold on
% plot(t, 1./x)
% hold off



%% real code starts here:

tstart = 0;
tend = ceil(6*cardiacCycle); % we at least n full cycles

%%% initial conditions
Vv0 = 150; % Vv(t) ml 
Up0 = 110; % Up(t) mmHg
Ua0 = 90; %  Ua(t) mmHg

%%% ordinary differential equation (ode) system
[t,y] = ode45( @(t, y) solve_system(t, y, cardiacCycle,...
    Emax, Emin, UleftAtria, UrightAtria, V0,...
    Rperi, Cart, Raorta, Laorta,...
    RmitralValve, RaorticValve),...
    [tstart tend], [Vv0 Up0 Ua0]);
Vv = y(:,1); % = Volume in the left ventricle
Up = y(:,2); % = pressure at the entry of capillaries? en tout cas je pense
Ua = y(:,3); % = pressure in the aorta (normally sin 80-120)

%%% other unknowns that we now know the values of
UleftVent = elastance(t, cardiacCycle, Emax, Emin) .* ((Vv)-V0);

%I is flow in aorta, so the cardiac output but in ml/s instead of ml/min 
%(i think hah)
I=(UleftVent-Ua)./RaorticValve .*...
    heaviside(elastance(t, cardiacCycle, Emax, Emin) .* (Vv-V0) - Ua);

%Imitval is the flow in the mitral valve (so into the left ventricle)
Imitval=(UleftAtria-UleftVent)./RmitralValve.*...
    heaviside(UleftAtria - elastance(t, cardiacCycle, Emax, Emin) .* (Vv-V0));

%Iperi is the flow that is at the input of the periperial resistances
%(Benoit said we can see it is very steady and doesn't oscillate much)
Iperi=(Up-UrightAtria)./Rperi;


figure(2)
subplot(2,1,1)
plot(t, UleftVent)
xlabel("time [s]");
ylabel("Pressure [mmHg]");
hold on
plot(t, Ua)
plot(t, Up)
legend('left ventricule','aorta','peri')
hold off


subplot(2,1,2)
plot(t,I)
xlabel("time [s]");
ylabel("Blood Flow [ml/s]");
hold on
plot(t, Imitval)
plot(t, Iperi)
legend('Aortic','Mitral','peri')
hold off

figure(3)
subplot(2,1,1)
plot(t,Vv)
xlabel("time [s]");
ylabel("Vent. volume [ml]");
subplot(2,1,2)
plot(t, UleftVent)
xlabel("time [s]");
ylabel("Ventricular pressure [mmHg]");



figure(4)
plot(Vv,UleftVent)
xlabel("Ventricular volume [ml]");
ylabel("Ventricular pressure [mmHg]");

%% part where we vary some stuff to mimic some physiological phenomena


%%
function e = elastance (t, cardiacCycle, Emax, Emin)
    % **
    Tmax = 0.2+0.15*cardiacCycle;
    tn = mod(t, cardiacCycle)/Tmax;
    e = (Emax - Emin) *1.55*((tn/0.7).^1.9) ./...
        (1+((tn/0.7).^1.9)) .* (1 ./ (1+(tn/1.17).^21.9))+ Emin;
end

%%
function dedt = delastance(t, cardiacCycle, Emax, Emin)
    % normally derivate the expression above but i just
    % put the equation into an online derivator https://www.wolframalpha.com/
    % and from this
    % 1.55 * (x/0.7)^1.9 / (1+((x/0.7)^1.9)) * (1/(1+(x/1.17)^21.9))
    % obtained its derivative the de/dt
    
    Tmax = 0.2 + 0.15*cardiacCycle;
    tn = mod(t, cardiacCycle)/Tmax;
    dEndt = (1449.81 * tn.^0.9 - 4.693711654661484.^-13 * tn.^2.8 - 490.138 * tn.^22.8 - 1056.93 * tn.^24.7)/((0.507792 + 1* tn.^1.9).^2 * (31.1365 + 1* tn.^21.9).^2 .* Tmax);
    dedt = (Emax - Emin)*dEndt + Emin;
end

%%
function [ydot] = solve_system(t, y, cardiacCycle,...
    Emax, Emin, UleftAtria, UrightAtria, V0,...
    Rperi, Cart, Raorta, Laorta,...
    RmitralValve, RaorticValve)
% holy moly
% y(1) = (Vv) Volume in the left ventricle
% y(2) = (Up) pressure at the entry of capillaries? en tout cas je pense
% y(3) = (Ua) pressure in the aorta (normally sin 80-120)

ydot(1,1) = (( UleftAtria - elastance(t, cardiacCycle, Emax, Emin) * (y(1)-V0))/RmitralValve ) *...
            heaviside(UleftAtria - elastance(t, cardiacCycle, Emax, Emin) * (y(1)-V0)) -...
            (( elastance(t, cardiacCycle, Emax, Emin) * (y(1)-V0) - y(3))/RaorticValve ) *...
            heaviside(elastance(t, cardiacCycle, Emax, Emin) * (y(1)-V0) - y(3));
        
ydot(2,1) = ( (UrightAtria - y(2))/(Rperi*Cart) ) + ...
            (( elastance(t, cardiacCycle, Emax, Emin) * (y(1)-V0) - y(3)) / (RaorticValve*Cart))*...
            heaviside(elastance(t, cardiacCycle, Emax, Emin) * (y(1)-V0) - y(3));
        
ydot(3,1) = -(RaorticValve/Laorta)*(y(3)-y(2))+...
            +(delastance(t, cardiacCycle, Emax, Emin)*(y(1)-V0) + elastance(t, cardiacCycle, Emax, Emin)*ydot(1))*...
            heaviside(elastance(t, cardiacCycle, Emax, Emin) * (y(1)-V0) - y(3))-...
            (Raorta/Laorta)*(elastance(t, cardiacCycle, Emax, Emin) * (y(1)-V0) - y(3))*...
            heaviside(elastance(t, cardiacCycle, Emax, Emin) * (y(1)-V0) - y(3));
end


