%% Computational Fluid Dynamics 2015 | Coursework 1: Section 4
%NAME: Comparison of phase and diffusion errors with varying frequency
%GOAL: For the advection-diffusion equation, using an implicit approximation centered in space for the diffusion and advection terms and first order in time, the behavior of the phase and diffusion errors with frequency will be plotted. For various values of the Courant (sigma) and Diffusive term (beta) numbers.
%Created by: SERGIU PETRE ILIEV. Date: 30.XI.2015
%Imperial College London | Aeronautical Engineering | Prof. Sherwin Spencer | AE3-414 Computational Fluid Dynamics | Coursework I

%% Clean and prepare workspace
clear all           % Clear all variables
clc                 % Clear the window
close all           % Close open figures

%% Initialise parameters
n=1800;                             % Number of phase increments to consider (i.e. resolution of the graph)
phi=0:(pi()/(n-1)):pi();            % Define the phase vector between 0 and pi with resolution (n) specified above

sigma=[0.2, 0.5, 1.0, 2.0, 3.0];    % Define the values to be considered for Courant number, the first term in the matrix is the one sigma will be held constant at when varying beta
beta=[0.2, 0.2, 0.3, 0.4, 0.5];    % Define the values to be considered for the Diffusive term, the first term in the matrix is the one beta will be held constant at when varying sigma

e_d=zeros(n,8);                     % Preallocate an array for the dispersion error (rows) for the 8 cases (collumns) to be considered
e_phi=zeros(n,8);                   % Preallocate the phase/diffision error (rows) for the 8 cases (collumns) to be considered

%% Calculate errors using for the advection-diffusion numerical approximation
% Holding sigma constant and varying beta
for j=2:1:5
      e_d(:,j-1)=(((1+4*beta(j)*(sin(phi/2)).^2).^2+((sigma(1).*sin(phi)).^2)).^(-1/2))./(exp(-beta(j)*(phi.^2)));        % Determine the diffusion error when beta is varied - algorith from equation (48) in the report
      e_phi(:,j-1)=(atan((sigma(1).*sin(phi))./(1+4*beta(j)*(sin(phi./2)).^2)))./(sigma(1).*phi);                         % Determine the phase error when beta is varied - algorith from equation (48) in the report
end

% Holding beta constant and varying sigma
for j=2:1:5
      e_d(:,4+j-1)=(((1+4*beta(1)*(sin(phi/2)).^2).^2+((sigma(j).*sin(phi)).^2)).^(-1/2))./(exp(-beta(1)*(phi.^2)));      % Determine the diffusion error when sigma is varied - algorith from equation (48) in the report
      e_phi(:,4+j-1)=(atan((sigma(j).*sin(phi))./(1+4*beta(1)*(sin(phi./2)).^2)))./(sigma(j).*phi);                       % Determine the phase error when sigma is varied - algorith from equation (48) in the report
end

%% Plot diffusion error graphs
figure(1)
hold all; grid on
title('Diffusion Error variation with phase for \sigma=0.2 and various values of the diffusive term \beta')
plot(phi, e_d(:,1),'r:', 'LineWidth',1.5);
plot(phi, e_d(:,2),'m-.', 'LineWidth',1.0);
plot(phi, e_d(:,3),'b--', 'LineWidth',1.0);
plot(phi, e_d(:,4),'k-', 'LineWidth',1.0);
ylabel('$\displaystyle\frac{\mid G \mid}{\mid \tilde{G} \mid}$','interpreter','latex'); % Label the y axis
xlabel('\phi (radians)');                                                      % Label the x axis
legend('\beta=0.2', '\beta=0.3', '\beta=0.4', '\beta=0.5', 'Location','northwest')      % Specify the legend and its position

figure(2)
hold all; grid on
title('Diffusion Error variation with phase for \beta=0.2 and various values of the Courant Number \sigma')
plot(phi, e_d(:,5),'r:', 'LineWidth',1.5);
plot(phi, e_d(:,6),'m-.', 'LineWidth',1.0);
plot(phi, e_d(:,7),'b--', 'LineWidth',1.0);
plot(phi, e_d(:,8),'k-', 'LineWidth',1.0);
ylabel('$\displaystyle\frac{\mid G \mid}{\mid \tilde{G} \mid}$','interpreter','latex'); % Label the y axis
xlabel('\phi (phase in radians)');                                                      % Label the x axis
legend('\sigma=0.5', '\sigma=1.0', '\sigma=2.0', '\sigma=3.0', 'Location','northwest')  % Specify the legend and its positio

%% Plot phase error graphs
figure(3)
hold all; grid on
title('Phase Error variation with phase for \sigma=0.2 and various values of the diffusive term \beta')
plot(phi, e_phi(:,1),'r:', 'LineWidth',1.5);
plot(phi, e_phi(:,2),'m-.', 'LineWidth',1.0);
plot(phi, e_phi(:,3),'b--', 'LineWidth',1.0);
plot(phi, e_phi(:,4),'k-', 'LineWidth',1.0);
ylabel('$\displaystyle\frac{\mid \phi \mid}{\mid \tilde{\phi} \mid}$','interpreter','latex'); % Label the y axis
xlabel('\phi (phase in radians)');                                                      % Label the x axis
legend('\beta=0.2', '\beta=0.3', '\beta=0.4', '\beta=0.5')                              % Specify the legend

figure(4)
hold all; grid on
title('Phase Error variation with phase for \beta=0.2 and various values of the Courant Number \sigma')
plot(phi, e_phi(:,5),'r:', 'LineWidth',1.5);
plot(phi, e_phi(:,6),'m-.', 'LineWidth',1.0);
plot(phi, e_phi(:,7),'b--', 'LineWidth',1.0);
plot(phi, e_phi(:,8),'k-', 'LineWidth',1.0);
ylabel('$\displaystyle\frac{\mid \phi \mid}{\mid \tilde{\phi} \mid}$','interpreter','latex'); % Label the y axis
xlabel('\phi (phase in radians)');                                                      % Label the x axis
legend('\sigma=0.5', '\sigma=1.0', '\sigma=2.0', '\sigma=3.0')                          % Specify the legend