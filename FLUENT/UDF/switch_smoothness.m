% Matlab script to depict how the *_smooth coefficient works towards
% changing the sensor region definition
%
% Michael Emory 4/5/2012
clear all; clc;

% smoothnes coefficient
SMOOTH = [010 020 050 080 200];

% vector of wall distance (maximum is 0.5h in square duct)
walld = 0:.01:0.5;

% user defined cut-off value
cutoff = 0.3;

myswitch = zeros(length(walld),length(SMOOTH));

for i=1:length(SMOOTH)
    
    % computation of *_switch function
    my_z = SMOOTH(i) * (walld - cutoff);

    % tanh() computation
    myswitch(:,i) = max(min(1.0, 0.5 - 0.5*tanh(my_z)),0.0);

    % plot profile
%     figure(1);hold on;
%     plot(walld,myswitch(:,i),'.-','MarkerSize',20,'LineWidth',2);
%     hold off;
end

figure(1);hold on;
plot(walld,myswitch,'.-','MarkerSize',15,'LineWidth',2);
M=[];
for i=1:length(SMOOTH)
   M = cat(1,M,num2str(SMOOTH(i),'%03.f'));
end
legend(M);
xlabel('d_{wall}/h');ylabel('f_{switch}');
axis([0 0.5 -.1 1.1]);grid on;
title('Influence of Smoothing Coefficient');
hold off;