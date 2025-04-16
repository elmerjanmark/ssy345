clc;  clear all ;  close all; 
%%
A = [1 -0.5; 0 1];
mu_q = [2;10];
sigma_q = [1 0; 0 8];


mu_z = A*mu_q;
sigma_z = A*sigma_q*A';

xy = sigmaEllipse2D(mu_z,sigma_z);

plot(xy(1,:), xy(2,:))

% %% How to prepare and save a nice figure
% figure('Position',[300  300  600  400]); % The size of the figure can be changed with the last two values of the position vector.
% plot(x,y) ;
% xlabel('Function input x') % Set x-axis label
% ylabel('Function output y') % Set y-axis label
% legend('y = 3x+2','Location','southeast') % Dont forget to add a legend. You can change its location by changing 'southeast'.
% xlim([-12 12]) % Adjust the xaxis if needed.
% ylim([-30 35]) % Adjust the yaxis if needed.
% grid on ; % Grids almost for all cases should be on
% print('Q1_result1.eps','-depsc') ; % Always save it as an .eps file so that it looks nice on your report!
% 
% %% How to export your source code as .txt file. 
% filename = fullfile('mainExample.m');  % You should change "mainExample.m" with the name of your source file!
% copyfile(filename,'mainExample.txt','f') % Here, 'mainExample.txt' is the output. You should upload the 'main.txt' (or whatever you name it). 
% 

