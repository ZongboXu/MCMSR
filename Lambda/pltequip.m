function pltequip(data,vp,vr,f)
% plot Rayleigh/P energy ratio changing with time
% data: equip3.dat from mcpray(fs).c
% vp: P-wave velocity, 5km/s 
% vr: Rayleigh-wave phase/group velocity, 2.7425 km/s
% f: frequency, 1Hz mcpray(fs).c

n = 2;% scatterer density 
sigR = 1.7154E-2;%cross section of Rayleigh due to delta lambda 
imft = sigR*n*vr; % inverse of mean free time

% plot
plot(data(:,1)*imft,data(:,5),'k','Linewidth',2);
hold on;
% the theoretical equipartition energy ratio
plot(1:1500,ones(1500,1)*pi*vp^3/f/2/pi/vr/vr,'r');
xlim([0 floor(data(end,1)*imft)]);
ylim([5 500]);
set(gca,'YScale','log');
xlabel('Time ($\tau_R$)','Interpreter','latex')
ylabel('$\overline{\overline{E_s}}/\overline{E_b}$','Interpreter','latex');
title('1Hz,10km-depth source, no free surface');
legend('Simulation','Theory');
end