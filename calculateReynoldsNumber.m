function calculateReynoldsNumber ()
prompt='Enter Temperature [°C] between 20-25:';
temperature = input(prompt);
if temperature==20
  kinViscosity = 1.0034; % [mm²/s]
elseif temperature==21
  kinViscosity = 0.9795; % [mm²/s]
elseif temperature==22
  kinViscosity = 0.9565; % [mm²/s]
elseif temperature==23
  kinViscosity = 0.9344; % [mm²/s]
elseif temperature==24
  kinViscosity = 0.9131; % [mm²/s]
elseif temperature==25
  kinViscosity = 0.8926; % [mm²/s]
end

load('sortTimeVoltage.mat');

D = 0.008; % inner diameter of hose in [m]
Re = sortTimeVoltage(:,3)*D/(kinViscosity/1000000);
% Re ... Reynolds number
% D  ... hydraulic diameter
%        For a circular tube this is simply the diameter of the tube

figure('Name','Reynoldsnumber vs. Voltage');
plot(sortTimeVoltage(:,2),Re(:));
xlabel('Supply voltage of pump [V]');
ylabel('Reynolds number');
xlim([min(sortTimeVoltage(:,2)) max(sortTimeVoltage(:,2))]);
ylim([min(Re) max(Re)]);
hold on
plot([min(sortTimeVoltage(:,2)) max(sortTimeVoltage(:,2))],[1800 1800]);
plot([min(sortTimeVoltage(:,2)) max(sortTimeVoltage(:,2))],[2100 2100]);
hold off

% Re < 1800             laminar flow
% 1800 < Re < 2100      transition range
% Re > 2100             turbulent flow
% Re == 2040            critical value

end