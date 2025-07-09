clear
clc
Record1 = [];
Record2 = [];
p1 = 0.9; 
p2 = 0.4; 
n = 10;
for k = 1:500
    %%%%%%%%
    alpha = 0.3;
    Simulation
    Record1 = [Record1; e1 e2 e4];
    %%%%%%%%
    alpha = 0.6;
    Simulation
    Record2 = [Record2; e1 e2 e4];
    %%%%%%%%
end
% save('ModelS')