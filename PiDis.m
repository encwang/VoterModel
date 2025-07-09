p1 = 0.1; 
p2 = 0.9; 
alpha = 0.1;
% n = 10; % number of vertices

A1 = [(1-alpha)* binopdf(0,3,p2), (1-alpha)*binopdf(1,3,p2), (1-alpha)*binopdf(2,3,p2), (1-alpha)*binopdf(3,3,p2); ...
    (1-alpha)* binopdf(0,3,p2), (1-alpha)*binopdf(1,3,p2), (1-alpha)*binopdf(2,3,p2), (1-alpha)*binopdf(3,3,p2); ...
    (1-alpha)* binopdf(0,3,p2), (1-alpha)*binopdf(1,3,p2), (1-alpha)*binopdf(2,3,p2), (1-alpha)*binopdf(3,3,p2); ...
    (1-alpha)* binopdf(0,3,p2), (1-alpha)*binopdf(1,3,p2), (1-alpha)*binopdf(2,3,p2), (1-alpha)*binopdf(3,3,p2)];

B1 = eye(4) * alpha;

a0 = (1- (p1+p2)/2)^2*(1-p2);
a1 = (p1+p2)/2*(1- (p1+p2)/2)*(1-p2) +(p1+p2)/2*(1- (p1+p2)/2)*(1-p2) + (1- (p1+p2)/2)^2*(p2) ;
a2 = (p1+p2)/2*(p1+p2)/2*(1-p2) + (p1+p2)/2*(1- (p1+p2)/2)*(p2) + (p1+p2)/2*(1- (p1+p2)/2)*(p2);
a3 = (p1+p2)/2*(p1+p2)/2*(p2);

A2 = [a0, a1, a2, a3; ...
    a0, a1, a2, a3; ...
    a0, a1, a2, a3; ...
    a0, a1, a2, a3] * (1-alpha);

B21 =  eye(4) * alpha/3; B22 =  eye(4) * alpha*2/3;

b0 = (1- (p1+p2)/2)^2*(1-p1);
b1 = (p1+p2)/2*(1- (p1+p2)/2)*(1-p1) +(p1+p2)/2*(1- (p1+p2)/2)*(1-p1) + (1- (p1+p2)/2)^2*(p1) ;
b2 = (p1+p2)/2*(p1+p2)/2*(1-p1) + (p1+p2)/2*(1- (p1+p2)/2)*(p1) + (p1+p2)/2*(1- (p1+p2)/2)*(p1);
b3 = (p1+p2)/2*(p1+p2)/2*(p1);

A3 = [b0, b1, b2, b3; ...
    b0, b1, b2, b3; ...
    b0, b1, b2, b3; ...
    b0, b1, b2, b3] * (1-alpha);

A4 = [(1-alpha)* binopdf(0,3,p1), (1-alpha)*binopdf(1,3,p1), (1-alpha)*binopdf(2,3,p1), (1-alpha)*binopdf(3,3,p1); ...
    (1-alpha)* binopdf(0,3,p1), (1-alpha)*binopdf(1,3,p1), (1-alpha)*binopdf(2,3,p1), (1-alpha)*binopdf(3,3,p1); ...
    (1-alpha)* binopdf(0,3,p1), (1-alpha)*binopdf(1,3,p1), (1-alpha)*binopdf(2,3,p1), (1-alpha)*binopdf(3,3,p1); ...
    (1-alpha)* binopdf(0,3,p1), (1-alpha)*binopdf(1,3,p1), (1-alpha)*binopdf(2,3,p1), (1-alpha)*binopdf(3,3,p1)];


P = [A1 B1 zeros(4,8); ...
    B21 A2 B22 zeros(4,4); ...
    zeros(4,4) B22 A3 B21; ...
    zeros(4,8) B1 A4];
Pnew = eye(16)-P;
Ppi = [Pnew(:,1:15), ones(16,1)];
b = [zeros(1,15) 1];
Pi = b*inv(Ppi);

[sum(Pi(1:4)), sum(Pi(5:8)), sum(Pi(9:12))  sum(Pi(13:16))]

[sum(Pi(1:4:end)), sum(Pi(2:4:end)), sum(Pi(3:4:end)), sum(Pi(4:4:end))]

C1 = [0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3; Pi(2)*P(2,:)+Pi(6)*P(6,:)+Pi(10)*P(10,:)+Pi(14)*P(14,:)];
C2 = [0 2 4 6 0 2 4 6 0 2 4 6 0 2 4 6; Pi(3)*P(3,:)+Pi(7)*P(7,:)+Pi(11)*P(11,:)+Pi(15)*P(15,:)];
C3 = [0 3 6 9 0 3 6 9 0 3 6 9 0 3 6 9; Pi(4)*P(4,:)+Pi(8)*P(8,:)+Pi(12)*P(12,:)+Pi(16)*P(16,:)];

sum(C1(1,:).*C1(2,:)) + sum(C2(1,:).*C2(2,:)) + sum(C3(1,:).*C3(2,:))

