clear
clc
p1 = 0.9; 
p2 = 0.4; 
alpha = 0.3;
n = 20; % number of vertices
NumberOfExperiments = 150;
Record = [];
RecordSim = [];
rng(1958)

for l = 1:NumberOfExperiments
    Simulation;
    RecordSim = [RecordSim; RecordEdges];
    fun = @(x)MomentEq2(x, e1, e2, n)
    options = optimoptions('lsqnonlin','FunctionTolerance',1e-20,'StepTolerance',1.0e-20,'OptimalityTolerance',1e-10,'MaxIterations',20000,'MaxFunctionEvaluations',6000000);
    while 1
        [Est1,f_star,~,flag] = lsqnonlin(fun,[rand rand],[0,0],[1,1],[],[],[],[],[],options); 
        if f_star <10^(-6)
            break;
        end
    end
    fun2 = @(x)mom_theo3(x,e4,min(Est1(1),Est1(2)), max(Est1(1),Est1(2)), n);
    % fun2 = @(x)mom_theo(x,e4,min(Est1(1),Est1(2)), max(Est1(1),Est1(2)), n,2);
    Est2 = fmincon(fun2,rand,[],[],[],[],0,1)
    if Est1(1) > Est1(2)
        Record = [Record; Est1 Est2]
    else
        Record = [Record; Est1(2) Est1(1) Est2]
    end
end
