p1 = 0.4; 
p2 = 0.9; 
alpha = 0.3;
n = 10; % number of vertices
indRecord = [];

% f = (p1+p2)/2;
f = p1*p2/(p1+p2);

K = 10^5+200; % number of observations in each run
RecordVertex = [];
RecordEdges = [];

VertexCurrent = binornd(n,0.5); % record the number of + vertices at time 0
PairOn = NchooseK(VertexCurrent,2);
PairOff = NchooseK(n-VertexCurrent,2); % n - VertexCurrent is the number of - vertices at time 0
PairOnOff = NchooseK(VertexCurrent,1)*NchooseK(n-VertexCurrent,1);
% PairOn = NumberOfOn * (NumberOfOn-1)/2;
% PairOff = (n-NumberOfOn) * (n-NumberOfOn-1)/2;
% PairOnOff = NumberOfOn * (n-NumberOfOn);
EdgesCurrent = binornd(PairOn,p1) + binornd(PairOff,p2) + binornd(PairOnOff,f);
%%%%%%%%%%%%%%%%%%%%%%%%
RecordVertex = [RecordVertex VertexCurrent];
RecordEdges = [RecordEdges EdgesCurrent];
%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:K
    if rand <= alpha %change opinion
        indRecord = [indRecord k];
        if rand <= (VertexCurrent/n) % + opion is chosen
            VertexCurrent = VertexCurrent - 1;
        else
            VertexCurrent = VertexCurrent +1;
        end
    else % sample links
        PairOn = NchooseK(VertexCurrent,2);
        PairOff = NchooseK(n-VertexCurrent,2);
        PairOnOff = NchooseK(VertexCurrent,1)*NchooseK(n-VertexCurrent,1);
        % PairOn = VertexCurrent * (VertexCurrent-1)/2;
        % PairOff = (n-VertexCurrent) * (n-VertexCurrent-1)/2;
        % PairOnOff = VertexCurrent * (n-VertexCurrent);
        EdgesCurrent = binornd(PairOn,p1) + binornd(PairOff,p2) + binornd(PairOnOff,f);
    end
RecordVertex = [RecordVertex VertexCurrent];
RecordEdges = [RecordEdges EdgesCurrent];
end

% RecordEdges = RecordEdges(5001:end); K = length(RecordEdges);
%%%%%%%%%%%%%%%%%%%
RecordEdges = RecordEdges(200:end);
K = length(RecordEdges);
e1 = mean(RecordEdges);
%%%%%%%%%%%%%%%%%%%
e2Record = [];
for k = 1:K
    e2Record = [e2Record, RecordEdges(k)^2];
end
e2 = mean(e2Record);
%%%%%%%%%%%%%%%%%%%
e3Record = [];
for k = 1:K-1
    e3Record = [e3Record, RecordEdges(k)*RecordEdges(k+1)];
end
e3 = mean(e3Record) - e1^2; 
%%%%%%%%%%%%%%%%%%%
e4Record = [];
for k = 1:K-1
    e4Record = [e4Record, (RecordEdges(k+1)-RecordEdges(k))^2];
end
e4 = mean(e4Record); 

% e5Record = [];
% for k = 1:K-1
%     e5Record = [e5Record, RecordEdges(k+1)-RecordEdges(k)];
% end
% e5 = mean(e5Record); 