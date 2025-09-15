function y = ProbOnVertex(n)

Vec = [1]; ind = 1;
for k = 1:n
    ind = ind * (n-k+1);
    Vec = [Vec ind/factorial(k)];
end

y = Vec/sum(Vec);