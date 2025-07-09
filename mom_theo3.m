function [y] = mom_theo3(alpha,mom_emp,a,b,n)
tic;  
%%% n = 20;  
%%% alpha = 0.3; 
%%% a = 0.4;  % \pi_+
%%% b = 0.9;  % \pi_-
% f = (a+b)/2;
f = a*b/(a+b);

% Define the state space dimensions
E1 = n + 1;  % N(t) ranges from 0 to n
E2 = nchoosek(n, 2) + 1;  % Y(t) ranges from 0 to nC2

% Initialize the transition matrix P (size: E1 * E2 x E1 * E2)
P = zeros(E1 * E2);

% Define the transition matrix P based on rules
for i = 0 : n
    %%%for j = 0 : binosafe(n, 2)
        % State (i, j) is at index (i * E2 + j)
        % % % idx1 = i * E2 + j + 1;
        % % % 
        % % % % Transition i -> i+1, i -> i-1 (adjusting for boundary conditions)
        % % % if i < n  % Transition from i to i + 1
        % % %     next_state = (i+1) * E2 + j + 1;
        % % %     P(idx1, next_state) = alpha * (n - i) / n;
        % % % end
        % % % if i > 0  % Transition from i to i - 1
        % % %     next_state = (i-1) * E2 + j + 1;
        % % %     P(idx1, next_state) = alpha * (i) / n;
        % % % end

        % Transition for the second component j -> j with S(t) transition
        maxK1 = 0;
        maxK2 = 0;

        % Check if we can safely call nchoosek
        if i >= 2
            maxK1 = nchoosek(i, 2);
        end
        if n - i >= 2
            maxK2 = nchoosek(n - i, 2);
        end
    %%%end
    for l = 0 : nchoosek(n, 2)
            sum_term = 0;
            for k1 = 0:min(l, maxK1)
                for k2 = 0:min(l-k1, maxK2)
                    if k1 <= maxK1 && k2 <= maxK2
                        % Ensure the arguments for binosafe are valid (non-negative integers)
                        term1 = max(0, i*(n-i));
                        term2 = l - k1 - k2;
                        if term2 >= 0 && term1 >= term2
                            sum_term = sum_term + binosafe(max(0, maxK1), k1) * binosafe(max(0, maxK2), k2) * ...
                                binosafe(term1, term2) * a^k1 * (1-a)^max(0, (maxK1-k1)) * ...
                                b^k2 * (1-b)^max(0, (maxK2-k2)) * (f)^(l-k1-k2) * ...
                                (1 - f)^max(0, (i*(n-i)-l+k1+k2));
                        end
                    end
                end
            end
            P(i * E2 + 1, i * E2 + l + 1) = (1 - alpha) * sum_term;
    end
    % P(idx1, i * E2 + l + 1) = (1 - alpha) * sum_term;
end



for s = 1 : n + 1
    if s <= n
         P((s - 1) * (nchoosek(n, 2) + 1) + 1 : s * (nchoosek(n, 2) + 1), s * (nchoosek(n, 2) + 1) + 1 : (s + 1) * (nchoosek(n, 2) + 1)) = ((alpha * (n - (s - 1))) / n) * eye(nchoosek(n, 2) + 1);
    P(s * (nchoosek(n, 2) + 1) + 1 : (s + 1) * (nchoosek(n, 2) + 1), (s - 1) * (nchoosek(n, 2) + 1) + 1 : s * (nchoosek(n, 2) + 1)) = ((alpha * s) / n) * eye(nchoosek(n, 2) + 1);
    end
    P((s - 1) * (nchoosek(n, 2) + 1) + 2 : s * (nchoosek(n, 2) + 1), (s - 1) * (nchoosek(n, 2) + 1) + 1 : s * (nchoosek(n, 2) + 1)) = repmat(P((s-1) * (nchoosek(n, 2) + 1) + 1, (s - 1) * (nchoosek(n, 2) + 1) + 1 : s * (nchoosek(n, 2) + 1)), nchoosek(n, 2), 1);
end

% Normalize rows of P to ensure each row sums to 1 (for valid transition matrix)
P = P ./ sum(P, 2);

% Display the transition matrix
%disp('Transition Matrix P:');
%disp(P);

A0 = cell(1, n + 1);
A_1 = cell(1, n);
A1 = cell(1, n);
 for i = 1 : n + 1 
     A0{1,i} = P(1 + (i-1) * (nchoosek(n, 2) + 1) : i * (nchoosek(n, 2) + 1), 1 + (i-1) * (nchoosek(n, 2) + 1) : i * (nchoosek(n, 2) + 1)); 
     if i <= n
     A_1{1, i} = P(1 + i * (nchoosek(n, 2) + 1): (i+1) * (nchoosek(n, 2) + 1), 1 + (i-1) * (nchoosek(n, 2) + 1) : i * (nchoosek(n, 2) + 1));
     A1{1, i} = P(1 + (i-1) * (nchoosek(n, 2) + 1) : i * (nchoosek(n, 2) + 1), 1 + i * (nchoosek(n, 2) + 1) : (i+1) * (nchoosek(n, 2) + 1));
     end
 end


U = cell(1, n);
R = cell(1, n);
G = cell(1, n);
 
U{n} = A0{n+1};
R{n} = (A1{n})/(eye((nchoosek(n, 2) + 1))-U{n});
G{n} = (eye((nchoosek(n, 2) + 1))-U{n})\(A_1{n});
 
for k = n-1:-1:1
    U{k} = A0{k+1} + A1{k+1}*G{k+1};
    R{k} = (A1{k})/(eye((nchoosek(n, 2) + 1))-U{k});
    G{k} = (eye((nchoosek(n, 2) + 1))-U{k})\(A_1{k});
end


H = A0{1} + A1{1} * G{1};

L = ones(nchoosek(n, 2) + 1, 1);
for i = 1 : n
    S = R{1};
    for j = 2 : i
        S = S * R{j};
     %L = ones(n+1, n) + R{1} * ones(n+1,1) + R{1}*R{2}*ones(n+1,1) + R{1}*R{2}*R{3}*ones(n+1,1)
    end
     L = L + S * ones((nchoosek(n, 2) + 1),1);
end

A = [H' - eye((nchoosek(n, 2) + 1)); ones(1,(nchoosek(n, 2) + 1))];
c = [zeros((nchoosek(n, 2) + 1),1); 1];

nu = cell(1, n + 1);

nu{1} = ((A \ c) ./ L)';

for i = 2 : n + 1 

nu{i} = nu{i-1} * R{i-1};

end

nu = [nu{:}];


 % % % % % 3. Cross Moment E[S(t)S(t+1)]: \sum_{i,j,\ell,r} \ell r P_{(i,\ell) (j,r)} \pi(i,\ell)
 % % % % third_moment = 0;
 % % % % for i = 0:n
 % % % %     for k = 0:nchoosek(n,2)
 % % % %         idx1 = i * E2 + k + 1;
 % % % %         for j = 0:n
 % % % %             for r = 0:nchoosek(n,2)
 % % % %                 idx2 = j * E2 + r + 1;
 % % % %                 third_moment = third_moment + r * k * (P(idx1, idx2))^5 * nu(idx1);
 % % % %             end
 % % % %         end
 % % % %     end
 % % % % end

 % 3. Cross Moment E[S(t)S(t+1)]: \sum_{i,j,\ell,r} \ell r P_{(i,\ell) (j,r)} \pi(i,\ell)
 third_moment = 0;
 for i = 0:n
     for k = 0:nchoosek(n,2)
         idx1 = i * E2 + k + 1;
         for j = 0:n
             for r = 0:nchoosek(n,2)
                 idx2 = j * E2 + r + 1;
                         third_moment = third_moment + r * k * P(idx1, idx2) * nu(idx1);
                     end
                 end
             end
 end

  E_S2 = 0;
 for i = 0:n
     for k = 1:nchoosek(n,2)+1
         E_S2 = E_S2 + (k-1)^2 * nu(i * E2 + k);
     end
 end


 moment = 2 * E_S2 - 2 * third_moment;

y = (moment - mom_emp)^2;

elapsedTime = toc;  
% fprintf('Simulation time: %.2f minutes\n', elapsedTime / 60);
end

function C = binosafe(n, k)
    C = exp(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1));
end