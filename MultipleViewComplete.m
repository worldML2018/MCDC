function [X_hat,Objs] = MultipleViewComplete(X, alpha, r)

% Input:
% X: a set of incomplete kernels
% alpha: regularized parameter
% r: the rank of low-rank representation

% Output:
% X_hat: a set of complete kernels
% Objs: objective value of function


if nargin < 2
    alpha  = 0.01 * ones(length(X),1);
end

if nargin < 3
    r = 40;
end




alpha  = 0.01 * ones(length(X),1);
epsilon = 1e-3;
MaxIter = 200;
Iter = 1;
Delta = 10^4;
Deltas = Delta;
 
numView = length(X);
Y= X;

for i = 1: numView
    temp = Y{i};
    index = find(temp== 0);
    temp(index) = 0.001*rand(length(index),1);
    Y{i} = temp;
end
    
    



%initialization
[row,column] = size(X{1});


F = cell(numView,1); % low-rank representation
P = cell(numView,1); % scale matrix

for i = 1:numView
    F{i} = rand(row,r);
    F{i} = F{i}/sqrt(trace(F{i}'*F{i}));
    P{i} = diag(sum(F{i}));
end

F_g = rand(row, r); % common low-rank representation across all view

J1 = obj_fun(Y, F, P, F_g, alpha);
Objs = J1;



% Alternating updatre rule



while Delta > epsilon && Iter <= MaxIter


    
    

	% update F{i}
    for i = 1:numView

    	F{i} = F{i}.* ((2*Y{i}*F{i} + alpha(i)* F_g * P{i}')./ (2*F{i}*F{i}'*F{i} + alpha(i)*F{i}*P{i}*P{i}')).^(0.5);
        P{i} = diag(sum(F{i}));

    end

    % update F_g
    num = 0;
	denom = 0;
	for i = 1: numView
		num = num + alpha(i)*F{i}*P{i};
		denom = denom + alpha(i)*F_g;
	end

    F_g = F_g.*(num./denom).^(0.5);

    
    J2 = obj_fun(Y, F, P, F_g, alpha);
    Delta = J1 -J2;
    Objs = [Objs, J2];
    Deltas = [Deltas, Delta];
    J1 = J2;
    Iter = Iter + 1;

end

for i = 1:numView
    Omega = X{i};
    Omega = logical(Omega);
	temp1 = F{i}*F{i}';
    temp2 = X{i};
    temp1(Omega) = temp2(Omega);
    X_hat{i} = temp1;
end

end


%% Objective function

function J = obj_fun(X, F, P, F_g, alpha)


J = 0;
numView = length(X);
for i = 1:numView
    J = J + norm(X{i}-F{i}*F{i}','fro')^2 +alpha(i) * norm(F{i}*P{i} - F_g,'fro')^2;
end

end
