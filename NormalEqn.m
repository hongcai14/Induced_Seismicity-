function[theta,V_theta] = NormalEqn(X,y,V)
%
%   Takes in three inputs: X, y, and V (if it exists). Outputs three 
%   variables: theta, and V_theta, using an appropriate variation of the
%   normal equation.
%
%   input(s):
%
%   X = a <m x (n + 1)> matrix of features that includes the x_0 feature
%   y = a <m x 1> vector of training examples
%   V = a <m x m> covariance/variance matrix (if it exists), else [] 
%
%   output(s):
%
%   theta = a <(n + 1) x 1> vector of model parameters  
%   V_theta = a < (n + 1) x (n + 1)> matrix of uncertainties(squared) 

    if isempty(V) 
        theta = inv(X'*X)*X'*y;
        %^(theta) ordinary least squares
        
        SSE = ((X*theta - y)'*(X*theta - y))/(size(X,1) - size(X,2) + 1);
        %^sum of squared errors of prediction
        
        V_theta = SSE*inv(X'*X);
        %^(V_theta) ordinary least squares
        
    else
        theta = inv(X'*inv(V)*X)*X'*inv(V)*y;
        %^(theta) weighted/generalized least squares
        
        V_theta = inv(X'*inv(V)*X);
        %^(V_theta) weighted/generalized least squares
    end
end