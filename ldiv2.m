% Long division of polynomials over GF(2)
% Highest degree of numerator and denominator should be made equal through
% zeros.
% NP is numerator polynomial. For 1+D^2, put NP=[1 0 1] and for D^2 put NP=[0 0 1]
% DP is denominator polynomial.
% 'num_terms' is the 'num_terms' first coefficients of the long division of
% NP by DP.
% input and output are strictly row vectors
% Same to the 'ldiv' function in SCILAB but in GF(2).
% written by Vineel Kumar Veludandi

function [quotient] = ldiv2(NP,DP,num_terms)
quotient = zeros(1,num_terms); % quotient 
for coef_cnt = 1:num_terms
quotient(coef_cnt) = NP(1);
temp = bitxor(NP,DP*NP(1));
NP = [temp(2:end) 0]; 
end
end