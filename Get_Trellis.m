% Trellis for the Convolutional encoder with generator matrix given by
% [1 (1+D^2)/(1+D+D^2)]
function [P_State,P_Ip_trans,Ga_Inx,P_State_trans,N_State,Gb_Inx]= Get_Trellis()
Ga_Inx = [1,4; 2,3; 1,4; 2,3];  % ex: row 1 corresponds to branch metrices that converge to state 1, similarly row 2,3,4
P_State = [1,2; 4,3; 2,1; 3,4]; % row 1 corresponds to the previous states corresponding to inputs 0 and 1
P_State_trans = [1,4,2,3;2,3,1,4]; % transpose of the P_State matrix
P_Ip_trans = [1,1,1,1;2,2,2,2] ; % column 1 corresponds to previous input order.
N_State = [1,3; 3,1; 4,2; 2,4];
Gb_Inx = [1,4; 1,4; 2,3; 2,3];
end
