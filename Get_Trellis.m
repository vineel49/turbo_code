% Trellis for the Convolutional encoder with generator matrix given by
% [1 (1+D^2)/(1+D+D^2)]
function [Prev_State,Prev_Ip_trans,Outputs_prev,Prev_State_trans,Next_State,Outputs_next]= Get_Trellis()
Outputs_prev = [1,4; 2,3; 1,4; 2,3];  % ex: row 1 corresponds to branch metrices that converge to state 1, similarly row 2,3,4
Prev_State = [1,2; 4,3; 2,1; 3,4]; % row 1 corresponds to the previous states corresponding to inputs 0 and 1
Prev_State_trans = [1,4,2,3;2,3,1,4]; % transpose of the P_State matrix
Prev_Ip_trans = [1,1,1,1;2,2,2,2] ; % column 1 corresponds to previous input order.
Next_State = [1,3; 3,1; 4,2; 2,4];
Outputs_next = [1,4; 1,4; 2,3; 2,3];
end
