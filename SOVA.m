% Bidirectional Soft Output Viterbi Algorithm 
% Notation: if soft_output >0, then decoded bit is 0.
% outputs extrinsic information
% Written by Vineel Kumar Veludandi

function[soft_output]=SOVA(apr_LLR,num_sym,branch_metric)
C = log(exp(-apr_LLR/2)./(1+exp(-apr_LLR)));

[Prev_State,Prev_Ip_trans,Outputs_prev,Prev_State_trans,Next_State,Outputs_next]= Get_Trellis();
num_states = 4;  % number of states 
soft_output = zeros(1,num_sym); % soft output
survivor_node = zeros(num_states,num_sym); % survivor nodes
survivor_ip = zeros(num_states,num_sym); % survivor inputs
F_path_metric = zeros(num_states,num_sym+1); % forward path metrics
B_path_metric = zeros(num_states,num_sym+1); % backward path metrics
index_temp = [0;1*2;2*2;3*2]; %for linear indexing.

for sym_cnt=  1:num_sym  
   [F_path_metric(:,sym_cnt+1),index] = min([F_path_metric(Prev_State(:,1),sym_cnt)+ branch_metric(Outputs_prev(:,1),sym_cnt)-(C(sym_cnt)+apr_LLR(sym_cnt)/2) ...
       F_path_metric(Prev_State(:,2),sym_cnt)+ branch_metric(Outputs_prev(:,2),sym_cnt)-((C(sym_cnt)-apr_LLR(sym_cnt)/2))],[],2);
    
   [B_path_metric(:,num_sym+1-sym_cnt),~] = min([B_path_metric(Next_State(:,1),num_sym+2-sym_cnt)+ branch_metric(Outputs_next(:,1),num_sym+1-sym_cnt)-(C(num_sym+1-sym_cnt)+apr_LLR(num_sym+1-sym_cnt)/2) ...
       B_path_metric(Next_State(:,2),num_sym+2-sym_cnt)+ branch_metric(Outputs_next(:,2),num_sym+1-sym_cnt)-(C(num_sym+1-sym_cnt)-apr_LLR(num_sym+1-sym_cnt)/2)],[],2);
   
   survivor_node(:,sym_cnt) = Prev_State_trans(index+index_temp);
   survivor_ip(:,sym_cnt) = Prev_Ip_trans(index+index_temp); 

end 

[ml_metric,trace_bf] = min(F_path_metric(:,num_sym+1));

for bk_cnt= num_sym:-1:1
ip = survivor_ip(trace_bf,bk_cnt);
com_ip = bitxor(ip-1,1)+1; % complementary input
com_metric = min(F_path_metric(:,bk_cnt)+branch_metric(Outputs_next(:,com_ip),bk_cnt) + B_path_metric(Next_State(:,com_ip),bk_cnt+1) );
soft_output(bk_cnt) = ((-1)^(ip+1))*(com_metric-ml_metric); 
trace_bf = survivor_node(trace_bf,bk_cnt);    
end

% normalizing (to avoid numerical instabilities)
soft_output(soft_output>50) = 50;
soft_output(soft_output<-50) = -50;


end % for function
