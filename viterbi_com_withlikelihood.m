function [max_seq,process_data] = viterbi_com_withlikelihood(hmm_para,log_likelihood)
%%% HMM structure, apply viterbi algorithm to inference hidden sequence
%%% with maximum likelihood
%%% Shuai Sun 25/01/2018
%%% Input: mea_seq: measurement sequence
%%%        hmm_para: structure including necessary HMM parameter
%%% Output:
%%%        max_seq: max sequence by viterbi algorithm
%%%        process_data: store some useful intermidiate result

T = size(log_likelihood,2);
% mean_u = hmm_para.mean;
% cov_u = hmm_para.cov;
Pt_log = hmm_para.log_Pt;
P_init = hmm_para.P_init;
N = size(Pt_log,2); % region_num

phi_s = zeros(N,T);
viterbi_seq = zeros(1,T);
index_record = zeros(N,T);


%%% ------------ viterbi algorithm
% Initialization

for i=1:N
    phi_s(i,1) = log(P_init(i)) + log_likelihood(i,1);
end
% Recursion
for t=2:T    
    for i=1:N
        [phi_s(i,t),index_record(i,t)] = max_path_finding(phi_s(:,t-1),i,log_likelihood(i,t),Pt_log);
    end
end
% Termination
[~,viterbi_seq(T)] = max(phi_s(:,T));
% Path Backtracking
for t=T-1:-1:1
    viterbi_seq(t) = index_record(viterbi_seq(t+1),t+1);
end

process_data.phi_s = phi_s;
process_data.loglikelihood = log_likelihood;
process_data.index = index_record;

max_seq = viterbi_seq;
end

function [s,index] = max_path_finding(phi_s,j,log_like,log_Pt)
index = 1;
s_max = phi_s(1) + log_Pt(1,j);

for i=2:length(phi_s)
    temp = phi_s(i) + log_Pt(i,j);
    if temp > s_max
        index = i;
        s_max = temp;
    end
end

s = s_max + log_like;
end