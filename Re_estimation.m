%%% ---------------------------------------------------------- %%%
%%%  Reestimation for HMM parameters
%%%  Ref:  https://en.wikipedia.org/wiki/Baum%E2%80%93Welch_algorithm
%%%        Rabiner: HMM tutorial paper
%%%    20/05/2018 Work for the whole Sunday, Shuai Sun
%%% ---------------------------------------------------------- %%%

clc; clear all; close all;

%%% parameter setting
N = 3;    % number of hidden states
Pt = [ .7 .3 0; .3 .4 .3; .2 .3 .5];      % transition probability
Pe = [0.3 0.5 0.7; 0.7 0.5 0.3]';         % measurement probability {0,1} finite measurement

% initial guess
ps = 1/N .* ones(1,N);
Pt0 = [ .6 .4 0; .35 .4 .25; .25 .3 .45];      
Pe0 = [0.4 0.5 0.6; 0.6 0.5 0.4]';

%%% Simulation data generation
true_seq = truth_gen(ps,Pt,100000); % based on true transition probability Pt

T = length(true_seq);
yt = zeros(1,T);
for i = 1:T    
    state_cur = true_seq(i);    
    yt(i) = binornd(1,Pe(state_cur,2)); % probability to generate     
end

%%% initialization for re-estimation parameters
likelihood = zeros(N,T);
count = 0;
max_iter = 20;

hmm_para.P_init = ps;
hmm_para.Pt = Pt0;
hmm_para.Pe = Pe0;

while count < max_iter
    
    count = count + 1
    % evaluate likelihood
    for i = 1:T        
        for k = 1:N
            likelihood(k,i) = hmm_para.Pe(k,(yt(i)+1));
        end        
    end    
    log_likelihood = log(likelihood);
    
    % forward process
    [forward_seq,scaling_factor] = forward_func(hmm_para.P_init,hmm_para.Pt,likelihood,1); % 1 is to scale the result
    % backward
    backward_seq = backward_func(hmm_para.Pt,likelihood,scaling_factor);    
    % for initial probability estimation
    yita_seq = (forward_seq.*backward_seq)./(sum(forward_seq.*backward_seq));  % Eq(27)
    
    hmm_para.P_init = yita_seq(:,1)';
    
    % estimation for a_ij (Pt) Eq(95)
    A = zeros(N,N);
    for i=1:N
        for j=1:N
            numer = 0;
            deno = 0;
            for t=1:T-1
                numer = numer + forward_seq(i,t)*hmm_para.Pt(i,j)*likelihood(j,t+1)*backward_seq(j,t+1);
                temp = 0;
                for k=1:N
                    temp = temp + hmm_para.Pt(i,k)*likelihood(k,t+1)*backward_seq(k,t+1);
                end
                deno = deno + forward_seq(i,t)*temp;
            end
            A(i,j) = numer/deno;
        end
    end    
    hmm_para.Pt = A;
    
    %%% estimation for bj(k)   Eq (40c)
    B = zeros(N,2);    
    for i=1:N
        for j=1:2  % 2 is measurement dimension
            numer = 0;
            deno = sum(yita_seq(i,:));
            for t=1:T
                if yt(t)+1 == j  % we need to check if the corresponding measurement matches
                    numer = numer + yita_seq(i,t);
                end
            end
            B(i,j) = numer/deno;
        end
    end
    hmm_para.Pe = B;
    %%% calculate total log likelihood
    log_p = 0;
    for t=1:T
        log_p = log_p + (-log(scaling_factor(t)));
    end
    log_p_all(count) = log_p;
end

figure;   % plot log_likelihood
plot(log_p_all,'r-.');
