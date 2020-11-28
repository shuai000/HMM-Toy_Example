%%% ---------------------------------------------------------- %%%
%%%    An example for HMM analysis, proposed by Bill
%%%    Hidden state number is 3 and measurement is {0,1}
%%%    19/05/2018 Shuai Sun
%%% ---------------------------------------------------------- %%%


clc; clear all; close all;

% parameter set up
N = 3;
ps = 1/N .* ones(1,N);
Pt = [ .7 .3 0; .3 .4 .3; .2 .3 .5];      % transition probability
Pe = [0.2 0.5 0.8; 0.8 0.5 0.2]';         % mesurement probability
% the reason to come up with a struct is purely due to the requirement of 
% vitebi function I wrote
hmm_para.log_Pt = log(Pt);    
hmm_para.P_init = ps;

true_seq = [3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 ];

T = length(true_seq);

yt = zeros(1,T);
likelihood = zeros(N,T);

for i = 1:T
    
    state_cur = true_seq(i);
    
    yt(i) = binornd(1,Pe(state_cur,2)); % probability to generate 1
    
    for k = 1:N
        likelihood(k,i) = Pe(k,(yt(i)+1));
    end
    
end
log_likelihood = log(likelihood);
%%% hmm based forward process

[alpha ,scaling_factor] = forward_func(ps,Pt,likelihood,1); % 1 is to scale the result
% backward
beta = backward_func(Pt,likelihood,scaling_factor);

forward_backward = alpha.*beta;

[max1,index_forward] = max(alpha);

[max2,index_forback] = max(forward_backward);

[v_seq,data ] = viterbi_com_withlikelihood(hmm_para,log_likelihood);

c = 0.2;
q = 0.3;
cb = 0.95;
r1 = ones(1,T);
r2 = 2*ones(1,T);
r3 = 3*ones(1,T);

for k=1:T
    figure(k);
    for i=1:k
        if i == k
            subplot(3,1,1);
            axis([0,30,0,4]);
            plot(r1,'r-.'); hold on; plot(r2,'b-.'); hold on; plot(r3,'g-.'); hold on;
            xp = [i-c,i+c,i+c,i-c,i-c];
            yp = [true_seq(i)-q, true_seq(i)-q, true_seq(i)+q,true_seq(i)+q, true_seq(i)-q];
            patch(xp,yp,'m'); hold on;
            set(gca,'ytick',[1,2,3]);
            set(gca,'YTickLabel',{'S_1','S_2','S_3'});
            set(gca,'color',[cb cb cb]);
            title('Ground Truth');
            
            subplot(3,1,2);
            axis([0,30,0,4]);
            plot(r1,'r-.'); hold on; plot(r2,'b-.'); hold on; plot(r3,'g-.'); hold on;
            xp = [i-c,i+c,i+c,i-c,i-c];
            forward_p = [index_forward(i)-q,index_forward(i)-q,index_forward(i)+q,index_forward(i)+q,index_forward(i)-q];
            patch(xp,forward_p,'b'); hold on;
            set(gca,'ytick',[1,2,3]);
            set(gca,'YTickLabel',{'S_1','S_2','S_3'});
            set(gca,'color',[cb cb cb]);
            title('Forward Algorithm');
            
            subplot(3,1,3);
            axis([0,30,0,4]);
            plot(r1,'r-.'); hold on; plot(r2,'b-.'); hold on; plot(r3,'g-.'); hold on;
            xp = [i-c,i+c,i+c,i-c,i-c];
            vitebi_p = [v_seq(i)-q,v_seq(i)-q,v_seq(i)+q,v_seq(i)+q,v_seq(i)-q];
            patch(xp,vitebi_p,'r'); hold on;
            set(gca,'ytick',[1,2,3]);
            set(gca,'YTickLabel',{'S_1','S_2','S_3'});
            set(gca,'color',[cb cb cb]);
            title('viterbi Algorithm');
        else
            subplot(3,1,1);
            axis([0,30,0,4]);
            plot(r1,'r-.'); hold on; plot(r2,'b-.'); hold on; plot(r3,'g-.'); hold on;
            xp = [i-c,i+c,i+c,i-c,i-c];
            yp = [true_seq(i)-q, true_seq(i)-q, true_seq(i)+q,true_seq(i)+q, true_seq(i)-q];
            patch(xp,yp,[.5,.5,.5]); hold on;
            set(gca,'ytick',[1,2,3]);
            set(gca,'YTickLabel',{'S_1','S_2','S_3'});
            set(gca,'color',[cb cb cb]);
            title('Ground Truth');
            
            subplot(3,1,2);
            axis([0,30,0,4]);
            plot(r1,'r-.'); hold on; plot(r2,'b-.'); hold on; plot(r3,'g-.'); hold on;
            xp = [i-c,i+c,i+c,i-c,i-c];
            forward_p = [index_forward(i)-q,index_forward(i)-q,index_forward(i)+q,index_forward(i)+q,index_forward(i)-q];
            patch(xp,forward_p,[.5,.5,.5]); hold on;
            set(gca,'ytick',[1,2,3]);
            set(gca,'YTickLabel',{'S_1','S_2','S_3'});
            set(gca,'color',[cb cb cb]);
            title('Forward Algorithm');
            
            subplot(3,1,3);
            axis([0,30,0,4]);
            plot(r1,'r-.'); hold on; plot(r2,'b-.'); hold on; plot(r3,'g-.'); hold on;
            xp = [i-c,i+c,i+c,i-c,i-c];
            vitebi_p = [v_seq(i)-q,v_seq(i)-q,v_seq(i)+q,v_seq(i)+q,v_seq(i)-q];
            patch(xp,vitebi_p,[.5,.5,.5]); hold on;
            set(gca,'ytick',[1,2,3]);
            set(gca,'YTickLabel',{'S_1','S_2','S_3'});
            set(gca,'color',[cb cb cb]);
            title('viterbi Algorithm');
        end
    end
    saveas(gcf,strcat('hmm_',num2str(k),'.jpg'));
    pause(0.1);
end
