function [alpha,scaling_coefficient] = forward_func(P_init,Pt,likelihood,scaling)

%%% To conduct forward process
%%% input: init_pt pt likelihood

%%% input: scaling  to control if the forward process in sclaed or not
%%% referenced from Rabiner 1989 implementation part
%%% to deal with the underflow of forward variable (alpha) as it is getting smaller
%%% with time

%%% output: alpha forward variable
%%%         scaling_coefficient scaling factor for each time index t

T = size(likelihood,2);	
region_num = length(P_init);
scaling_coefficient = zeros(1,T);

if scaling
    %%% --------  scaling
    % initialization
    alpha_init = P_init'.*likelihood(:,1);
    scaling_coefficient(1) = 1/sum(alpha_init);
    alpha(:,1) = alpha_init.*scaling_coefficient(1);
    % induction
    for t=2:T
        for j=1:region_num
            temp = 0;
            for i=1:region_num
                temp = temp + Pt(i,j)*alpha(i,t-1);
            end
            alpha(j,t) = temp*likelihood(j,t);
        end
        % scaling
        scaling_coefficient(t) = 1/sum(alpha(:,t));
        alpha(:,t) = alpha(:,t).*scaling_coefficient(t);
    end
    
else
    %%% --------  non-scaling
    alpha(:,1) = P_init'.*likelihood(:,1);
    
    for t=2:T
        for j=1:region_num
            temp = 0;
            for i=1:region_num
                temp = temp + Pt(i,j)*alpha(i,t-1);
            end
            alpha(j,t) = temp*likelihood(j,t);
        end
    end
    
end

end