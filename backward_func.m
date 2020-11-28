function beta = backward_func(Pt,likelihood,scaling_coefficient)

%%% To conduct backward process
%%% input: pt likelihood scaling_coefficient
%%%        scaling_coefficient is calculated from backward process

%%% input: scaling  to control if the forward process in sclaed or not
%%% referenced from Rabiner 1989 implementation part
%%% to deal with the underflow of forward variable (alpha) as it is getting smaller
%%% with time

%%% output: beta  backward variable

[region_num,T] = size(likelihood);
beta = zeros(region_num,T);

if nargin == 2     % no scalling required
    beta(:,T) = 1;    % initialization
    for t=T-1:-1:1
        for i=1:region_num
            temp = 0;
            for j=1:region_num
                temp = temp + Pt(i,j)*likelihood(j,t+1)*beta(j,t+1);
            end
            beta(i,t) = temp;
        end
    end
else
    if nargin == 3 % scaling is required
        beta(:,T) = 1*scaling_coefficient(T);    % initialization
        for t=T-1:-1:1
            for i=1:region_num
                temp = 0;
                for j=1:region_num
                    temp = temp + Pt(i,j)*likelihood(j,t+1)*beta(j,t+1);
                end
                beta(i,t) = temp*scaling_coefficient(t);
            end
        end
    end
end

end