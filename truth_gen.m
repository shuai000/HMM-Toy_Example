function y = truth_gen(p_init,pt,N)

% p: transition matrix
% N: number of samples
y = zeros(1,N);

y(1) = rand_gen(p_init);

for i = 2:N    
    y(i) = rand_gen(pt(y(i-1),:));    
end

end


function x = rand_gen(p)

K = length(p);
num = rand;

accumulation = 0;

for i = 1:K
    accumulation = accumulation + p(i);
    if num <= accumulation
        x = i;
        break;
    end    
end

end