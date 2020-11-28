function demo(true_seq,mea_seq,detection)

N = length(true_seq);

%%% para

r = .3;
center = [0 0; 2 2; -2 2];

style = ['c';'g';'b'];
textdata = ['S_1'; 'S_2'; 'S_3'];

for i = 1:N
    clf
    xlim([-3 3]);
    plot_circle(0,0,r,'k','S_1');   hold on;
    plot_circle(2,2,r,'k','S_2');   hold on;
    plot_circle(-2, 2,r,'k','S_3'); hold on;
    
    s = true_seq(i);
    dets = detection(i);
    plot_circle(center(s,1),center(s,2),r,style(s)); hold on;
    pause(0.5);
    plot_circle(center(dets,1),center(dets,2),r,style(s),textdata(dets,:),'r'); hold on;
    pause(0.5);
end





% count_fig = 1;
% saveas(gcf,strcat('d',num2str(count_fig),'.jpg'));
end

