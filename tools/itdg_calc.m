clear
close all


% dir_dist = 0.3:0.1:10;
% rb_dist = 0.95; % reflective boundary min distance from the speaker - head axis
% ref_dist = 2 * sqrt((dir_dist/2).^2 + rb_dist^2);
% itdg = (ref_dist - dir_dist) * 2.9; % in ms
% itdg_fs = ceil(itdg .* (48000 / 1000));
% 
% plot(dir_dist,itdg_fs)
% ylabel('Initial Time Delay Gap (samples@48kHz)')
% xlabel('speaker - head distance (m)')

% how far is the closest boundary? 
% Standard ceiling height: 240 cm
% Tall peron height: 180 cm
% let's assume 0.6 m


figure
hold on

rb_dist = 0.6:0.1:1.2; % reflective boundary min distance from the speaker - head axis
direct_dist = 0.6:0.1:10; % speaker - head distance vector

c = 343;

for i = 1:length(rb_dist)
   for j = 1:length(direct_dist)
      ref_dist = 2 * sqrt((direct_dist(j)/2)^2 + rb_dist(i)^2);
      itdg(i,j) = (ref_dist - direct_dist(j)) * (1000/c); % in ms
   end
end

for i = 1:length(rb_dist)
    plot(direct_dist,itdg(i,:))
    text(0.3,itdg(i,1)+0.2,[num2str(rb_dist(i)) ' m'])
%     x = [0.2/10  (itdg(i,1)+0.2)/7];
%     y = [1.2/10 (itdg(i,1)+1.2)/7];
%     annotation('textarrow',x,y,'String',[num2str(rb_dist(i)) ' m'])
end

xlim([0 4])
ylim([0 7])
ylabel('Initial Time Delay Gap (ms)')
xlabel('Speaker - head distance d_d (m)')

text(0.17,3.1,'Room boundary distance d_b','rotation',90)

ff_limit = 1.5;
xline(ff_limit,'--k')
text(ff_limit+0.2,5,'Far field boundary','rotation',90)

win_length = (140/48000) * 1000;
yline(win_length,'--k')
text(2.8,win_length+0.3,'Time window length')

% save figure
figlen = 4;
width = 4*figlen;
height = 3*figlen;
set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
saveas(gcf,'itdg_calc.pdf')




