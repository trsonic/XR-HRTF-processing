function plotAzElM(az,el,val,lim,clbl,texta,textb)
    step = 0.5; % 0.5
    azimuth = -180:step:180;
    elevation = -90:step:90;
    azimuth = azimuth + randn(size(azimuth))*0.01;
    elevation = elevation + randn(size(elevation))*0.01;
    azimuth = min(max(azimuth,-180),180);
    elevation = min(max(elevation,-90),90);

    for i = 1:length(azimuth)
        for j = 1:length(elevation)
            dist = distance(elevation(j),azimuth(i),el,az);
            idx = find(dist <= min(dist) + 0.001);
            value(j,i) = mean(val(idx));
        end    
    end
    
    axesm('MapProjection','mollweid','MapLatLimit',[-90 90],'Gcolor','black','GLineWidth',1.0,'MLineLocation',45,'PLineLocation',30); 
    axis off
    gridm('on')

    % plot data
    s = pcolorm(elevation,azimuth,value);
    s.EdgeColor = 'none';

    % color bar
    caxis(lim)
%     c=colorbar('location','southoutside','position',[0.2 0.0 0.6 0.05],'box','on','color',[0 0 0]); %[left, bottom, width, height]
    c=colorbar('location','southoutside');
    c.Label.String=clbl;
    c.Label.FontSize=14;
    c.FontSize=14;

    % labels
    fsize = 16;
    fcolor = 'black';
    vertshift = -5;
    horishift = 5;
    textm(vertshift,horishift-135,'-135','color',fcolor,'fontsize',fsize);
    textm(vertshift,horishift-90,'-90','color',fcolor,'fontsize',fsize);
    textm(vertshift,horishift-45,'-45','color',fcolor,'fontsize',fsize);
    textm(vertshift,horishift+0,'0','color',fcolor,'fontsize',fsize);
    textm(vertshift,horishift+45,'45','color',fcolor,'fontsize',fsize);
    textm(vertshift,horishift+90,'90','color',fcolor,'fontsize',fsize);
    textm(vertshift,horishift+135,'135','color',fcolor,'fontsize',fsize);
    textm(vertshift-60,horishift+0,'-60','color',fcolor,'fontsize',fsize);
    textm(vertshift-30,horishift+0,'-30','color',fcolor,'fontsize',fsize);
    textm(vertshift+30,horishift+0,'30','color',fcolor,'fontsize',fsize);
    textm(vertshift+60,horishift+0,'60','color',fcolor,'fontsize',fsize);
    
    text(-2,1.25,texta,'color',fcolor,'fontsize',fsize,'rotation',0,'horizontalalignment','center','verticalalignment','middle');
    text(2,1.25,textb,'color',fcolor,'fontsize',fsize,'rotation',0,'horizontalalignment','center','verticalalignment','middle');
end