close all
clear



% ls = getLebedevSphere(350);
ls = getLebedevSphere(4334);
[dirs(:,1), dirs(:,2), ~] = cart2sph(ls.x,ls.y,ls.z);
dirs = rad2deg(dirs);


figure('Name','Mollweide test','NumberTitle','off','WindowStyle','docked');

% data
load moonalb % a 540x1080 matrix of values is loaded along with moonalbrefvec=[3;90;0]
moonalbrefvec(1)=3; % this is the subdivisions per degree - thus 180*3=540 & 3*360=1080
moonalbrefvec(2)=90; % NW lat value
moonalbrefvec(3)=180; % NW long value

% setup the plot
mymap=colormap('jet');
mymap(1,1:3)=1;
colormap(mymap);
caxis([0 500]); 
axesm('MapProjection','mollweid','MapLatLimit',[-90 90],'Gcolor','black','GLineWidth',2.0,'MLineLocation',[-135 -90 -45 0 45 90 135],'PLineLocation',30); 
axis off;
grid on;
gridm('on');


plabel('LabelFormat','none');
plabel('on')
mlabel('equator');
mlabel('FontColor','white');
mlabel('on');

% geoshow(moonalb, moonalbrefvec, 'DisplayType', 'texturemap');


scatterm(dirs(:,2),dirs(:,1))


% color bar
c=colorbar('location','southoutside','box','on','color',[0 0 0]);
c.Label.String='T_{sky} (K)';
c.Limits=[0 500];
c.Ticks=0:50:500;
c.FontSize=12;

% labels
fsize = 10;
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



%  figure(1);clf(1);
%  dec=-90:1:90;
%  ra=180:-1:-180;
%  mydata=rand(181,361);
%  m_proj('mollweide','clongitude',0);
%  obj=m_pcolor(ra,dec,mydata); 
%  shading interp; 
%  m_grid('xaxislocation','middle','xtick',8,'ytick',7,'fontsize',7,'color','white');
% 
%  mymap=colormap('jet');
%  colormap(mymap);
%  h=colorbar('h');
%  caxis([0 1]);
% 
% 

