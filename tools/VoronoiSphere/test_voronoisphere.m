% Sctipt to test voronoisphere
% Random data
caseid = 0;
if caseid == 2
    xyz = [0  0;
           0  0;
           1 -1];
elseif caseid == 3
    theta = (0:2)*(2*pi/3);
    xyz = bsxfun(@plus, [0; 0; 1], 0.1*[cos(theta); sin(theta); zeros(size(theta))]);
else
    n = 100;
    xyz = randn(3,n);
end

xyz = bsxfun(@rdivide, xyz, sqrt(sum(xyz.^2,1)));
n = size(xyz,2);

[P, K, voronoiboundary, s] = voronoisphere(xyz);
sum(s)

%% Graphic
f = figure(1);
clf(f);
set(f,'Renderer','zbuffer');
ax = axes('Parent', f);
hold(ax, 'on');
axis(ax,'equal');

plot3(ax, xyz(1,:),xyz(2,:),xyz(3,:),'ro');
clmap = cool();
ncl = size(clmap,1);
if n < 20
    [X,Y,Z] = sphere();
    cl = clmap(1,:);
    surf(X,Y,Z,'FaceColor',cl,'Parent',ax,'EdgeColor','none');
end
for k = 1:n
    X = voronoiboundary{k};
    cl = clmap(mod(k,ncl)+1,:);
    fill3(X(1,:),X(2,:),X(3,:),cl,'Parent',ax,'EdgeColor','w');
end
axis(ax,'equal');
axis(ax,[-1 1 -1 1 -1 1]);

f = figure(2);
clf(f);
set(f,'Renderer','zbuffer');
ax = axes('Parent', f);
hold(ax, 'on');
%input('Type <CR>: ', 's');
cla(ax);
plot3(ax, P(1,:),P(2,:),P(3,:), 'ro');
for k = 1:n
    V = P(:,K{k});
    if ~isempty(V)
        V = V(:,[1:end 1]);
        for i=1:length(V)-1
            plot3(ax, V(1,i:i+1),V(2,i:i+1),V(3,i:i+1), 'r', 'Linewidth', 1);
        end
    end
end
axis(ax,'equal');
axis(ax,[-1 1 -1 1 -1 1]);
