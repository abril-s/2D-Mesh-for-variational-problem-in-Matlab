% ejemplo de malla dando puntos iniciales
x = rand(10,1);
y = rand(10,1);
dt = delaunayTriangulation(x,y);
% plot(x,y,'o:');  plot(x,y,'o:') con lineas
triplot(dt,x,y,'r');
dt.Points;
triplot(dt, dt.Points(:,1), dt.Points(:,2), 'r');
dt.ConnectivityList;
triplot(dt.ConnectivityList, dt.Points(:,1), dt.Points(:,2), 'r')
hold on;
vxlabels = arrayfun(@(n) {sprintf('P%d', n)}, (1:10)');
Hpl = text(x,y,vxlabels,'FontWeight','bold','HorizontalAlignment',...
   'center','BackgroundColor','none');
ic = incenter(dt);
numtri = size(dt,1);
trilabels = arrayfun(@(x) {sprintf('T%d',x)}, (1:numtri)');
Htl = text(ic(:,1),ic(:,2),trilabels,'FontWeight','bold', ...
   'HorizontalAlignment','center','Color','blue');
hold off
dt.ConnectivityList(1,:) % puntos triangulo 1