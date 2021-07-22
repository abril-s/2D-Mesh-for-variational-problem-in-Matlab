% ------------------------------------------------
% -------------------MALLA------------------------
% ------------------------------------------------

function[dt, rec] = malla(R,m,l,d)
% Inputs -----------------------------------------

% R = 1; R = radio del circulo
% m = 50;  num divisiones en y
% l = 5;  longitud del lado en x
% d = 4;   longitud del lado en y 

% Outputs ----------------------------------------

% dt - the triangulation produced
% rect - points without boundaries

% ------------------------------------------------
% -------------- main variables ------------------
% ------------------------------------------------

dy = d/m; % step in y
dx = dy; % step in x
n = round(l/(dx)); % numero de divisiones en x del rectangulo
nc = 2*ceil((2*pi)/(dy)); % num divisiones en el circulo

% ------------------------------------------------
% --------------- inner circle -------------------
% ------------------------------------------------

alpha = linspace(0,2*pi,nc)'; %  nc valores del angulo
circle = [R*cos(alpha),R*sin(alpha)] ; % coordenadas del circulo (0,1), (1,0)
numc = (1:nc)'; % numeramos los puntos del circulo
numc = [numc,numc+1]; % vertices
numc(end,2) = 1;
ss = size(numc,1); % numero de puntos en el circulo
% plot(circle(:,1),circle(:,2),'r*')
% axis equal
% hold on

% ------------------------------------------------
% ------------- inner rectangles -----------------
% ------------------------------------------------

dim = (n+2)*(m+2); % total number of points of all rectangles
cry = (linspace(-d/2,d/2, m+2)); % coordenada en y
crx = (linspace(-l/2,l/2, n+2)); % coordenada en x
rec0 = zeros(dim,2); % arreglo inicial coord de todos los rectangulos
% tienen dif numero de puntos pero el mismo step
% generacion de coordenadas
for j = 0:n+1
    for i=1:m+2
        rec0(j*(m+2)+i,2) = cry(i);
        rec0(j*(m+2)+i,1) = crx(j+1);
    end
end
%plot(rec0(:,1),rec0(:,2),'g*');
%axis equal
  
% ------------------------------------------------
% ---------- boundary identification -------------
% ------------------------------------------------

bx1 = zeros(n+2,2);
bx1(:,1) = crx;
bx1(:,2) = cry(1);
bx2 = zeros(n+2,2);
bx2(:,1)= crx;
bx2(:,2) = cry(m+2);
boundaryx = [bx1; bx2];
%plot(boundaryx(:,1),boundaryx(:,2),'g*');
by1 = zeros(m+2,2);
by1(:,1)= crx(1);
by1(:,2) = cry;
by2 = zeros(m+2,2);
by2(:,1)= crx(n+2);
by2(:,2) = cry;
boundaryy = [by1; by2];
%plot(boundaryy(:,1),boundaryy(:,2),'g*');

boundary = [boundaryx;boundaryy];

% ------------------------------------------------
% ---------- rectangles w/o boundary -------------
% ------------------------------------------------
rec = []; % rect inside w/o boundary
col1notin = ~ismember(rec0(:, 2), boundaryx);
col2notin = ~ismember(rec0(:, 1), boundaryy);
rec = rec0(col1notin & col2notin, :);

% ------------------------------------------------
% --------------- removing points ----------------
% ------------------------------------------------

auxc = polyshape(R*cos(alpha),R*sin(alpha)); % circulo generado para excluir puntos

q = isinterior(auxc,rec0(:,1),rec0(:,2)); % verificacion
rec0(q,:) = []; % pasamos solo la fila donde q diferente a 0
% here rec0 has all points between circle and outer rect including outer
% rect boundary
tot = [circle;rec0]; % total points
% plot(rec0(:,1),rec0(:,2),'g*');
% axis equal
% hold on
tri = delaunayTriangulation(tot);
% triplot(tri)
% axis equal

% ------------------------------------------------
% ------------ triangulacion final ---------------
% ------------------------------------------------

% remove triangles where all 3 vertices are on the circular boundary

TriVert = tri.ConnectivityList; % new array of conn.list
TriVert(all(TriVert <= ss,2),:) = []; % erase vertices (2 indicates rows, 1 columns)
dt = triangulation(TriVert,tri.Points);
triplot(dt);
axis equal;
% s9 = size(tri.Points,1);
% vxlabels = arrayfun(@(n) {sprintf('P%d', n)}, (1:s9)');
% Hpl = text(tri.Points(:,1),tri.Points(:,2),vxlabels,'FontWeight','bold','HorizontalAlignment',...
%    'center','BackgroundColor','none');
% ic = incenter(dt);
% numtri = size(dt,1);
% trilabels = arrayfun(@(x) {sprintf('T%d',x)}, (1:numtri)');
% Htl = text(ic(:,1),ic(:,2),trilabels,'FontWeight','bold', ...
%    'HorizontalAlignment','center','Color','blue');
end

