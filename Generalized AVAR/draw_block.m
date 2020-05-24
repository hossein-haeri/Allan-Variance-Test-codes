function draw_block(center,w,h,d)

e = 0.05;
w = w-e;
h = h-e;
d = d-e;


verts = [-w -h -d; ...
          w -h -d; ...
         -w  h -d; ...
          w  h -d; ...
         -w -h  d; ...
          w -h  d; ...
         -w  h  d; ...
          w  h  d]./2;
verts = verts + center;
faces = [3 4 8 7; ...
         4 2 6 8; ...
         2 1 5 6; ...
         1 3 7 5; ...
         7 8 6 5; ...
         1 2 4 3];
g = hgtransform;     
patch('Vertices',verts,'Faces',faces,'FaceColor',[.9 .35 .35],'Parent',g)
alpha(0.5)
view(3)
box on
axis vis3d
% daspect([1 1 1])
x = 0;
y = 0;
z = 0;


g.Matrix = makehgtform('xrotate',x,'yrotate',y,'zrotate',z);
hold on


end

