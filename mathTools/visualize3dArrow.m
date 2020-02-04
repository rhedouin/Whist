function out = visualize3dArrow(directions,varargin)

if nargin == 3
    createFig = varargin{1};
    color = varargin{2};
elseif nargin == 2 
    createFig = varargin{1};
    color = 's';
elseif nargin == 1;
    createFig = 1;
    color = 's';
end

if createFig == 1
    figure;
end

axis([-2 2 -2 2 -2 2])
CameraPosition = [5 0 0]
view([0 90])

hold on
arrow3([0 0 0], [0 0 1],'--k',0.4,[],0)
arrow3([0 0 0], [0 1 0],'--k',0.4,[],0)     % Y
arrow3([0 0 0], [1 0 0],'--k',0.4,[],0)     %X

for j = 1:size(directions, 1)
    current_dir = directions(j,:);
    if j == 1
        arrow3([0 0 0],current_dir, ['-' color '5'],1,[],0)
    else
        arrow3([0 0 0],current_dir, ['-' color '5'],1,[],0)
    end    
end


hold off, axis off, camlight left
hold off, camlight left
set(gca,'CameraViewAngle',4)
text(1.4,0,0,'X','FontSize',14), text(0,1.4,0,'Y','FontSize',14)
text(0,0,1.4,'Z','VerticalAlignment','bottom',...
    'HorizontalAlignment','center','FontSize',14)

out = 0;
end




