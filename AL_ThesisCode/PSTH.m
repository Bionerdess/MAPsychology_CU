function [out,f] = PSTH(mat,varargin)

%%% PSTH takes an M x N matrix and calculates a peri-stimulus time
%%% histogram.
%%% To compare multiple types of stimuli, X can be a cell of MxN matrixes
%%% X := M observations x N time bins
%%% [X] = PSTH(mat) returns a matrix X of size 2 x N.
%%% X(1,:) := mean PSTH
%%% X(2,:) := SEM PSTH
%%% [X,F] = PSTH(X) returns a figure plotting the mean PSTH +/- SEM
%%% [...] = PSTH(X,'vertline') plots vertical lines at values indicated in
%%% vertline
%%% [...] = PSTH(X,...,'color') plots each stimulus in a specified color
%%% (stimuli x 3 -> stimuli x R G B)
%%% (matrix dimension between X (class = 'cell') and color have to agree
%%%
%%% Anna-Lena Schlenner, July 2021
%%% al.schlenner@gmail.com



if isa(mat,'cell')
    out = cell(length(mat),1);
    for i=1:length(mat)
        n = size(mat{i},1);
        out{i}(1,:) = mean(mat{i},1);
        out{i}(2,:) = std(mat{i},1)/sqrt(n);
    end
elseif isa(mat,'double')
    n = size(mat,1);
    out(1,:) = mean(mat,1);
    out(2,:) = std(mat,1)/sqrt(n);
end
if nargout >1 
    if ~isempty(varargin)
        [vLine,Color] = getargs(mat,varargin);
    else
        vLine = [];
        Color = [];
    end
    [f] = plot_PSTH(out,vLine,Color);
end


%----------------------------------------------
function [vLine,color] = getargs(mat,args)
p = inputParser;
for a=1:2:size(args,2)
    addParameter(p,args{a},[]);
end

m=length(mat);
parse(p,args{:});
vLine = [];
color=[];
for b=1:length(p.Parameters)
    if strcmp(p.Parameters{b},'vertline')
        vLine = p.Results.(p.Parameters{b});
    end
    if strcmp(p.Parameters{b},'color')
        color = p.Results.(p.Parameters{b});
        a = size(color,1);
        
        if m~=a
            error(message('signal:PSTH:MatrixDimensionMismatch', num2str(m),num2str(a)));
        end
        if size(color,2)~=3
            error(message('signal:PSTH:ColorHasToProvide3Values',num2str(size(color,2))));
        end
    end
end



%-------------------------------------------------
function [fig] = plot_PSTH(dat,vLine,Color)

fig = figure;
plot([0 size(dat{1},2)],[0 0],'k-')
hold on
for j=1:length(dat)
    x = 1:size(dat{j},2);
    figure(fig)
    if ~isempty(Color)
        fill([x fliplr(x)], [dat{j}(1,:)-dat{j}(2,:) fliplr(dat{j}(1,:)+dat{j}(2,:))],Color(j,:)*1.5)
        hold on
        plot(x,dat{j}(1,:),'Color',Color(j,:),'LineWidth',2)
        hold on
    else
        fill([x fliplr(x)], [dat{j}(1,:)-dat{j}(2,:) fliplr(dat{j}(1,:)+dat{j}(2,:))],[173 172 170]/255)
        hold on
        plot(x,dat{j}(1,:),'LineWidth',2)
        hold on
    end
    
end

if ~isempty(vLine)
    figure(fig)
    for w=1:length(vLine)
        plot([vLine(w) vLine(w)],[min(dat{j}(1,:))-1 max(dat{j}(1,:))+1],'k-','LineWidth',1.5)
        hold on
    end
end
