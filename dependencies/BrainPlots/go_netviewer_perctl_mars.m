function go_netviewer_perctl_mars(C,thresh, varargin)

% varargin{1} : maximum sphere width if scaling outputs for comparison
% varargin{2} : limit if scaling outputs for comparison
% varargin{3} : labels for regions

C(eye(size(C))==1) = 0;
if thresh < 1; thresh = thresh*100; end

if nargin>2 & ~isnan(varargin{1})
    limit = varargin{1};
else
    limit = prctile(C(triu(ones(size(C)),1)==1),thresh);
end

mask = abs(C) >= limit;

cLims = [-max(abs(C(:))) max(abs(C(:)))];
C(~mask) = NaN;
edgeLims = sum(~isnan(C));
sphereCols = repmat([0 0 0]/255, 82, 1);
edgeLims = [1 10];
% Set the sphere widths to 1
% sphereWidths = repmat(1,78,1);
% Set them to only visible if a node survives thersholding

sphereWidths = sum(~isnan(C))*0.5;

if nargin>2 & ~isnan(varargin{1})
    sphereWidths = 5*sphereWidths./varargin{2};
else
    sphereWidths = 5*sphereWidths./max(sphereWidths(:));
end

load aalviewer
load sourcepos_mars_mni
mnipos = aalviewer.centroids;
mnipos_mars = round(sourcepos_mars.*1000);
% figure;clf
% Plot the cortical mesh
set(gcf,'color',[1 1 1]);
axis off
p = patch('faces',aalviewer.faces,'vertices',aalviewer.vertices,'edgecolor','none','facecolor','k','facealpha',0.05);
hold on


cmap      = colormap(RdBu);
if cLims(1) > 0
    cmap = cmap(129:end,:);
    isSingleColour = true;
elseif cLims(2) < 0
    cmap = cmap(1:128,:);
    isSingleColour = true;
else
    
    isSingleColour = false;
end

emap = (linspace(-3,3,length(cmap))).^2;

[i,j] = find(~isnan(C));
for p=length(i):-1:1,
    colorInd(p) = closest(C(i(p),j(p)),linspace(cLims(1),cLims(2),size(cmap,1)));
end

for p=1:length(i),
    edgecolour  = cmap(colorInd(p),:);
    edgeWeight = emap(colorInd(p));
    line(mnipos_mars([i(p) j(p)],1),mnipos_mars([i(p) j(p)],2),mnipos_mars([i(p) j(p)],3),'color',edgecolour,'linewidth',edgeWeight)
end


set(gca,'clim',cLims);
colormap(RdBu);
hc = colorbar;
FONTSIZE = 14;
if isSingleColour,
    YTicks = [cLims(1) cLims(2)];
else
    YTicks = [cLims(1) 0 cLims(2)];
end%if


% YTL = get(hc,'yticklabel');
% set(hc,'yticklabel',[repmat(' ',size(YTL,1),1), YTL]);



% Plot Spheres

for ii = 1:length(mnipos_mars);
    [x y z] = sphere(10);
    sw = sphereWidths(ii);
    surf(sw*x+mnipos_mars(ii,1),sw*y+mnipos_mars(ii,2),sw*z+mnipos_mars(ii,3),'edgecolor','none','facecolor','b','facealpha',1);
end

if nargin>4
    text(mnipos_mars(sphereWidths>0,1), mnipos_mars(sphereWidths>0,2), mnipos_mars(sphereWidths>0,3), strrep(varargin{3}(sphereWidths>0), '_', ' '))
end

axis vis3d
axis equal
axis off
hold off
rotate3d('on')
set(gcf,'renderer','opengl')
colorbar off

end

function i = closest(a,k)
%CLOSEST finds index of vector a closest to k
assert(isscalar(k) | isscalar(a));

[~,i] = min(abs(a-k));
end


function [cmap] = RdBu(varargin)

switch nargin
    case 0
        ncols = 256;
        pn = 'div';
        deep = 0;
    case 1
        ncols = varargin{1};
        pn = 'div';
        deep = 0;
    otherwise
        ncols = varargin{1};
        if sum(strcmp('type',varargin));
            pn = varargin{find(strcmp('type',varargin))+1};
        else
            pn = 'div';
        end
        if sum(strcmp('deep',varargin));
            deep = 1;
        end
end


if rem(ncols,2)~=0;
    error('Can only accept even numbers');
end

ncols = ncols./2;


% colours
lo        = [5 48 97] / 255;
bottom    = [5 113 176] / 255;
botmiddle = [146 197 222] / 255;
middle    = [247 247 247] / 255;
topmiddle = [244 165 130] / 255;
top       = [202   0  32] / 255;
hi        = [103 0 31] / 255;

% Find ratio of negative to positive
if strncmp(pn,'div',3) || strncmp(pn,'neg',3)
    
    
    % Just negative
    if deep
        neg = [lo; bottom; botmiddle; middle];
    else
        neg = [bottom; botmiddle; middle];
    end
    len = length(neg);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, ncols);
    neg128 = zeros(ncols, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        neg128(:,i) = min(max(interp1(oldsteps, neg(:,i), newsteps)', 0), 1);
    end
    
    cmap = neg128;
    
end

if strncmp(pn,'div',3) || strncmp(pn,'pos',3)
    % Just positive
    if deep
        pos = [middle; topmiddle; top; hi];
    else
        pos = [middle; topmiddle; top];
    end
    len = length(pos);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, ncols);
    pos128 = zeros(ncols, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        pos128(:,i) = min(max(interp1(oldsteps, pos(:,i), newsteps)', 0), 1);
    end
    cmap = pos128;
end

if strmatch(pn,'div')
    % And put 'em together
    cmap = [neg128; pos128];
end
end
