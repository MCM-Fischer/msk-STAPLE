function [ U, Uridge ,LowestPoints_End ] = LSSLFitPatellaRidge( TR,U,nbSlice,StartDist,EndDist, debug_plot)
%Least Square Straight Line Fit on Patellar Ridge
%To obtain the direction of the ridge a Least Square Straight Line is
%fitted on the lowest points on slices of normal U, and the normal U is
%updated until convergence (of U or number of iterations > 100)


%% inputs Tests
if nargin<3
    StartDist = 0.025*range(TR.Points*U); %Offset distance from first points
    EndDist = StartDist; %Offset distance from last points
    nbSlice = 50;
end

if nargin<4
    StartDist = 0.025*range(TR.Points*U); %Offset distance from first points
    EndDist = StartDist; %Offset distance from last points
end

if nargin<6
    debug_plot=0;
end
%% Code
U0 = U;
reslts = [0;0;0];
k=0;
test = true;
while test && k < 100 % Cut at 100 iterations (convergence generally occuring before 10)
    k=k+1;
    Alt = linspace( min(TR.Points*U)+StartDist , max(TR.Points*U)-EndDist, nbSlice);
    LowestPoints = zeros(length(Alt),3);
    i=0;
    for d = -Alt
        i=i+1;
        [ Curves , ~ , ~ ] = TriPlanIntersect( TR, U , d );
        EdgePts = vertcat(Curves(:).Pts);
        [~,lowestPointID] = min(EdgePts(:,3));
        LowestPoints(i,:) = EdgePts(lowestPointID(1),:);
    end
    
    [Vridge_all,~] = eig(cov(LowestPoints));
    U = sign(U0'*Vridge_all(:,3))*Vridge_all(:,3);
    
    U(3) = 0; U = normalizeV( U );
    reslts(1:3,end+1) = U;
    test = reslts(1:3,end-1)'*reslts(1:3,end) > (1 - 100*eps);
end


if nargout == 2
Uridge = sign(U0'*Vridge_all(:,3))*Vridge_all(:,3);

elseif nargout == 3
        Uridge = sign(U0'*Vridge_all(:,3))*Vridge_all(:,3);
        LowestPoints_End = LowestPoints;
end

%% Plots
if debug_plot == 1
    figure()
    trisurf(TR,'facecolor','red','faceAlpha',0.5)
    hold on
    grid off
    axis equal
    pl3t(LowestPoints,'gs')
    plotArrow( U, 1.5, mean(LowestPoints), 30, 1, 'k')
    Uridge = sign(U0'*Vridge_all(:,3))*Vridge_all(:,3);
    plotArrow( Uridge, 1.5, mean(LowestPoints), 50, 1, 'c')
    hold off
end
end

