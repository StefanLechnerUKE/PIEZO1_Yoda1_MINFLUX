%%  This script loads raw localistation data and calculates probability densities shown in Fig. 4D, E and F in Verkest et al. (doi???) 
%   Lechner SG, 10.7.2025, Version 1.0
%   Dependencies: Fig4_examples_FINAL.mat must be available in the same folder as script

%% load data
load('Fig4_examples_FINAL.mat');
fns = fieldnames(Fig4_examples);

for k = 1:size(fns,1)
    %% extract protomer coordinates    
        TrimA = Fig4_examples.(fns{k}).A;
        TrimB = Fig4_examples.(fns{k}).B;
        TrimC = Fig4_examples.(fns{k}).C;
    
    %% create 2D in-plane projection
        [TrimAtrans, TrimBtrans, TrimCtrans, AllCenters] = CreateInPlaneProjection(TrimA, TrimB, TrimC, k);

    %% fit gaussians to protomers in 2D in-plane projected data
        [X, Y, AllZ] = Fit2DGaus_2_protomer(TrimAtrans, TrimBtrans, TrimCtrans);

    %% create or add to plot
    if k==1
        figure('Name','Control', 'Position',[100 100 1400 500]);
        set(gcf,'renderer','Painters');
        tiledlayout(2,5);
    elseif k==11
        figure('Name','P1_A2094W', 'Position',[100 100 1400 500]);
        set(gcf,'renderer','Painters');
        tiledlayout(2,5);
    elseif k==21
        figure('Name','P1_V1714A_F1715A', 'Position',[100 100 1400 500]);
        set(gcf,'renderer','Painters');
        tiledlayout(2,5);
    end
    nexttile
        hold on;
        s=pcolor(X,Y,AllZ);
        s.EdgeColor = 'none';
        colormap hot; axis tight; axis equal;
        xlabel('(nm)'); ylabel('(nm)');
        ax = gca; ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14;
        ax.Box = 'on'; ax.BoxStyle = 'full';
        colorbar; c = colorbar; c.Label.String = 'probability'; c.Label.FontSize = 14; c.FontSize = 14;  
end

%% clean up
    clearvars -except Fig4_examples;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, Y, AllZ] = Fit2DGaus_2_protomer(TrimAtrans, TrimBtrans, TrimCtrans)
        x = linspace(-30, +30, 61);
        y = linspace(-30, +30, 61);
        [X, Y] = meshgrid(x, y);
        XY = [X(:) Y(:)];
        TrimAtrans2D = TrimAtrans(:,1:2);
        TrimBtrans2D = TrimBtrans(:,1:2);
        TrimCtrans2D = TrimCtrans(:,1:2);
        GausA = fitgmdist(TrimAtrans2D,1);
        Za = reshape(pdf(GausA, XY), size(X));
        GausB = fitgmdist(TrimBtrans2D,1);
        Zb = reshape(pdf(GausB, XY), size(X));
        GausC = fitgmdist(TrimCtrans2D,1);
        Zc = reshape(pdf(GausC, XY), size(X));
        AllZ = Za+Zb+Zc;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TrimAtrans, TrimBtrans, TrimCtrans, AllCenters] = CreateInPlaneProjection(TrimA, TrimB, TrimC,k);
    
    AddRotation = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -140 140 145 0 0 0 0 0 0 0 0 0 0 0];
    AddROT = AddRotation(k);

    % shift coordinates to center of first protamer
    CenterA = [mean(TrimA(:,1)) mean(TrimA(:,2)) mean(TrimA(:,3))];
    TrimAtrans = TrimA-CenterA;
    TrimBtrans = TrimB-CenterA;
    TrimCtrans = TrimC-CenterA;
    MyTrimer = [TrimAtrans; TrimBtrans;TrimCtrans];

    % calculate center of each protamer
    ProtA = [mean(TrimAtrans(:,1)) mean(TrimAtrans(:,2)) mean(TrimAtrans(:,3))];
    ProtB = [mean(TrimBtrans(:,1)) mean(TrimBtrans(:,2)) mean(TrimBtrans(:,3))];
    ProtC = [mean(TrimCtrans(:,1)) mean(TrimCtrans(:,2)) mean(TrimCtrans(:,3))];
    TrimerOrig = [ProtA; ProtB; ProtC];
    
    %% rotate data to create projection to plane formed by ABC
        % rotate around Z
        a = ProtB(:,1);
        b = ProtB(:,2);
        c = sqrt(a^2+b^2);
        angle = asind(b/c);
        if a>0 && b>0
            angle = angle;
        elseif a>0 && b<0
            angle = 360+angle;
        elseif a<0 && b>0
            angle = 180-abs(angle);
        elseif a<0 && b<0
            angle = 180 + abs(angle);
        end
        
        RotMatrixZ = rotz(angle);
        ProtA = ProtA*RotMatrixZ;
        ProtB = ProtB*RotMatrixZ;
        ProtC = ProtC*RotMatrixZ;
        
    % rotate around y
        a = ProtB(:,1);
        b= ProtB(:,3);
        c = sqrt(a^2+b^2);
        angle = -asind(b/c);
        if a>0 && b>0
            angle = angle;
        elseif a>0 && b<0
            angle = 360+angle;
        elseif a<0 && b>0
            angle = 180-abs(angle);
        elseif a<0 && b<0
            angle = 180 + abs(angle);
        end
        RotMatrixY = roty(angle);
        ProtA = ProtA*RotMatrixY;
        ProtB = ProtB*RotMatrixY;
        ProtC = ProtC*RotMatrixY;

    % rotate around X
        a = ProtC(:,2);
        b = ProtC(:,3);
        c = sqrt(a^2+b^2);
        angle = asind(b/c);
        if a>0 && b>0
            angle = angle;
        elseif a>0 && b<0
            angle = 360-abs(angle);
        elseif a<0 && b>0
            angle = 180 - abs(angle);
        elseif a<0 && b<0
            angle = 180 + abs(angle);
        end
        RotMatrixX = rotx(angle);
        ProtA = ProtA*RotMatrixX;
        ProtB = ProtB*RotMatrixX;
        ProtC = ProtC*RotMatrixX;
    
    % create array with rotated protamer means
        TrimerRot = [ProtA; ProtB; ProtC];
    
    % apply rotation matrices to individual localisations
        TrimAtrans = TrimAtrans*RotMatrixZ;
        TrimAtrans = TrimAtrans*RotMatrixY;
        TrimAtrans = TrimAtrans*RotMatrixX;
        
        TrimBtrans = TrimBtrans*RotMatrixZ;
        TrimBtrans = TrimBtrans*RotMatrixY;
        TrimBtrans = TrimBtrans*RotMatrixX;
        
        TrimCtrans = TrimCtrans*RotMatrixZ;
        TrimCtrans = TrimCtrans*RotMatrixY;
        TrimCtrans = TrimCtrans*RotMatrixX;
    
    % calculate center positions of final rotated trimer (projection)
        CenterTrimTransA = [mean(TrimAtrans(:,1)) mean(TrimAtrans(:,2)) mean(TrimAtrans(:,3))];
        CenterTrimTransB = [mean(TrimBtrans(:,1)) mean(TrimBtrans(:,2)) mean(TrimBtrans(:,3))];
        CenterTrimTransC = [mean(TrimCtrans(:,1)) mean(TrimCtrans(:,2)) mean(TrimCtrans(:,3))];
        AllCenters =[CenterTrimTransA; CenterTrimTransB; CenterTrimTransC];
   
    % additional rotation for figure
        RotMatrixZ = rotz(AddROT);
        TrimAtrans = TrimAtrans*RotMatrixZ;
        TrimBtrans = TrimBtrans*RotMatrixZ;
        TrimCtrans = TrimCtrans*RotMatrixZ;
        AllCenters = AllCenters*RotMatrixZ;

    % shift final coordinates such that the trimer center is in the origin
        FinalCenter = mean(AllCenters,1);
        TrimAtrans = TrimAtrans-FinalCenter;
        TrimBtrans = TrimBtrans-FinalCenter;
        TrimCtrans = TrimCtrans-FinalCenter;
        AllCenters = AllCenters-FinalCenter;


end
