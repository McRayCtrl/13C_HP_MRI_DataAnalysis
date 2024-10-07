%% AUC pyruvate/lactate overlay plotting
clear
load ROImask.mat
close all

Nt = 90; % Number of time points
Ms = 16; % Matrix size
NHPA = 2; % Number of HP agents. 2 if only focusing on Lactate and pyruvate 

t = linspace(1,Nt,Nt);
FOV_HP = zeros(Ms,Ms,NHPA,Nt); % HP map of FOV, 3rd dim 1 for lactate 2 for pyruvate

for jj = 1:Ms
    for ll = 1:Ms

    FOV_lac = squeeze(uu(jj,ll,1,1,:));
    FOV_pyr = squeeze(uu(jj,ll,1,2,:));

    FOV_AUC_l = trapz(t,FOV_lac); %AUC of lactate
    FOV_AUC_p = trapz(t,FOV_pyr); %AUC of pyruvate

    FOV_pl = angle(FOV_AUC_l); %phase of AUC lactate
    FOV_pp = angle(FOV_AUC_p); %phase of AUC pyruvate

    FOV_HP(jj,ll,1,:) = real(FOV_lac*exp(-1i*FOV_pl)); %phase corrected HP lactate
    FOV_HP(jj,ll,2,:) = real(FOV_pyr*exp(-1i*FOV_pp)); %phase corrected HP pyruvate
    
    end
end

lac_dyn = squeeze(FOV_HP(:,:,1,:));
pyr_dyn = squeeze(FOV_HP(:,:,2,:));

pyr_AUC = sum(pyr_dyn,3);
lac_AUC = sum(lac_dyn,3);

pyr_AUC_mask = zeros(Ms,Ms,Nt);
lac_AUC_mask = zeros(Ms,Ms,Nt);

for jj = 1:Nt
    pyr_AUC_mask(:,:,jj) = pyr_AUC(:,:);
    lac_AUC_mask(:,:,jj) = lac_AUC(:,:);
end

%% FOV pyruvate AUC overlay
viewover(pyr_AUC_mask,bkg);

min_valp = min(pyr_AUC(:));
max_valp = max(pyr_AUC(:));
% colorbar
clim([min_valp, max_valp]);

%% FOV lactate AUC overlay
viewover(lac_AUC_mask,bkg);

min_vall = min(lac_AUC(:));
max_vall = max(lac_AUC(:));
% colorbar
clim([min_vall, max_vall]);

%% ROI pyruvate AUC overlay

ptr = squeeze(pyr_AUC_mask(:,:,1)).*mask1; % 13C pyruvate map on tumor region:16x16
ptr3 = repmat(ptr,[1,1,Nt]); % 13C pyruvate map on tumor region: 16x16x90
viewover(ptr3,bkg)
clim([0 max(ptr3(:))])
% make a figure that has a colorbar following ROI pyruvate AUC overlay
% color limit

figure;imagesc(randn(Ms,Ms))
colormap('hot');
clim([0 max(ptr3(:))])

c = colorbar;
c.LineWidth = 1;
c.FontSize = 15;
c.FontWeight = 'bold';

%% ROI lactate AUC overlay
ltr = squeeze(lac_AUC_mask(:,:,1)).*mask1; % 13C lactate map on tumor region:16x16
ltr3 = repmat(ltr,[1,1,Nt]); % 13C lactate map on tumor region: 16x16x90
viewover(ltr3,bkg)
clim([0 max(ltr3(:))])
% make a figure that has a colorbar following ROI pyruvate AUC overlay
% color limit

figure;imagesc(randn(Ms,Ms))
colormap('hot');
clim([0 max(ltr3(:))])

c = colorbar;
c.LineWidth = 1;
c.FontSize = 15;
c.FontWeight = 'bold';

%% lactate AUC overlay has pyruvate AUC overlay limit

new_lac_AUC_mask = lac_AUC_mask/max(pyr_AUC(:));
new_pyr_AUC_mask = pyr_AUC_mask/max(pyr_AUC(:));

viewover(new_lac_AUC_mask,bkg);
clim([0 1])

viewover(new_pyr_AUC_mask,bkg);
clim([0 1])

% heatmap colorbar plotting
figure;imagesc(randn(Ms,Ms))
colormap('hot');
clim([0 max(pyr_AUC_mask(:))])
% clim([0 max(lac_AUC_mask(:))])
% clim([0 1])
c = colorbar;
c.LineWidth = 1;
c.FontSize = 15;
c.FontWeight = 'bold';
