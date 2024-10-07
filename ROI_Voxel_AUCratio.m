clear; close all

if exist('img.mat')
    load img.mat
else
    readBrukerStudy
    uu = recoBrukerKspace(9);
    bb = recoBrukerKspace(10);
    bkg = sum(abs(bb),3);
    save img.mat
end
load img.mat
load ROImask.mat

Nt = 90; % Number of time points
Ms = 16; % Matrix size

% Background noise setup: Last ten slices of whole FOV

ss = Nt-9; % Backgound noise start slice
es = Nt; % Backgound noise end slice

B = zeros(1,2);

C1 = squeeze(uu(2:end-1,2:end-1,1,1,ss:es)); % lactate noise: middle frame (avoid outermost voxels) of last ten slices
C2 = squeeze(uu(2:end-1,2:end-1,1,2,ss:es)); % pyruvate noise: middle frame (avoid outermost voxels) of last ten slices

B(1) = std(reshape(C1,1,[]));

B(2) = std(reshape(C2,1,[]));

% Calculate voxel AUC

a = sum(mask1,1);mc = nnz(a); % nnz: number of nonzero matrix elements
b = sum(mask1,2);mr = nnz(b);
s = sum(mask1(:));
N=1;

OVLP = zeros(Nt,s); % one voxel lactate phase corrected
OVPP = zeros(Nt,s); % one voxel pyruvate phase corrected
fa = find(a);
fb = find(b);
fc = fa(1);
fr = fb(1);

for jj = 1:Ms
    for ll = 1:Ms
        if mask1(jj,ll)==1

            OVL = squeeze(uu(jj,ll,1,1,:)); % one voxel lactate
            OVP = squeeze(uu(jj,ll,1,2,:)); % one voxel pyruvate

            OVLP(:,N) = OVL.*exp(-1i*(angle(sum(OVL(:))))); % phase corrected OVL
            OVPP(:,N) = OVP.*exp(-1i*(angle(sum(OVP(:))))); % phase corrected OVP
         
            N = N+1;
        end
    end
end

% Plot voxel AUC
f(322) = figure(322);
LineWidth = 1;
rlw = 0.5; %reference line width
set(gcf,'Position', [2000,10,700,400])

mz = mask1(fr:fr+mr-1,fc:fc+mc-1); %mask1 zoomed
r = reshape(mz',mr*mc,1);
for ii=1:mr*mc
    if r(ii) == 1
    cs = sum(r(1:ii));
    subplot(mr,mc,ii)
    plot(real(OVPP(:,cs)),'g','LineWidth',LineWidth)
    hold on
    plot(real(OVLP(:,cs)),'b','LineWidth',LineWidth)

    plot(linspace(1,90,90), B(2)*ones(90,1),'g','LineWidth',LineWidth*rlw)
    plot(linspace(1,90,90), -B(2)*ones(90,1),'g','LineWidth',LineWidth*rlw)

    plot(linspace(1,90,90), B(:,1)*ones(90,1),'b','LineWidth',LineWidth*rlw)
    plot(linspace(1,90,90), -B(:,1)*ones(90,1),'b','LineWidth',LineWidth*rlw)
    hold off  
    else 
        continue
    end
end

savefig(f(322),'voxel AUC.fig');saveas(f(322),'voxel AUC.png');

%% SNR and AUC Ratio Calculation

stp = 10; etp = 80; % Start Time Point; End Time Point

NStd = B(2);% noise standard deviation

SNRtcV = zeros(1,s);% SNR total carbon in one voxel
SNRV = zeros(1,s);% SNR in one voxel: pyr/std.
AUCp = zeros(1,s);% AUC pyruvate in one voxel
AUCl = zeros(1,s);% AUC lactate in one voxel

for mm = 1:s

    AUCp(mm) = sum(real(OVPP(stp:etp,mm)));
    AUCl(mm) = sum(real(OVLP(stp:etp,mm)));

    SNRtcV(mm) = (AUCp(mm)+AUCl(mm))/(NStd*sqrt(2*(etp-stp+1)));
    SNRV(mm) = max(real(OVPP(stp:etp,mm)))/NStd;

end

AUCr = zeros(1,s);% AUC ratio in one voxel
AUCr = AUCl./AUCp;