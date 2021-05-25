%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% comparePlanes.m
%
% J. Demas 09/22/2019
%
% This routine collates and thresholds all neurons found within the MAxiMuM
% data set.
%
% Inputs:
% path: directory with output files from planarSegmentation routine
% (char)
% z0: axial depth of the shallowest plane in the volume (1,1) (um)
% r_thr: threshold for spatial correlation (1,1) (0-1)
% pixel_resolution: size of each pixel in microns (1,1) (um)
% min_snr: minimum signal to noise ratio acceptable for found components
% (1,1)
% frameRate: recording frame/volume rate (1,1) (Hz)
% FOVx: size of the field of view in X (1,1) (um)
% FOVy: size of the field of view in Y (1,1) (um)
%
% Outputs:
% T_all: neuronal time series (Km,T) (single)
% C_all: denoised time series (Km,T) (single)
% nx: centroid in x direction for each neuron (1,Km) (single)
% ny: centroid in y direction for each neuron (1,Km) (single)
% nz: centroid in z direction for each neuron (1,Km) (single)

function comparePlanes(path,z0,r_thr,pixel_resolution,min_snr,frameRate,FOVx,FOVy)

tau = ceil(7.5/pixel_resolution);

merge_thr = 0.8;
Kms = zeros(30,1);
num_corr_no_ovp = 0; 
num_ovlpd = 0;

p = load([path 'caiman_output_plane_1.mat']); % Load data from first plane
rVals = p.rVals;

Tinit = p.T_keep; 
decay_time = 0.5;
Nsamples = ceil(decay_time*frameRate);
min_fitness = log(normcdf(-min_snr))*Nsamples;
[fitness] = compute_event_exceptionality(Tinit,Nsamples,0); % Calculate the fitness of all the traces for SNR thresholding
clear Tinit

if size(rVals)>0 % If neurons were found, unpack the relevant fields
    kp = logical(rVals>r_thr & fitness<min_fitness); % keep only neurons with specified r values and SNRs

    T = p.T_keep(kp,:);
    C = p.C_keep(kp,:);
    Y = p.Ym;
    A = p.Ac_keep(:,:,kp);
    K = size(T,1);
    Kms(1) = K;
    N = zeros(K,4);
    N(:,1) = p.acm(kp)';
    N(:,2) = p.acy(kp)';
    N(:,3) = p.acx(kp)';
    N(:,4) = 1;
    
else % If no neurons found, make dummy variables
    fff = p.f;
    bbb = p.b;
    T = NaN(1,size(fff,2));
    C = NaN(1,size(fff,2));
    Y = p.Ym;
    A = NaN(size(bbb,1),1);
    K = 1;
    Kms(1) = K;
    N = zeros(K,4);
end

T_all = T;
N_all = N;
C_all = C;

c = load([path 'three_neuron_mean_offsets.mat'],'offsets'); % Load offsets
offsets = round(c.offsets);

xo = cumsum(-offsets(:,2));
xo = xo-min(xo);

yo = cumsum(-offsets(:,1));
yo = yo-min(yo);

for ijk = 2:30 % Cycle through the rest of the planes
    
    disp(['Beginning calculation for plane ' num2str(ijk) ' of 30...'])
    
    pm = load([path 'caiman_output_plane_' num2str(ijk) '.mat']);
    
    Tinit = pm.T_keep;
    [fitness] = compute_event_exceptionality(Tinit,Nsamples,0);
    
    clear Tinit
    
    rValsm = pm.rVals;

    kpm = logical(rValsm>r_thr & fitness<min_fitness);
    
    Tm = pm.T_keep(kpm,:);    
    Cm = pm.C_keep(kpm,:);
    Ym = pm.Ym;
    Am = pm.Ac_keep(:,:,kpm);
    Km = size(Tm,1);
    Kms(ijk) = Km;
    Nm = zeros(Km,4);
    Nm(:,1) = pm.acm(kpm)';
    Nm(:,2) = pm.acy(kpm)';
    Nm(:,3) = pm.acx(kpm)';
    Nm(:,4) = ijk;

    Nm(:,2) = Nm(:,2) + cumsum(xo(ijk));
    Nm(:,3) = Nm(:,3) + cumsum(yo(ijk));
    
    if size(T,1)>0 && size(Tm,1)>0 % Find correlation between time traces for neurons in plane n and plane n+1
        RR = corr(T',Tm');
        MM = RR;
        MM(RR<merge_thr) = 0; % If correlation is below threshold, ignore
    else
        MM = 0;
    end
    
    if sum(sum(MM))>0 % If there are significant correlations, look for spatial overlap
        
        inds = find(MM(:));
        [ys,xs] = ind2sub(size(MM),inds);
        mm = MM(MM>0);
        [mm,sinds] = sort(mm,'ascend');
        ys = ys(sinds);
        xs = xs(sinds);
        
        Nk = numel(ys);
        
        for xyz = 1:Nk
            k = ys(xyz);
            km = xs(xyz);
            
            distance = sqrt(abs(N(k,2)-Nm(km,2)).^2 + abs(N(k,3)-Nm(km,3)).^2); % Determine distance between centroids of the components
            
            overlapped = distance<3*tau; % If centroids are within 3tau, define as overlapped

            if overlapped
                if ijk>2
                    indbuffer = sum(Kms(1:(ijk-2)));
                else indbuffer = 0;
                end
                
                T_all(indbuffer+k,:) = NaN(1,size(T_all,2));
                C_all(indbuffer+k,:) = NaN(1,size(T_all,2));
                N(indbuffer+k,:) = NaN(1,4);
                
                
                % Use weighted averaging to update trace and coordinates
                new_T = (T(k,:).*N(k,1) + Tm(km,:).*Nm(km,1))./(N(k,1) + Nm(km,1)); 
                new_C = (C(k,:).*N(k,1) + Cm(km,:).*Nm(km,1))./(N(k,1) + Nm(km,1));
                new_x = round((N(k,2)*N(k,1) + Nm(km,2)*Nm(km,1))./(N(k,1) + Nm(km,1)));
                new_y = round((N(k,3)*N(k,1) + Nm(km,3)*Nm(km,1))./(N(k,1) + Nm(km,1)));
                new_z = (N(k,4)*N(k,1) + Nm(km,4)*Nm(km,1))./(N(k,1) + Nm(km,1));
                new_sum = N(k,1) + Nm(km,1);
                
                Tm(km,:) = new_T;
                Cm(km,:) = new_C;
                Nm(km,:) = [new_sum new_x new_y new_z];
                
                num_ovlpd = num_ovlpd+1;
                
            else
                num_corr_no_ovp = num_corr_no_ovp+1;
            end
        end
    end

    % Append kept components to catalogue of all components in the data set
    T_all = cat(1,T_all,Tm);
    N_all = cat(1,N_all,Nm);
    C_all = cat(1,C_all,Cm);
    T = Tm;
    C = Cm;
    Y = Ym;
    A = Am;
    N = Nm;

    clear pm Tm Ym Cm RR MM XX Am Nm xxx XXX yyy YYY

end

is = ~isnan(sum(T_all,2));
T_all = T_all(is,:);
C_all = C_all(is,:);
N_all = N_all(is,:);

C_all = zeros(size(T_all),'single');

% Re-deconvolve raw traces to get updated component traces
disp('De-convolving raw traces...')
parfor j = 1:size(T_all,1);
    spkmin = 0.5*GetSn(T_all(j,:));
    [cc, spk, opts_oasis] = deconvolveCa(T_all(j,:),'ar2','optimize_b',true,'method','thresholded',...
        'optimize_pars',true,'maxIter',100,'smin',spkmin);    
    cb = opts_oasis.b;
    
    C_all(j,:) = full(cc(:)' + cb);
    
end

% If deconvolution fails, use raw trace
if sum(isnan(C_all(:)))>0 
    inds = find(isnan(sum(C_all,2)));
    disp(['Replacing ' num2str(numel(inds)) ' traces where de-convolution failed...'])
    
    for ijk = 1:numel(inds)
       C_all(inds(ijk),:) = T_all(inds(ijk),:); 
    end
end

%% Z plane correction

% Use Z positions from most recent pollen calibration to correct locations

try
    open([path 'pollen_calibration_Z_vs_N.fig'])
    fig = gcf;
    do = findobj(fig,'-property','Ydata');
    x = [do(3,1).XData do(2,1).XData];
    y = [do(3,1).YData do(2,1).YData];
    ftz = fit(x',y','cubicspline');
    close(fig)
catch % If calibration figure does not exist, use naive calibration
    x = 1:30;
    y = linspace(0,450,30);
    ftz = fit(x',y','cublicspline');
end

%% X, Y positions and Z field curvature correction

% Move coordinates from pixel space to physical space, correct for field
% curvature
ny = (FOVy./size(Y,1)).*N_all(:,2);
nx = (FOVx./size(Y,2)).*N_all(:,3);
nz = N_all(:,4);

nz = ftz(nz);
curvz = 158/2500^2;
nz = nz - curvz.*((ny-FOVy/2).^2 + (nx-FOVx/2).^2);
nz = nz+z0;

keep = logical(nz>0); % Ignore "neurons" found above the surface of the brain

T_all = T_all(keep,:);
C_all = C_all(keep,:);
N_all = N_all(keep,:);
nx = nx(keep);
ny = ny(keep);
nz = nz(keep);

%%
disp('Planes collated. Saving data...')
savefast([path 'collated_caiman_output_minSNR_' strrep(num2str(min_snr),'.','p') '.mat'],'T_all','nx','ny','nz','C_all')
disp('Routine complete.')

NN = std(T_all - movmean(T_all,Nsamples,2),[],2);
MM = max(movmean(T_all,Nsamples,2),[],2);
Z = MM./NN;

figure; histogram(Z,0:0.2:20)
xlabel('Z-score')
ylabel('Neurons')
saveas(gcf,[path 'all_neuron_Zscore.fig'])

figure
histogram(100.*max(T_all,[],2))
xlabel('Max \DeltaF/F_0 (%)')
ylabel('Neurons')
saveas(gcf,[path 'all_neuron_maxDF.fig'])