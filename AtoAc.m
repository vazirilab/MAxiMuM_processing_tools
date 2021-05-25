function [Ac_keep,acx,acy,acm] = AtoAc(A_keep,tau,d1,d2)

tau = tau(1);

x = 1:d2;
y = 1:d1;
[X,Y] = meshgrid(x,y);
Ac_keep = zeros(4*tau+1,4*tau+1,size(A_keep,2),'single');

acx = zeros(1,size(A_keep,2));
acy = acx;
acm = acx;

parfor ijk = 1:size(A_keep,2)
    
    AOI = reshape(single(full(A_keep(:,ijk))),d1,d2);
    cx = round(trapz(trapz(X.*AOI))./trapz(trapz(AOI)));
    cy = round(trapz(trapz(Y.*AOI))./trapz(trapz(AOI)));
    
    acx(ijk) = cx;
    acy(ijk) = cy;
    acm(ijk) = sum(AOI(:));
    
    sx = max([cx-2*tau 1]); % handle cases where neuron is closer than 3*tau pixels to edge of FOV
    sy = max([cy-2*tau 1]);
    ex = min([cx+2*tau d2]);
    ey = min([cy+2*tau d1]);
    
    AOIc = nan(4*tau+1,4*tau+1);
    AOIc(1:(ey-sy+1),1:(ex-sx+1)) = AOI(sy:ey,sx:ex);
    Ac_keep(:,:,ijk) = single(AOIc);
    
end