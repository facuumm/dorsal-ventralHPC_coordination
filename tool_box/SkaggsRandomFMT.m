% This function asign random positions to the spikes 100 times, each time 
% skaggs was calculated and then 90-quantile was defined.
%
% --- Inputs ---
% spks: vector, spikes times
% pos: matrix, positions in x and their respective time stamps (column1: time / column2: position in x axis)
% smooth: int for smooth the firing curve
% Xedges: int for spatial bins for the firing curve
%
% --- OUTPUT ---
% q: float, 50-quantile from random skaags values.
% Silva Azul 2023 - edited by Morici Juan Facundo

function q = SkaggsRandom(spks, pos, smooth, Xedges)

info = zeros(100,1);
for S = 1:100
    a = min(spks) + rand*10;% 1° tiempo
    isi = diff(spks);
    isi_per = isi(randperm(length(isi)));%shuffle isi
    
    % generation of spiking activity by keeping the isi
    spks_R = [a; isi_per];
    spks_R =cumsum(spks_R);
    
    [curve , stats] = FiringCurve(pos , spks_R , 'smooth' , smooth , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.2);
    info(S,1) = stats.specificity;
end
    q = quantile(info, 0.5);
end 