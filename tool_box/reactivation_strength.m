function [R] = reactivation_strength(patterns , cond , SpikeTrain , Is , th , type , config)
% This function calculate Reactivation Strength (van de Ven et al (2016)).
% Find the peaks of the zscored assemblies strength and calculate the mean
% during Pre sleep and Post sleep. Then, it calculates Post-Pre.
%
% [R] = reactivation_strength(patterns , cond , SpikeTrain , is , th , type)
%
% --- Inputs ---
% patterns: matrix with the weigths for each cell in each assembly.
%           Structure: Single-Units x Assemblies (rows x column)
%
% cond: logical, to select which pattern will be used (1, include / 0, not include)
%       Row vector with the same elements as number of assemblies.
%       Example,
%       patterns:     A1    A2    A3
%                 SU1 0.6   0.1   0.1       ---> cond: 0   1   0
%                 SU2 0.3   0.6   0.1       In this case, only the second
%                 SU3 0.1   0.3   0.3       assembly (A2) will be selected.
%                 SU4 0.1   0.1   0.6
%
% SpikeTrain: matrix, First column, time bins, Rest columns Spike Counts
%             Example:   t    SU1   SU2   SU3   SU4
%                        0     3     5     1     3    This Spike Train was
%                       0.1    1     2     1     0    constructed using a
%                       0.2    3     5     1     3    0.1 sec time bin.
%                       ...   ...   ...   ...   ...
%
% Is: structure, it contains the periods of time (in sec) of interest.
%     Struture:
%                 Is.baseline     Should be a vector with the same length
%                 Is.aversive     of the SpikeTrain with logical values to
%                 Is.reward       include or not the time bin.
%                 Is.runaversive
%                 Is.runreward
%
%     Example:
%             is.aversive:
%                         b1  b2  b3  b4  b5  b6  b7  b8  b9
%                         0   0   1   1   1   1   0   0   0
%
% th: float/int, threshold value for peaks detection. Note that this
%     function zscored the assemblies activity.
%
% type: String, 'A' if is Aversive assemblies 'R' if they are Reward
%
% config: int, if is 1 Aversive occurs first, if is 2, Reward was first
%
% --- OUTPUT ---
% R: column vector storing all the Reactivation Strength values for sleep
%    and awake periods.
%
%    S1 ---- *C1*  ---- S2 ----  ---- C2 ---- S3
%
%    * Assemblies detected during this session. It could do the same if C2
%    is the condition of interest. Just be sure you are choosing correctly
%    the type and config inputs.
%
%    Structure:  mean(S1)-mean(S2)   mean(C1)   mean(C2)
%                         R1           A11        A21
%                         R2           A12        A22
%                         R3           A13        A23
%                         ...          ...        ...
%
% requirments:
%       assembly_activity.m from Lopes-dos-Santos et al 2013 (*, see below)
%
% Morici Juan Facundo 09/2023
% also see DOI: 10.1016/j.neuron.2016.10.020
% also see DOI:10.1016/j.jneumeth.2013.04.010 (*)

bins = SpikeTrain(:,1);
spks = SpikeTrain(:,2:end);

a = assembly_activity(patterns(:,cond) , spks');
a = zscore(a,1,2);
R = [];
for i = 1:size(a,1)
    [pks.baseline,loc.baseline] = findpeaks(a(i,Is.baseline),bins(Is.baseline),'MinPeakHeight',th);
    [pks.reward,loc.reward] = findpeaks(a(i,Is.reward),bins(Is.reward),'MinPeakHeight',th);
    [pks.aversive,loc.aversive] = findpeaks(a(i,Is.aversive),bins(Is.aversive),'MinPeakHeight',th);
    
    [pks.runreward,loc.reward] = findpeaks(a(i,Is.runreward),bins(Is.runreward),'MinPeakHeight',th);
    [pks.runaversive,loc.aversive] = findpeaks(a(i,Is.runaversive),bins(Is.runaversive),'MinPeakHeight',th);
    
    if type == 'A' % check if is aversive assembly
        if config == 1
            R = [R ; (mean(pks.aversive) - mean(pks.baseline)) mean(pks.runaversive) mean(pks.runreward)];
            clear tmp
        else
            R = [R ; (mean(pks.aversive) - mean(pks.reward)) mean(pks.runaversive) mean(pks.runreward)];
            clear tmp
        end
        clear pks loc
        
    elseif type == 'R' % check if is reward assembly
        if config == 2
            R = [R ; (mean(pks.reward) - mean(pks.baseline)) mean(pks.runreward) mean(pks.runaversive)];
            clear tmp
        else
            R = [R ; (mean(pks.reward) - mean(pks.aversive)) mean(pks.runreward) mean(pks.runaversive)];
            clear tmp
        end
    end
end

end