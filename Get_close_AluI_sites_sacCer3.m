% AluI sites that are <50bp apart will produce short fragments that are not
% sequenced. Correct this problem using the function:
% CorrectedRatios = CorrectRatios_sacCer3(Ratios, Sites_with_problems_Right, Sites_with_problems_Left, Sites_with_problems_Both_Sides);

%%
load('sacCer3_genome.mat', 'genome')
noChr = numel(genome);
chrLen = [genome.chrLen];

Motif_Filter = cell(1, noChr);

for chr = 1:noChr
    Motif_Filter{chr} = zeros(1, chrLen(chr));
    ind = strfind(upper(genome(chr).Seq), 'AGCT');
    
    noSites = numel(ind);
    for s = 1:noSites
        Motif_Filter{chr}(ind(s)+[1,2]) = 1;
    end
end

%% Check for motifs that are very close to each other
oneSitePattern = [0,1,1,0]; 

d = 0;
for chr = 1:noChr
    Sites_with_problems_Right{chr} = zeros(1, chrLen(chr));
    Sites_with_problems_Left{chr} = zeros(1, chrLen(chr));
    
    ind = strfind(Motif_Filter{chr}, [oneSitePattern, oneSitePattern]);
    Sites_with_problems_Right{chr}(ind+2) = 1;
    Sites_with_problems_Left{chr}(ind+d+5) = 1;
end

for d = 1:46
    for chr = 1:noChr
        ind = strfind(Motif_Filter{chr}, [oneSitePattern, zeros(1, d), oneSitePattern]);
        Sites_with_problems_Right{chr}(ind+2) = 1;
        Sites_with_problems_Left{chr}(ind+d+5) = 1;
    end
end

for chr = 1:noChr
    Sites_with_problems_Both_Sides{chr} = zeros(1, chrLen(chr));
    
    ind_Right = strfind(Sites_with_problems_Right{chr}, [0, 1]);
    ind_Left = strfind(Sites_with_problems_Left{chr}, [1, 0]);
    
    ind = intersect(ind_Left, ind_Right);
    
    Sites_with_problems_Both_Sides{chr}(ind) = 1;
    Sites_with_problems_Both_Sides{chr}(ind+1) = 1;
end

save('AluI_sites_closer_than_50bp.mat', 'Sites_with_problems_Right', 'Sites_with_problems_Left', 'Sites_with_problems_Both_Sides');