function Write_cut_ratios_to_WIG_file(AluI_cleavages_filename)
% Function to write the cut ratios to a WIG file that can be viewed in IGV
% INPUT: AluI_cleavages_filename - contains the cell arrays generated by
% Compute_Cuts_and_Occ_sacCer3.m, e.g. 'AluI_cleavages_111P_4_AluI_400U.mat';

load('sacCer3_genome.mat', 'genome')
noChr = numel(genome);
chrLen = [genome.chrLen];
chrName = {genome.chrName};

Motif_Filter = cell(1, noChr);

for chr = 1:noChr
    Motif_Filter{chr} = zeros(1, chrLen(chr));
    ind = strfind(upper(genome(chr).Seq), 'AGCT');
    
    noSites = numel(ind);
    for s = 1:noSites
        Motif_Filter{chr}(ind(s)+[1,2]) = 1;
    end
end

load('AluI_sites_closer_than_50bp.mat', 'Sites_with_problems_Left', 'Sites_with_problems_Right', 'Sites_with_problems_Both_Sides')

%%
load(AluI_cleavages_filename, 'Cuts', 'Occ')
rawRatios = cellfun(@(x,y) x./y, Cuts, Occ, 'un', 0);

% Correct the problem produced by short fragments that are not sequenced
correctedRatios = CorrectRatios_sacCer3(rawRatios, Sites_with_problems_Right, Sites_with_problems_Left, Sites_with_problems_Both_Sides);
filteredRatios = cellfun(@(x,y) x.*y, correctedRatios, Motif_Filter, 'un', 0);

for c = 1:noChr
    filteredRatios{c}(isnan(filteredRatios{c})) = 0;
end

WIG_filename = strrep(AluI_cleavages_filename, '.mat', '.wig');
WIG_filename = strrep(WIG_filename, 'AluI_cleavages_', 'Cut_ratios_');

fileID = fopen(WIG_filename,'w');
for chr = 1:noChr
    fprintf(fileID,'fixedStep  chrom=%s  start=1  step=1\n', chrName{chr});
    fprintf(fileID,'%0.3f\n', filteredRatios{chr}(:));
end
fclose(fileID);

