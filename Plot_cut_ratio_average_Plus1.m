function Plot_cut_ratio_average_Plus1(AluI_cleavages_filename)

load('sacCer3_genome.mat', 'genome')
noChr = numel(genome);
chrLen = [genome.chrLen];
chrName = {genome.chrName};

Filter = cell(1, noChr);

for chr = 1:noChr
    Filter{chr} = nan(1, chrLen(chr));
    ind = strfind(upper(genome(chr).Seq), 'AGCT');
    
    noSites = numel(ind);
    for s = 1:noSites
        Filter{chr}(ind(s)+[1,2]) = 1;
    end
end

load(AluI_cleavages_filename, 'Cuts', 'Occ')
rawRatios = cellfun(@(x,y) x./y, Cuts, Occ, 'un', 0);
load('AluI_sites_closer_than_50bp.mat', 'Sites_with_problems_Left', 'Sites_with_problems_Right', 'Sites_with_problems_Both_Sides')
correctedRatios = CorrectRatios_sacCer3(rawRatios, Sites_with_problems_Right, Sites_with_problems_Left, Sites_with_problems_Both_Sides);
filteredRatios = cellfun(@(x,y) x.*y, correctedRatios, Filter, 'un', 0);

%% Align all TSSs
load('Annotations_sacCer3.mat', 'Chr', 'Plus1', 'Watson')
% Eliminate genes for which the +1 nuc. annotations are missing
Chr(isnan(Plus1)) = [];
Watson(isnan(Plus1)) = [];
Plus1(isnan(Plus1)) = [];

Ref = Plus1;
beforeRef = 1000;
afterRef = 1000;

noGenes = numel(Ref);
AlignedCutRatios = nan(noGenes, 1 + beforeRef + afterRef);
for g = 1:noGenes
    if Watson(g)
        leftEdge = max([Ref(g) - beforeRef, 1]);
        rightEdge = min([Ref(g) + afterRef, chrLen(Chr(g))]);
        AlignedCutRatios(g, beforeRef + 1 - (Ref(g) - leftEdge)...
            : beforeRef + 1 + (rightEdge - Ref(g))) = ...
            filteredRatios{1,Chr(g)}(leftEdge : rightEdge);
    else
        leftEdge = max([Ref(g) - afterRef, 1]);
        rightEdge = min([Ref(g) + beforeRef, chrLen(Chr(g))]);
        AlignedCutRatios(g, beforeRef + 1 - (rightEdge - Ref(g))...
            : beforeRef + 1 + (Ref(g) - leftEdge)) = ...
            fliplr(filteredRatios{1,Chr(g)}(leftEdge : rightEdge));
    end
end
meanCutRatios = smooth(nanmean(AlignedCutRatios), 21);

%% Plot figure
figure
plot([-beforeRef:afterRef]/1000, meanCutRatios);
ylabel('Cut fraction')
xlabel('Position relative to +1 nuc. (kb)')
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', 'Average_cut_ratios_+1.eps');