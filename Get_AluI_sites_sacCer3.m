% Find all AluI sites in the yeast genome
% Yeast DNA sequence is stored in the file sacCer3_genome.mat

load('sacCer3_genome.mat', 'genome')
noChr = numel(genome);

% Initialization
Chr_AluI = [];
Loc_AluI = [];
TotalNoSites = 0;

% Find the AluI motifs, AGCT, in each chromosome, and keep the locations
for chr = 1 : noChr
    ind = strfind(upper(genome(chr).Seq), 'AGCT');
    noSitesPerChr = numel(ind);
    
    Chr_AluI = [Chr_AluI; chr*ones(noSitesPerChr, 1)];
    Loc_AluI = [Loc_AluI; ind';];
    TotalNoSites = TotalNoSites + noSitesPerChr;
end

save('AluI_sites_sacCer3.mat', 'Chr_AluI', 'Loc_AluI', 'TotalNoSites');