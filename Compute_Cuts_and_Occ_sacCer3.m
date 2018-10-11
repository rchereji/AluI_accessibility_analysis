function [Cuts, Occ] = Compute_Cuts_and_Occ_sacCer3(bamFilename)

% Minimum and maximum DNA fragment sizes to be used in the analysis (we keep all reads)
Lmin = 0;
Lmax = 5000;

% Chromosome names and lengths (sacCer3)
chrName = {'chrI';'chrII';'chrIII';'chrIV';'chrV';'chrVI';'chrVII';'chrVIII';...
    'chrIX';'chrX';'chrXI';'chrXII';'chrXIII';'chrXIV';'chrXV';'chrXVI'};
chrLen = [230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643,...
    439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066];

% Total number of chromosomes
noChr = numel(chrName);

% Initialize the distribution of left/right ends of all DNA fragments
LeftEnd = cell(1, noChr);
RightEnd = cell(1, noChr);

% Initialize the occupancy profile (number of fragments overlapping each bp)
Occ = cell(1, noChr);

% Process all chromosomes
for chr = 1 : noChr
    
    % Create a BioMap object
    bm = BioMap(bamFilename, 'SelectReference', chrName{chr});
    
    % Discard the reads with low quality
    bm = getSubset(bm, getMappingQuality(bm) > 0);
    
    % Make sure that no read falls outside the known chromosome size
    bm = getSubset(bm, getStop(bm) <= chrLen(chr));
    
    % Get the reads mapped on the Watson strand
    Indices = filterByFlag(bm, 'pairedInMap', true, 'strandQuery', 0);
    bm_filtered_Watson = getSubset(bm, Indices);
    
    % Get the reads mapped on the Crick strand
    Indices = filterByFlag(bm, 'pairedInMap', true, 'strandQuery', 1);
    bm_filtered_Crick = getSubset(bm, Indices);
    
    % Keep only the properly paired reads, which have identical headers
    [~, Watson_Idx, Crick_Idx] = intersect(bm_filtered_Watson.Header, bm_filtered_Crick.Header);

    % Compute fragment length histogram
    leftBP = getStart(bm_filtered_Watson, Watson_Idx);
    rightBP = getStop(bm_filtered_Crick, Crick_Idx);
    
    % Compute the fragment lengths
    fragmentLengths = rightBP - leftBP + 1;
    
    % Optional size selection; we keep all reads (0 bp <= L <= 5000 bp)
    goodInd = ((fragmentLengths >= Lmin) & (fragmentLengths <= Lmax));
    leftBP = leftBP(goodInd);
    rightBP = rightBP(goodInd);
   
    % Construct distributions for left/right ends of DNA fragments
    uniqueLeftBP = unique(leftBP);
    [~, Index] = ismember(leftBP, uniqueLeftBP) ;
    NumberUniqueLeftBP = histc(Index, 1:numel(uniqueLeftBP)) ;
    
    LeftEnd{1,chr} = zeros(1, chrLen(chr));
    LeftEnd{1,chr}(uniqueLeftBP) = NumberUniqueLeftBP;
    
    uniqueRightBP = unique(rightBP);
    [~, Index] = ismember(rightBP, uniqueRightBP) ;
    NumberUniqueRightBP = histc(Index, 1:numel(uniqueRightBP)) ;
    
    RightEnd{1,chr} = zeros(1, chrLen(chr));
    RightEnd{1,chr}(uniqueRightBP) = NumberUniqueRightBP;
    
    % Compute the occupancy profile (stack all the reads)
    OccDerivative = [LeftEnd{1,chr}, 0]; % add 1 position for the case when some reads have the right end exactly at the end of chr
    OccDerivative(2:end) = OccDerivative(2:end) - RightEnd{1,chr};
    tmp = cumsum(OccDerivative);
    Occ{1,chr} = tmp(1 : end-1);
    
    fprintf('Chromosome %s was processed.\n', chrName{chr})
end

% Compute the number of AluI cuts at each position along the genome
Cuts = cellfun(@(x,y) x+y, LeftEnd, RightEnd, 'UniformOutput', 0);
