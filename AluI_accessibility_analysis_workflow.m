% Compute the following 2 profiles (cell arrays, one cell for each chromosome):
%   1) Cuts - number of cuts (fragment ends) at each genomic position 
%   2) Occ  - number of fragments (cut or uncut) that cover each position
%             (i.e. the usual occupancy profile obtained by stacking all reads)
filename = '111P_4_AluI_400U.bam';
[Cuts, Occ] = Compute_Cuts_and_Occ_sacCer3(filename);

% Optional: save results in a separate file
save('AluI_cleavages_111P_4_AluI_400U.mat', 'Cuts', 'Occ')

% Write the genome-wide cut ratio (Cuts / Occ) into a WIG file
Write_cut_ratios_to_WIG_file('AluI_cleavages_111P_4_AluI_400U.mat')
% The WIG file can be further converted into BIGWIG track using the
% instructions from https://genome.ucsc.edu/goldenpath/help/bigWig.html

% Plot the cut ratios near TSS
Plot_cut_ratio_average_TSS('AluI_cleavages_111P_4_AluI_400U.mat')

% Similarly, one can align the +1 nucleosomes and plot the average cut
% ratio relative to the typical +1 nuc. positions
Plot_cut_ratio_average_Plus1('AluI_cleavages_111P_4_AluI_400U.mat')
