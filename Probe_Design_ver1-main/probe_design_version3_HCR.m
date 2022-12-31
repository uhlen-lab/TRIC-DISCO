function [Export_list] = probe_design_version3_HCR(FASTA_sequence, hairpin_type, seq_position_FASTA)

%% Hairpin sequence

B1_I1_odd = "gAggAgggCAgCAAACgg";
B2_I1_odd = "CCTCgTAAATCCTCATCA";
B3_I1_odd = "gTCCCTgCCTCTATATCT";
B4_I1_odd = "CCTCAACCTACCTCCAAC";
B5_I1_odd = "CTCACTCCCAATCTCTAT";

B1_I1_even = "gAAgAgTCTTCCTTTACg";
B2_I1_even = "ATCATCCAgTAAACCgCC";
B3_I1_even = "CCACTCAACTTTAACCCg";
B4_I1_even = "TCTCACCATATTCgCTTC";
B5_I1_even = "CTACCCTACAAATCCAAT";

Spacer = "AA";

if hairpin_type == "B1" 
    hairpin_odd = B1_I1_odd; hairpin_even = B1_I1_even;
elseif hairpin_type == "B2"
    hairpin_odd = B2_I1_odd; hairpin_even = B2_I1_even;
elseif hairpin_type == "B3"
    hairpin_odd = B3_I1_odd; hairpin_even = B3_I1_even;
elseif hairpin_type == "B4"
    hairpin_odd = B4_I1_odd; hairpin_even = B4_I1_even;
elseif hairpin_type == "B5"
    hairpin_odd = B5_I1_odd; hairpin_even = B5_I1_even;
else
    disp('Hairpin type is not collect. Write from "B1" to "B5".')   
end

% 5.8S ribosomal RNA (mouse)
% %FASTA_sequence = char(['GACTCTTAGCGGTGGATCACTCGGCTCGTGCGTCGATGAAGAA' ....
%                 'CGCAGCTAGCTGCGAGAATTAATGTGAATTGCAGGACACATTGATCATC' ....
%                 'GACACTTCGAACGCACTTGCGGCCCCGGGTTCCTCCCGGGGCTACGCCTGTCTGAGCGTCGGTTG']);

%% Sequence search parameter

RNA_sequence = FASTA_sequence;
RNA_length = length(RNA_sequence);

Probe_pair_number = RNA_length / (25 + 2 + 25 + 2); % (odd probe + gap1 even probe + gap2)

%% Preallocation

odd_probe  = strings(floor(Probe_pair_number), 1);
even_probe = strings(floor(Probe_pair_number), 1);
seq_position = zeros(floor(Probe_pair_number), 1);
%polyA_tail = repmat('A', 1, 25);

    for i = 1:Probe_pair_number
    
        seq_position(i) = 1 + (i-1) * 52;
        Target_seq_odd  = RNA_sequence (seq_position(i):(seq_position(i) + 25 - 1));
        Target_seq_even = RNA_sequence ((seq_position(i) + 27):(seq_position(i) + 52 - 1));
        
        check_polyA = [find(Target_seq_even == 'T'), find(Target_seq_even == 'C'), find(Target_seq_even == 'G')];        
        
        if sum(check_polyA) == 0 
        % if Target seq is not equal to polyA_tail, use the this sequence to probes
        
        seq_position(i) = [];
        odd_probe(i) = [];
        even_probe(i) = [];
        
        else
        
        Probe_seq_odd  = seqrcomplement(Target_seq_odd);
        Probe_seq_even = seqrcomplement(Target_seq_even);
        
        odd_probe(i)  = hairpin_odd + Spacer + Probe_seq_odd;
        even_probe(i) = Probe_seq_even + Spacer + hairpin_even;
        
        end
        
    end

    
%% Make export list
    
    seq_position_probe = seq_position + seq_position_FASTA - 1;
    seq_position_odd = [seq_position_probe, seq_position_probe + 25 - 1];
    seq_position_even = [seq_position_probe + 27, seq_position_probe + 52 - 1];

    Export_list = cat(2, odd_probe, seq_position_odd, even_probe, seq_position_even);

end
