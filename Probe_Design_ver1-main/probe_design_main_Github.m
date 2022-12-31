function probe_design_main_Github(FASTA_dir, Probelist, out_dir, max_probe_number)

%% Make gene list from folder structure

% FASTA_dir should has folders for each gene. The folders for gene should have one text FASTA file.
% FASTA_dir = 'E:\Shigeaki\Probe_design\Probes';

% Probelist is a table of [Gene, Length, hairpin, Start1, End1, Start2, End2 .... StartN, EndN];
% [Gene_name, length_of_gene, hairpin(B1-B5), the area you want to remove showed by pairs of (start-end)]

% set max_probe_number = -1 if you don't want to reduce the number of probes

% Error log
% 1. FASTA_dir should be the folder which contain gene folder such as "Mouse_Ang2" etc. directly
% 2. Mismatch of gene list between Probelist and folders in FASTA_dir
% 3. Hairpin type is strange (B14 etc.)

file_list = struct2table(dir(FASTA_dir));

folder_list = file_list(file_list.isdir, :);
folder_list = folder_list(3:end, :);
num_of_gene = length(folder_list{:,1});

%% preallocation

All_probe_list= {}; 
hairpin_type = strings(1, num_of_gene);
Fasta_header = strings(1, num_of_gene);
Fasta_seq = strings(1, num_of_gene);

%% Probe sequence design (Tiling)

for i=1:num_of_gene

    target_name = folder_list.name{i};
    check_vector = Probelist.Gene == target_name;
    target_row = Probelist(check_vector, :);  % extract target_row
    hairpin_type(i) = char(target_row.hairpin);
    
    gene_dir = [FASTA_dir '\' char(target_name)]; 
    gene_fastatxt = dir([gene_dir '\' '*.fasta']);
    fname = [gene_fastatxt.folder '\' gene_fastatxt.name]; 
    
    Fasta_struct = fastaread(fname);
    Fasta_cell = struct2cell(Fasta_struct);
 
    length_sequence = target_row{1, 4:end};

    if isstring(length_sequence) ==1 
        length_sequence = str2double(length_sequence);
    end
        
    length_sequence(isnan(length_sequence)) = [];

    trial_sequence = length(length_sequence) ./ 2;
    
    Fasta_seq(i) = Fasta_cell{2, 1};
    Fasta_header(i) = Fasta_cell{1, 1};
    
    Probe_list = [];
    
        for s=1:trial_sequence
    
            start_seq = target_row{1, (2 * s + 2)}; 
            end_seq =  target_row{1, (2 * s + 3)}; 
                     
            Fasta_seq_part = char(Fasta_seq(i));
            Fasta_seq_use = Fasta_seq_part(start_seq:end_seq);
            
            [Export_list] = probe_design_version3_HCR(Fasta_seq_use, hairpin_type(i), start_seq);

            Probe_list = cat(1, Probe_list, Export_list);
              
        end
    
     All_probe_list{i} = Probe_list;
     
end


%% save data
    
    target_list = string(folder_list.name);
    hairpin_type_list = hairpin_type;
    probe_pair_number = zeros(length(target_list),1);
        
for i=1:num_of_gene
    
      Odd_header = strcat('Odd probe_', hairpin_type_list(i));
    Even_hearder = strcat('Even probe_', hairpin_type_list(i));
    Export_list_Header = cat(2, Odd_header, 'Start', 'End', Even_hearder, 'Start', 'end');
    
    Fasta_header_cell = {Fasta_header(i), '', '', '', '' , ''};  % to make compatible to cat
    Fasta_header_string = string(Fasta_header_cell);
    
    probe_pair_number(i) = size(All_probe_list{i}, 1);
    
    Export_list = cat(1, Fasta_header_string, Export_list_Header, All_probe_list{i});
      
    writematrix(Export_list, [out_dir '\Probe_list_'  folder_list.name{i} '.csv']);
    
    save ([out_dir '\All_probe_list.mat'], 'All_probe_list');  
        
end

%% Save summary

     Header_table = ["Target_name", "Sequence_length", "Probe_pair_number"];
     ProbeSequence_length_list = probe_pair_number .* 90 ; % one pair length is 45*2 = 90
     Summary_table_pre = cat(2, target_list, ProbeSequence_length_list, probe_pair_number);
     
     Summary_table = cat(1, Header_table, Summary_table_pre);
     writematrix(Summary_table, [out_dir '\' 'Probe_design_summary.csv']);

%% Make excel sheet for opool oligo
     
    opool_header = ["Pool name", "Sequence"];
        
    opool_hairpin_list = transpose(hairpin_type_list); % these two should be from the same source
    opool_name_list = strcat(target_list, '_', opool_hairpin_list);  % these two should be from the same source
    
    opool_excel_list = {};
   
    for i=1:length(opool_name_list)
        
        Gene_probe_list = All_probe_list{i};
        
        probe_list_length = size(Gene_probe_list, 1);
        text_store = strings(probe_list_length .* 2, 1);
        text_store(:) = opool_name_list(i);
        
        opool_sequence_list = cat(1, Gene_probe_list(:,1), Gene_probe_list(:,4));
        opool_export_list = cat(2, text_store, opool_sequence_list);

        opool_excel_list = [opool_excel_list; opool_export_list];
        
    end

     opool_excel_list_save = cat(1, opool_header, opool_excel_list);
     
     save ([out_dir '\opool_excel_list_save.mat'], 'opool_excel_list_save');  
     
     xlswrite([out_dir '\' 'opool_probe_sheet.xlsx'], opool_excel_list_save);
     %writematrix(opool_excel_list_save, [out_dir '\' 'opool_probe_sheet.xlsx']); % writematrix has a bug.
    
     
%% reduce probe size if necessary
    if max_probe_number > 0       
        All_probe_list_reduced = probe_set_reduce_ver1(All_probe_list, max_probe_number);
        mkdir([out_dir '\Reduced']);
    
    %% save data
    
              
        for i=1:num_of_gene
        
            probe_pair_number(i) = size(All_probe_list_reduced{i}, 1);
            Export_list = cat(1, Fasta_header_string, Export_list_Header, All_probe_list_reduced{i});
      
            writematrix(Export_list, [out_dir '\Reduced\Probe_list_'  folder_list.name{i} '.csv']);
    
            save ([out_dir '\Reduced\All_probe_list_reduced.mat'], 'All_probe_list_reduced');  
        
        end

%% Save summary

       Header_table = ["Target_name", "Sequence_length", "Probe_pair_number"];
       ProbeSequence_length_list = probe_pair_number .* 90 ; % one pair length is 45*2 = 90
       Summary_table_pre = cat(2, target_list, ProbeSequence_length_list, probe_pair_number);
     
       Summary_table = cat(1, Header_table, Summary_table_pre);
       writematrix(Summary_table, [out_dir '\Reduced\' 'Probe_design_summary_reduced.csv']);

%% Make excel sheet for opool oligo
     
    opool_header = ["Pool name", "Sequence"];
        
    opool_hairpin_list = transpose(hairpin_type_list); % these two should be from the same source
    opool_name_list = strcat(target_list, '_', opool_hairpin_list);  % these two should be from the same source
    
    opool_excel_list = {};
   
    for i=1:length(opool_name_list)
        
        Gene_probe_list = All_probe_list_reduced{i};
        
        probe_list_length = size(Gene_probe_list, 1);
        text_store = strings(probe_list_length .* 2, 1);
        text_store(:) = opool_name_list(i);
        
        opool_sequence_list = cat(1, Gene_probe_list(:,1), Gene_probe_list(:,4));
        opool_export_list = cat(2, text_store, opool_sequence_list);

        opool_excel_list = [opool_excel_list; opool_export_list];
        
    end

     opool_excel_list_save = cat(1, opool_header, opool_excel_list);
     
     save ([out_dir '\Reduced\opool_excel_list_save_reduced.mat'], 'opool_excel_list_save');  
     
     xlswrite([out_dir '\Reduced\' 'opool_probe_sheet_reduced.xlsx'], opool_excel_list_save);
     %writematrix(opool_excel_list_save, [out_dir '\' 'opool_probe_sheet.xlsx']); % writematrix has a bug.
               
        
    end
    
     
     
end     