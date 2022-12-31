function All_probe_list_reduced = probe_set_reduce_ver1(All_probe_list, max_probe_number)

probelist_number = cellfun(@(c) size(c,1),All_probe_list);
%All_probe_list_reduced = cell(size(All_probe_list));
All_probe_list_reduced = All_probe_list;
Approx_mat = [1	0.5	0.333333333	0.25	0.2	0.166666667	0.142857143	0.125 0.111111111	0.1 ];

% [1	2	3	4	5	6	7	8	9	10]
% [1	0.5	0.333333333	0.25	0.2	0.166666667	0.142857143	0.125 0.111111111	0.1]


for i=1:length(All_probe_list)
        
    if  probelist_number(i) > max_probe_number 
    
        remove_number = probelist_number(i)- max_probe_number;        % calculate how many probe set need to be removed
        percent_value = remove_number / probelist_number(i);               % calculate how many every probe set should be removed
       
        [c index] = min(abs(Approx_mat-percent_value));         
        All_probe_list_reduced{1,i}(1:index:end,:) = [];   % every index-th data changed to empty 
    
        if index == 1
        
            selection_percent_value = 1 - percent_value;
            selection_value = ceil(1 / selection_percent_value);
            
            All_probe_list_reduced{1,i} = All_probe_list{1,i}(1:selection_value:end,:);   % every index-th data changed to empty 
    
            
        end
    end 
end 