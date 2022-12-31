function [output_list, target_seq_column] = target_selection_ver1(sequence_length, rem_start, rem_end)

% This script remove the sequence area which overlapping with other genes
% example, rem_start = [5;20;400];
% example, rem_end = [25;40;600]; 
% then result is that sequence we can target is start[1;41;601], end[4;399;800], length[4; 359; 200]
% sequence_length is the lengeth of mRNA sequence (1000 for example)

position = true(sequence_length,1); % position is logical (now all true and change the seq to false later)

for i=1:length(rem_start)
    position(rem_start(i):rem_end(i)) = 0; % change to false if sequence is removed
    
end

a = find(position==1); % get sequence area which is true

D = diff([0;diff(a)==1;0]); % get subtraction from next position and convert =1 (start position), =-1 (end position)
first = a(D>0);  % to get start position
last = a(D<0);   % to get end position
target_seq = [first last];
target_length = 1 + find(D<0) - find(D>0); % Calculate sequence length

output_list = [target_seq target_length];
target_seq_column = reshape(target_seq.',1,[]);  %  array of Start1, End1, Start2, End2, ..

end
