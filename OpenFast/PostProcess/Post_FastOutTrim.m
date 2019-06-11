function fastout = Post_FastOutTrim(fastout,t0,t1)
% This function trims the solution of a fast simulation to the desired
% time-steps. 
%
% Inputs: fastout - structure of FAST output data
%         t0      - start time
%         t1      - end time
%
% Nikhar Abbas - February 2019

% Load fields
if isfield(fastout,'headers')
    fastout = rmfield(fastout,'headers');
    fastout = rmfield(fastout,'units');
end

fields = fieldnames(fastout);

% find time indices
t0_ind = find(fastout.Time == t0);
t1_ind = find(fastout.Time == t1);

if isempty(t0_ind) || isempty(t1_ind)
    error('The desired time indices do not exist in this simulation')
end

for f_ind = 1:length(fields)
    if isnumeric(fastout.(fields{f_ind}))
        fastout.(fields{f_ind}) = fastout.(fields{f_ind})(t0_ind:t1_ind);
    end
end

fastout.Time = fastout.Time - t0;




end