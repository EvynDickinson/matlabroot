function [ts,T,t0]=f(file_name,channel_name)

% ts is an array of times in seconds, where t==0 means right at the start
%    of the recording.  It's a col vector.
% T is the total duration of the recording, where T==dt*N, where dt is the
%   sampling interval in seconds, and N is the number of samples.  (Note
%   that the time of the first sample is t0, and the time of the last sample
%   is t0+(N-1)*dt).  T is in seconds.
% t0 is the start time of the recording.  

% load the data
fid=fopen(file_name,'r');
channel_number=channel_name_to_channel_number(fid,channel_name);
[ts,h]=SONGetChannel(fid,channel_number);
file_header=SONFileHeader(fid);
N=file_header.maxFTime+1;  % number of samples, since first sample
                           % is at time 0 in clock ticks
T=SONTicksToSeconds(fid,N);  % in s
t0=0;  % s
fclose(fid);
