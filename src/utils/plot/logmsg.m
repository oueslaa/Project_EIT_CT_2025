function logmsg(fmt, varargin)
fprintf('[%s] ', datestr(now, 'HH:MM:SS'));
fprintf(fmt, varargin{:}); fprintf('\n');
end