function y = istrue(x)
y = ~isempty(x) && ( (islogical(x) && any(x)) || (isnumeric(x) && any(x(:)~=0)) );
end
