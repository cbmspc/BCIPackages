function hex = convert_to_double_memory_format (num)

d = typecast(cast(num, 'double'), 'uint8'); 
if strcmpi('little', 'big') 
   d = d(end:-1:1); 
end 
hex = reshape(dec2hex(d).', 1, []);
