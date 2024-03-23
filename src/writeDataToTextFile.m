function writeDataToTextFile(filename, header, array1, array2, array3)
    % Check if the arrays have the same length
    if ~isequal(numel(array1), numel(array2), numel(array3))
        error('Arrays must have the same length.');
    end
    
    % Open the file for writing
    fid = fopen(filename, 'w');
    
    % Check if the file opened successfully
    if fid == -1
        error(['Could not open file ' filename ' for writing.']);
    end
    
    % Write header
    fprintf(fid, '%s\t%s\t%s\n', header{1}, header{2}, header{3});
    
    % Write data from arrays to the file
    for i = 1:numel(array1)
        fprintf(fid, '%e\t%e\t%e\n', array1(i), array2(i), array3(i));
    end
    
    % Close the file
    fclose(fid);
end