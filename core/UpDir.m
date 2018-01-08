function dir = UpDir(dir,n)
    for i = 1:n
        dir = fileparts(dir);
    end
end
