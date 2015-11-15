function out = l0norm(x, f)

    out = size(x, 1);
    maxX = max(x);
    
    for i=1:size(x, 1)
        
        if abs(x(i)) < maxX*f
            out = out - 1;
        end
        
    end

end