function gout = gFunc(z, t)

    gout = 0;
    
    for i=1:size(z, 1)
        
        gout = gout + max(z(i)-t, 0) + max(-z(i)-t, 0);
        
    end

end