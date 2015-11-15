function xio = xi(x, t)

    xio = zeros(size(x, 2), 1);
    
    xio(x < -t) = -1;
    xio(x == -t) = -0.5;
    xio(-t < x & x < t) = 0;
    xio(x == t) = 0.5;
    xio(x > t) = 1;

end