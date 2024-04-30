function vecsk = skewsymmatrix(vec)

% Take a 3 x 1 vector "vec" and turn it into a skew symmetric matrix vecsk

% This is fine

vecsk = [0 -vec(3) vec(2); 
         vec(3) 0 -vec(1); 
        -vec(2) vec(1) 0];

end