function result = sphericalZYXSubspace(q)
    % Extract q1, q2, and q3 from the 3x1 vector q
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    
    % Calculate sine and cosine values for q1, q2, and q3
    s1 = sin(q1);
    c1 = cos(q1);
    s2 = sin(q2);
    c2 = cos(q2);
    s3 = sin(q3);
    c3 = cos(q3);
    
    % Construct the 3x3 matrix S
    S = [  -s2,        0,         1;
            c2*s3,     c3,       0;
            c2*c3,    -s3,       0 ];
    
    % Create a 6x3 matrix with the top 3 rows as zeros and bottom 3 rows as S
    result = [zeros(3, 3); S];
end
