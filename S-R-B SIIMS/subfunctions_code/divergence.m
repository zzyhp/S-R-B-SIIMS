function div_p = divergence(p1, p2)
    [N,M]=size(p1);
    

    p1_x = zeros(N,M);
    for i = 2:N-1
        p1_x(i, :) = (p1(i+1, :) - p1(i-1, :)) / 2;
    end


    p2_y = zeros(N,M);
    for j = 2:M-1
        p2_y(:, j) = (p2(:, j+1) - p2(:, j-1)) / 2;
    end

    div_p = p1_x + p2_y;
end