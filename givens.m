function [R, QTy] = givens(A, b) % Givens-rotation QR for overdetermined least squares
    [m,n] = size(A);
    R = A;
    QTy = b;
    for j = 1:n
        for i = m:-1:(j+1)
            a = R(j,j);
            g = R(i,j);
            if g ~= 0
                r = hypot(a, g);
                c = a / r;
                s = -g / r;
                Gblock = [c -s; s c];
                R([j i], j:n) = Gblock * R([j i], j:n);
                QTy([j i]) = Gblock * QTy([j i]);
            end
        end
    end
    R   = R(1:n, :);
    QTy = QTy(1:n);
end