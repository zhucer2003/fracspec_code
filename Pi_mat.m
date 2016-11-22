function P = Pi(m, r, s, Pi1)

if ( any(s) )
    error
end

m1 = ceil(m) + .5;
m2 = floor(m) + 1;
P = blkdiag(Pi1(r, m1), Pi1(r, m2));

% if ( round(m) == m )
%     P = blkdiag(Pi1(r, m+.5), Pi1(r, m));
% else
%     P = blkdiag(Pi1(r, m+1), Pi1(r, m+.5));    
% end

end