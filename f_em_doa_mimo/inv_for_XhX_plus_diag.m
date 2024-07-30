function result = inv_for_XhX_plus_diag(   xi, X, lambda   )
% return the result of inv(   xi*X'*X + diag(lambda)   )

if size(X,2) <= size(X,1) % a direct inverse would be ok
    result = inv(   xi*(X'*X) + diag(lambda)   );
else
    inv_Lambda = diag(   1./lambda   );
    result = inv_Lambda - xi*inv_Lambda*X'* (   (eye(size(X,1))+xi*X*inv_Lambda*X')...
                                \X*inv_Lambda   );
end


end