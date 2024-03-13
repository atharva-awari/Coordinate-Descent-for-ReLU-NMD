#This function computes an update for a single entry of the matrix H, say H(i,j)


using LinearAlgebra

function updateH(X,W,H,WH)
    
    (r,n) = size(H)

    Threads.@threads for j=1:n
        #c = X[:,j]
        c=view(X,:,j)
    
        for i = 1:r
            a = W[:,i]
            b = WH[:,j] - W[:,i]*H[i,j]
            Hijold = H[i,j]
            H[i,j] = findmin_cd(a,b,c)
            @views WH[:,j]  .+= (H[i,j] - Hijold) .* W[:,i] 
        end

    end
    return nothing
end
