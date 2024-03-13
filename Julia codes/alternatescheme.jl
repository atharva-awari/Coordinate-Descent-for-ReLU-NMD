# Computes an approximate solution of the following non linear matrix
# decomposition problem (NMD)
#
#       min_{W,H} ||X-max(0,WH)||_F^2
#
# using coordinate descent algorithm
#
#****** Input ******
#   X       : m-by-n matrix, sparse and non negative
#   W  : m-by-r matrix
#   H  : r-by-n matrix; r is the desired rank of the factorization
#   maxiter  = maximum number of iterations (default: 1000)

#****** Output ******
#   W       : m-by-r low-rank factor W
#   H       : r-by-n low-rank factor H
#  err      : vector containing evolution of relative error along
#            iterations ||X-max(0,WH)||_F / || X ||_F

using LinearAlgebra

function alternatescheme(X,W,H,maxiter)
    
    #Scale W,H
    normX = norm(X)
    Z = max.(0, W*H)
    optscal = sum(Z.*X) / (norm(Z)^2)
    W = optscal*W

    err = [norm(X - (optscal*Z))./normX]
    WH = W*H
    for iter = 1:maxiter
        updateH(X,W,H,WH)
        updateH(X',H',W',WH')
        err = push!(err, norm(X - max.(0, WH))./normX)
        if err[end]<=1e-4
            maxiter = iter
            break
        end
    end
    print("total it :",maxiter)

    #err = err./ normX

    return W,H,err
end 