using LinearAlgebra
using PlotlyLight
using MAT


function test()
    m = 50
    n = 50
    r = 5
    maxiter = 100

    
        W = randn(m,r)
        H = randn(r,n)

        Wt = randn(m,r)
        Ht = randn(r,n)

        X = max.(0, Wt*Ht)

       # filename = "data500" * string(i) * ".mat"
        #filename = "data501.mat"
       # file = matread(filename)

        #X = file["X"]
        #time1 = file["time"]

        (U,S,V) = svd(X)
        Ur = U[:, 1:r]
        Sr = Diagonal(S[1:r])
        Vr = V[:, 1:r]

        W = Ur
        H = Sr * Vr'

        time = @elapsed begin (W,H,err_julia) = alternatescheme(X,W,H,maxiter) end
        println(" - error = ",err_julia[end]," - time = ",time)
        

        layout = Config(yaxis=Config(type="log"))
        Plot(y=err_julia, layout=layout)

       

    
end

test()