using LinearAlgebra
using PlotlyLight
using MAT
# using JLD2 #@save ----.jld2 @load ---.jld2

function test()
    m = 50
    n = 50
    r = 5
    maxiter = 100

    #for i = 1:5
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
        #@profview alternatescheme(X,W,H,maxiter)

        layout = Config(yaxis=Config(type="log"))
        Plot(y=err_julia, layout=layout)

        #Checking error from MATLAB and Julia are same or not

        #mat_data = matread("err1_sparse.mat")
        #error_matlab = mat_data["err"]'


        #tolerance = 1e-5

        # Compare arrays element-wise
        #result = isapprox(error_matlab, err_julia; atol=tolerance)

        # Display the result
        #if result
        #   println("Arrays are approximately equal within 5 digits of accuracy.")
        #else
        #   println("Arrays are not approximately equal within 5 digits of accuracy.")
        #end

    #end
end

test()