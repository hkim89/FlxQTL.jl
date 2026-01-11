

###Permutation test : all permutations are implemented without loco
function trans2iid(Y::Array{Float64,2},τ2_nul::Float64,Σ_nul::Array{Float64,2},λg::Array{Float64,1},λc::Array{Float64,1}
        ,trnsback::Bool=false)

    Y_t=similar(Y);
   # Std=Array{Float64}(undef,m,m)
    Λc=Diagonal(τ2_nul*λc)

    @fastmath @inbounds @views for j in eachindex(λg)
        Std=sqrt(λg[j]*Λc+Σ_nul)
        if (trnsback)
            Y_t[:,j]=Std*Y[:,j]
        else
            Y_t[:,j]=Std\Y[:,j]
        end
    end
    return Y_t
end
##MVLMM
function trans2iid(Y::Array{Float64,2},Vc_nul::Array{Float64,2},Σ_nul::Array{Float64,2},λg::Array{Float64,1},trnsback::Bool=false)

    Y_t=similar(Y);
   # Std=Array{Float64}(undef,m,m)
   @fastmath @inbounds @views for j in eachindex(λg)
         Std=sqrt(Vc_nul*λg[j]+Σ_nul)
        if (trnsback)
            Y_t[:,j]=Std*Y[:,j]
        else
            Y_t[:,j]=Std\Y[:,j]
        end
    end
    return Y_t
end

## permuteY : shuffling a Y matrix columnwise (by individuals)
function permutY!(Y1::Array{Float64},Y::Array{Float64,2},τ2_nul::Float64,Σ_nul::Array{Float64,2}
        ,λg::Array{Float64,1},λc::Array{Float64,1})
#    n=size(Y,2);
    ### permutation by individuals
    #transforming to iid Y1
    # Y_t=trans2iid(Y,τ2_nul,Σ_nul,λg,λc);
    rng=shuffle(axes(Y,2));
    # transforming back to correlated Y1 after permuting
    Y1[:,:]=trans2iid(Y[:,rng],τ2_nul,Σ_nul,λg,λc,true);

end

#MVLMM
function permutY(Y::Array{Float64,2},Vc_nul::Array{Float64,2},Σ_nul::Array{Float64,2},λg::Array{Float64,1})
    n=size(Y,2);
    ### permutation by individuals
    #transforming to iid Y1
    Y_t=trans2iid(Y,Vc_nul,Σ_nul,λg);
    rng=shuffle(1:n);
    # transforming back to correlated Y1 after permuting
    Y2=trans2iid(Y_t[:,rng],Vc_nul,Σ_nul,λg,true);
    return Y2
end

#MVLMM
function permutY!(Y2::Matrix{Float64},Y::Array{Float64,2},Vc_nul::Array{Float64,2},Σ_nul::Array{Float64,2},λg::Array{Float64,1})
    #  n=size(Y,2);
    ### permutation by individuals
    #transforming to iid Y1
    # Y_t=trans2iid(Y,Vc_nul,Σ_nul,λg);
    rng=shuffle(axes(Y,2));
    # transforming back to correlated Y1 after permuting
    Y2[:,:]=trans2iid(Y[:,rng],Vc_nul,Σ_nul,λg,true);

end


## finding distribution of max lod's for a multivariate model by permutation for 4waycross/intercross
function permutation(nperm::Int64,cross::Int64,p::Int64,q::Int64,Y::Array{Float64,2},X::Union{Array{Float64,2},Array{Float64,3}},
        tNul::TNull,Nullpar::InitKc,λg::Array{Float64,1},λc::Array{Float64,1},ν₀
        ;tol0=1e-3,tol::Float64=1e-4,ρ=0.001)

    # n=size(Y,2);    q=size(Z,2);
    lod=zeros(nperm);H1par=[]; Y2=similar(Y)

    init=Init(Nullpar.B,Nullpar.τ2,Nullpar.Σ)
    
    for l= 1:nperm
        ### permuting a phenotype matrix by individuals
        permutY!(Y2,Y,1.0,tNul.Σ,λg,λc)
        #initial parameter values for permutations are from genome scanning under the null hypothesis.
        perm_est0=nulScan(init,1,λg,λc,Y2,tNul.Xnul,tNul.Z,tNul.Σ,ν₀,tNul.Ψ;ρ=ρ,itol=tol0,tol=tol)
        LODs,H1par0=marker1Scan(p,q,1,cross,perm_est0,λg,λc,Y2,tNul.Xnul,X,tNul.Z,ν₀,tNul.Ψ;tol0=tol0,tol1=tol,ρ=ρ)
         lod[l]= maximum(LODs);  H1par=[H1par; H1par0]

            if (mod(l,100)==0)
              println("Scan for $(l)th permutation is done.")
            end
        end

    return lod, H1par
end


#MVLMM
function permutation(nperm::Int64,cross::Int64,p::Int64,Y::Array{Float64,2},X::Union{Array{Float64,2},Array{Float64,3}},
        Nullpar::Result,λg::Array{Float64,1},Xnul_t,ν₀,Ψ;tol0=1e-3,tol::Float64=1e-4,ρ=0.001)

     m=size(Y,1);
    lod=zeros(nperm);H1par=[]
   
    init=Init0(Nullpar.B,Nullpar.Vc,Nullpar.Σ)
    
     for l= 1:nperm
        ### permuting a phenotype matrix by individuals
        Y2=permutY(Y,Nullpar.Vc,Nullpar.Σ,λg);

        #initial parameter values for permutations are from genome scanning under the null hypothesis.
         perm_est0=nulScan(init,1,λg,Y2,Xnul_t,ν₀,Ψ;itol=tol0,tol=tol,ρ=ρ)
        LODs,H1par0=marker1Scan(p,m,1,cross,perm_est0,λg,Y2,Xnul_t,X,ν₀,Ψ;tol0=tol0,tol1=tol,ρ=ρ)
    
         lod[l]= maximum(LODs);  H1par=[H1par; H1par0]
             if (mod(l,100)==0)
              println("Scan for $(l)th permutation is done.")
             end
        end


    return lod, H1par
end



#no loco-null scan for permutation
function Scan0(cross::Int64,Tg,Λg,Y::Array{Float64,2},XX::Markers,Z::Array{Float64,2};
            Xnul::Array{Float64,2}=ones(1,size(Y,2)),m=size(Y,1),df_prior=m+1,
            Prior::Matrix{Float64}=cov(Y,dims=2)*3,itol=1e-3,tol::Float64=1e-4,ρ=0.001)

    
    # q=size(Z,2);  
    # p=Int(size(XX.X,1)/cross); 

    ## picking up initial values for parameter estimation under the null hypothesis
    if (Z!=diagm(ones(m)))
         init0=initial(Xnul,Y,Z,false)  
       else #Z0=I
         init0=initial(Xnul,Y,false)
      end

   λc, tNul, NulKc = getKc(Y,Tg,Λg,init0;Xnul=Xnul,m=m,Z=Z,df_prior=df_prior,Prior=Prior,itol=itol,tol=tol,ρ=ρ)  
    
     if (cross!=1)
        X0=mat2array(cross,XX.X)
     end
#        
             if (cross!=1)
               X1=transForm(Tg,X0,cross)
               else
               X1=transForm(Tg,XX.X,cross)
             end
             
        Y0 = trans2iid(tNul.Y,1.0,tNul.Σ,Λg,λc)
        
   
    return λc, tNul, NulKc, Y0, X1
end

##MVLMM
function Scan0(cross::Int64,Tg,Λg,Y::Array{Float64,2},XX::Markers;Xnul::Array{Float64,2}=ones(1,size(Y,2)),
    m=size(Y,1), df_prior=m+1,Prior::Matrix{Float64}=cov(Y,dims=2)*3,itol=1e-3,tol::Float64=1e-4,ρ=0.001)

   
    # p=Int(size(XX.X,1)/cross);

    #check the prior
    if (!isposdef(Prior))
        println("Error! Plug in a postivie definite Prior!")
     end

     #initialization
       init=initial(Xnul,Y,false)
        if (cross!=1)
            X0=mat2array(cross,XX.X)
         end
    
#             Xnul_t=Xnul*Tg';
             Xnul_t= BLAS.gemm('N','T',Xnul,Tg)
                if (cross!=1)
                   Y,X=transForm(Tg,Y,X0,cross)
                   else
                   Y,X=transForm(Tg,Y,XX.X,cross)
                 end

                  est0=nulScan(init,1,Λg,Y,Xnul_t,df_prior,Prior;itol=itol,tol=tol,ρ=ρ)
      
        return est0,Xnul_t,Y,X

end


include("PpermutationTest1.jl")
include("permutationTest.jl")


