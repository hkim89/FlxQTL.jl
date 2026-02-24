
## marker1Scan : CPU 1D-genome scanning under H1 only (with/without loco)
function marker1Scan(nmar,q,kmin,cross,Nullpar::Approx,λg,λc,Y1,Xnul_t,X1,Z1,ν₀,Ψ;ρ=0.001,tol0=1e-3,tol1=1e-4,nchr=0)

        # nmar=size(X1,1);
    if (cross==1) ## scanning genotypes
        B0=hcat(Nullpar.B,zeros(Float64,q));

        lod=@distributed (vcat) for j=1:nmar
            XX=vcat(Xnul_t,@view X1[[j],:])
        B0,τ2,Σ,loglik0 =ecmLMM(Y1,XX,Z1,B0,Nullpar.τ2,Nullpar.Σ,λg,λc,ν₀,Ψ;tol=tol0)
            lod0= (loglik0-Nullpar.loglik)/log(10)
        est1=ecmNestrvAG(lod0,kmin,Y1,XX,Z1,B0,τ2,Σ,λg,λc,ν₀,Ψ;ρ=ρ,tol=tol1,numChr=nchr,nuMarker=j)
            [(est1.loglik-Nullpar.loglik)/log(10) est1]
                 end

    else # cross>1
        ## scanning genotype probabilities

            B0=hcat(Nullpar.B,zeros(Float64,q,cross-1))

          lod=@distributed (vcat) for j=1:nmar
                XX= vcat(Xnul_t, @view X1[2:end,:,j])
                B0,τ2,Σ,loglik0 =ecmLMM(Y1,XX,Z1,B0,Nullpar.τ2,Nullpar.Σ,λg,λc,ν₀,Ψ;tol=tol0)
                  lod0= (loglik0-Nullpar.loglik)/log(10)
                est1=ecmNestrvAG(lod0,kmin,Y1,XX,Z1,B0,τ2,Σ,λg,λc,ν₀,Ψ;ρ=ρ,tol=tol1,numChr=nchr,nuMarker=j)
            [(est1.loglik-Nullpar.loglik)/log(10) est1]
                                  end

     end

    return lod[:,1],lod[:,2]
end

#Z=I
function marker1Scan(nmar,m,kmin,cross,Nullpar::Approx,λg,λc,Y1,Xnul_t,X1,ν₀,Ψ;ρ=0.001,tol0=1e-3,tol1=1e-4,nchr=0)

        # nmar=size(X1,1);
    if (cross==1) ## scanning genotypes
        B0=hcat(Nullpar.B,zeros(Float64,m));
#      f= open(homedir()*"/GIT/fmulti-lmm/result/test_ecmlmm.txt","w")
        lod=@distributed (vcat) for j=1:nmar
            XX=vcat(Xnul_t,@view X1[[j],:])
        B0,τ2,Σ,loglik0 =ecmLMM(Y1,XX,B0,Nullpar.τ2,Nullpar.Σ,λg,λc,ν₀,Ψ;tol=tol0)
                lod0= (loglik0-Nullpar.loglik)/log(10)
        est1=ecmNestrvAG(lod0,kmin,Y1,XX,B0,τ2,Σ,λg,λc,ν₀,Ψ;ρ=ρ,tol=tol1,numChr=nchr,nuMarker=j)
            [(est1.loglik-Nullpar.loglik)/log(10) est1]
#             f=open(homedir()*"/GIT/fmulti-lmm/result/test_ecmlmm.txt","a")
#               writedlm(f,[loglik0 est1.loglik Nullpar.loglik])
#             close(f)
                 end

    else # cross>1
        ## scanning genotype probabilities

        #initialize B under the alternative hypothesis
        B0=hcat(Nullpar.B,zeros(Float64,m,cross-1))

          lod=@distributed (vcat) for j=1:nmar
                XX=vcat(Xnul_t, @view X1[2:end,:,j])
                B0,τ2,Σ,loglik0 =ecmLMM(Y1,XX,B0,Nullpar.τ2,Nullpar.Σ,λg,λc,ν₀,Ψ;tol=tol0)
                 lod0= (loglik0-Nullpar.loglik)/log(10)
                est1=ecmNestrvAG(lod0,kmin,Y1,XX,B0,τ2,Σ,λg,λc,ν₀,Ψ;ρ=ρ,tol=tol1,numChr=nchr,nuMarker=j)
            [(est1.loglik-Nullpar.loglik)/log(10) est1]
                                  end

     end

    return lod[:,1],lod[:,2]
end



##MVLMM
function marker1Scan(nmar,m,kmin,cross,Nullpar::Result,λg,Y1,Xnul_t,X1,ν₀,Ψ;ρ=0.001,tol0=1e-3,tol1=1e-4)

        # nmar=size(X1,1);
    if (cross==1)
        B0=hcat(Nullpar.B,zeros(m))

             lod=@distributed (vcat) for j=1:nmar
               XX= vcat(Xnul_t,@view X1[[j],:])
               B0,Vc,Σ,loglik0 = ecmLMM(Y1,XX,B0,Nullpar.Vc,Nullpar.Σ,λg,ν₀,Ψ;tol=tol0)
                     lod0= (loglik0-Nullpar.loglik)/log(10)
               est1=ecmNestrvAG(lod0,kmin,Y1,XX,B0,Vc,Σ,λg,ν₀,Ψ;ρ=ρ,tol=tol1)
               [(est1.loglik-Nullpar.loglik)/log(10) est1]
                           end

    else #cross>1

        B0=hcat(Nullpar.B,zeros(m,cross-1))

        lod=@distributed (vcat) for j=1:nmar
                XX= vcat(Xnul_t,@view X1[2:end,:,j])
            B0,Vc,Σ,loglik0 = ecmLMM(Y1,XX,B0,Nullpar.Vc,Nullpar.Σ,λg,ν₀,Ψ;tol=tol0)
                 lod0= (loglik0-Nullpar.loglik)/log(10)
                est1=ecmNestrvAG(lod0,kmin,Y1,XX,B0,Vc,Σ,λg,ν₀,Ψ;tol=tol1,ρ=ρ)
                     [(est1.loglik-Nullpar.loglik)/log(10) est1]
                          end

    end
    return lod[:,1], lod[:,2]

end

###### estimating Kc with prior

struct TNull
Y::Matrix{Float64}
Xnul::Matrix{Float64}
Z::Matrix{Float64}
Σ::Matrix{Float64}
Ψ::Matrix{Float64}
end

function transByTrait(m,Tc,Y,Z,Prior,init::Result)

#  if (λc!= ones(m))
        if (Prior!= diagm(ones(m)) && Z!= diagm(ones(m)) )
            Z1, Σ1, Ψ =transForm(Tc,Z,init.Σ,Prior)
        elseif (Z != diagm(ones(m)))
            Z1,Σ1 =  transForm(Tc,Z,init.Σ,true)
            Ψ =Prior
        else #Z & prior =I 
            Σ1 =  transForm(Tc,init.Σ,Z)
            Z1=Z; Ψ = Prior
        end
        
        Y1= transForm(Tc,Y,init.Σ) # transform Y only by row (Tc)

    return Y1,Z1,Σ1,Ψ
end
#estimate Kc with prior
function nul1Scan(init::Init0,kmin,λg,Y,Xnul,Z,m,ν₀,Ψ;ρ=0.001,itol=1e-3,tol=1e-4)
       
      # n=size(Y,2); 

    if (Z!=diagm(ones(m)))   
        B0,Kc_0,Σ1,_=ecmLMM(Y,Xnul,Z,init.B,init.Vc,init.Σ,λg,ν₀,Ψ;tol=itol)
        nulpar=NestrvAG(kmin,Y,Xnul,Z,B0,Kc_0,Σ1,λg,ν₀,Ψ;tol=tol,ρ=ρ)
        
       else #Z=I
        nulpar = nulScan(init,kmin,λg,Y,Xnul,ν₀,Ψ;ρ=ρ,itol=itol,tol=tol)
     end
    return nulpar #Result
end

#H0 MVLMM for Kc estimation with prior
function getKc(Y::Array{Float64,2},Tg::Matrix{Float64},λg::Vector{Float64},init::Init0;m=size(Y,1),Z=diagm(ones(m)), df_prior=m+1,
     Prior::Matrix{Float64}=cov(Y,dims=2)*3,
     Xnul::Array{Float64,2}=ones(1,size(Y,2)),itol=1e-2,tol::Float64=1e-3,ρ=0.001)
     
     Y1,Xnul_t = transForm(Tg,Y,Xnul,1) #null model transformation
     
 
     est0= nul1Scan(init,1,λg,Y1,Xnul_t,Z,m,df_prior,Prior;ρ=ρ,itol=itol,tol=tol)
      Tc, λc = K2eig(est0.Vc)
     
      #trait-wise transformation
      Y2,Z1,Σ1,Ψ = transByTrait(m,Tc,Y1,Z,Prior,est0)
      τ² = 1.0
    
    return λc, TNull(Y2,Xnul_t,Z1,Σ1,Ψ),InitKc(est0.Vc,est0.B,est0.Σ,τ²,est0.loglik)
 
 end
 
 

include("geneScan1.jl")

