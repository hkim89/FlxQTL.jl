##### permutation for no prior 
#actual null parameters by genescan from the data (not permuted) to permute 
function scan0loco(cross::Int64,n::Int64,m::Int64,Tg::Array{Float64,3},Λg::Matrix{Float64},Y::Array{Float64,2},XX::Markers
        ;Xnul::Array{Float64,2}=ones(1,n),Z=diagm(ones(m)),itol=1e-3,tol::Float64=1e-4,ρ=0.001)

         # LODs=zeros(p);
        Chr=unique(XX.chr); nChr=length(Chr);NulKc=[];tNuls=[];
        Λc= Array{Array{Float64,1}}(undef,nChr);fill!(Λc,zeros(m))
        Y0= similar(Y)

        #    tbyt, init= transByTrait(m,Tc,λc,Y,Z,Xnul,Prior)
        if (cross!=1)
            X1=mat2array(cross,XX.X)
         else
            X1=similar(XX.X) #pre-assigned 
         end
       
        
      if (Z!=diagm(ones(m)))
         init0=initial(Xnul,Y,Z,false)  
       else #Z0=I
         init0=initial(Xnul,Y,false)
      end
  
     for i= eachindex(Chr)
         maridx=findall(XX.chr.==Chr[i]);
            Λc[i], tt,init = getKc(init0,Y,Tg[:,:,i],Λg[:,i];Xnul=Xnul,m=m,Z=Z,itol=itol,tol=tol,ρ=ρ)
        if (cross!=1)
           @fastmath @inbounds X1[:,:,maridx]=transForm(Tg[:,:,i],X1[:,:,maridx],cross)
         else
           @fastmath @inbounds X1[maridx,:]=transForm(Tg[:,:,i],XX.X[maridx,:],cross)
        end
          Y0[:,:]= trans2iid(tt.Y,1.0,tt.Σ,Λg[:,i],Λc[i]) # preparing transformed Y to be iid for permutation
            
          tNuls=[tNuls;TNul(Y0,tt.Xnul,tt.Z,tt.Σ)];NulKc=[NulKc;init]

     end
        
         return Λc, tNuls, NulKc, X1
 
end

#no loco-null scan for permutation
function Scan0(cross::Int64,n::Int64,m::Int64,Tg,Λg,Y::Array{Float64,2},XX::Markers;Z=diagm(ones(m)),
            Xnul::Array{Float64,2}=ones(1,n),itol=1e-3,tol::Float64=1e-4,ρ=0.001)

    
    # q=size(Z,2);  
    # p=Int(size(XX.X,1)/cross); 

    ## picking up initial values for parameter estimation under the null hypothesis
    if (Z!=diagm(ones(m)))
         init0=initial(Xnul,Y,Z,false)  
       else #Z0=I
         init0=initial(Xnul,Y,false)
      end

   λc, tNul, NulKc = getKc(init0,Y,Tg,Λg;Xnul=Xnul,m=m,Z=Z,itol=itol,tol=tol,ρ=ρ)  
    
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
      
   
    return λc,TNul(Y0,tNul.Xnul,tNul.Z,tNul.Σ), NulKc, X1
end


## finding distribution of max lod's for a multivariate model by permutation for 4waycross/intercross
function locoPermutation(cross::Int64,p::Int64,q::Int64,n::Int64,m::Int64,X::Union{Array{Float64,2},Array{Float64,3}},
        chr::Array{Any,1},tnuls,nulKc,λg::Array{Float64,2},λc::Array{Array{Float64,1},1}
        ;tol0=1e-3,tol::Float64=1e-4,ρ=0.001)

#     n=size(Y,2); p=size(X,1);
    Chr=unique(chr)
    kmin=1; lod=zeros(p);H1par=[]

    rng=shuffle(1:n);
    Y1= Matrix{Float64}(undef,m,n)
    
    for l=eachindex(Chr)
         maridx=findall(chr.==Chr[j]);nmar=length(maridx)
       #transformation back to correlated Y after permutation
         Y1 = trans2iid(tnuls[l].Y[:,rng],1.0,tnuls[l].Σ,λg[:,l],λc[l],true) 
        
         if (cross!=1)
             X1=X[:,:,maridx]
            else 
             X1=X[maridx,:]
           end    
       
        #initial parameter values for permutations are from genome scanning under the null hypothesis.
        perm_est0=nulScan(nulKc[l],kmin,λg[:,l],λc[l],Y1,tnuls[l].Xnul,tnuls[l].Z,tnuls[l].Σ,true;itol=tol0,tol=tol,ρ=ρ)    
        LODs,H1par0=marker1Scan(nmar,q,kmin,cross,perm_est0,λg[:,l],λc[l],Y1,tnuls[l].Xnul,X1,tnuls[l].Z;tol0=tol0,tol1=tol,ρ=ρ)
         
          lod[maridx]= LODs;  H1par=[H1par; H1par0]
           
        end

    return maximum(lod), H1par
end

function permutation(nperm::Int64,cross::Int64,p::Int64,q::Int64,X::Union{Array{Float64,2},Array{Float64,3}},
        tNul::TNul,Nullpar::InitKc,λg::Array{Float64,1},λc::Array{Float64,1}
        ;tol0=1e-3,tol::Float64=1e-4,ρ=0.001)

    # n=size(Y,2);    q=size(Z,2);
    lod=zeros(nperm);H1par=[]; Y2=similar(Y)

    # init=Init(Nullpar.B,Nullpar.τ2,Nullpar.Σ)
    
    for l= 1:nperm
        ### permuting a phenotype matrix by individuals
        permutY!(Y2,tNul.Y,1.0,tNul.Σ,λg,λc)
        #initial parameter values for permutations are from genome scanning under the null hypothesis.
        perm_est0=nulScan(Nullpar,1,λg,λc,Y2,tNul.Xnul,tNul.Z,tNul.Σ;ρ=ρ,itol=tol0,tol=tol)
        LODs,H1par0=marker1Scan(p,q,1,cross,perm_est0,λg,λc,Y2,tNul.Xnul,X,tNul.Z,true;tol0=tol0,tol1=tol,ρ=ρ)
         lod[l]= maximum(LODs);  H1par=[H1par; H1par0]

            if (mod(l,100)==0)
              println("Scan for $(l)th permutation is done.")
            end
        end

    return lod, H1par
end

#MVLMM
function permutation(nperm::Int64,cross::Int64,Y::Array{Float64,2},X::Union{Array{Float64,2},Array{Float64,3}},
        Nullpar::Result,λg::Array{Float64,1},Xnul_t;tol0=1e-3,tol::Float64=1e-4,ρ=0.001)

     m=size(Y,1);
    kmin=1;lod=zeros(nperm);H1par=[]
   
    init=Init0(Nullpar.B,Nullpar.Vc,Nullpar.Σ)
    
     for l= 1:nperm
        ### permuting a phenotype matrix by individuals
        Y2=permutY(Y,Nullpar.Vc,Nullpar.Σ,λg);

        #initial parameter values for permutations are from genome scanning under the null hypothesis.
         perm_est0=nulScan(init,kmin,λg,Y2,Xnul_t;itol=tol0,tol=tol,ρ=ρ)
        LODs,H1par0=marker1Scan(m,kmin,cross,perm_est0,λg,Y2,Xnul_t,X;tol0=tol0,tol1=tol,ρ=ρ)
    
         lod[l]= maximum(LODs);  H1par=[H1par; H1par0]
             if (mod(l,50)==0)
              println("Scan for $(l)th permutation is done.")
             end
        end


    return lod, H1par
end



function permutationTest(Kg::Union{Matrix{Float64},Array{Float64,3}},Y,XX::Markers,nperm::Int64,cross::Int64,LOCO::Bool=false;
       pval=[0.05 0.01],m=size(Y,1),Z=diagm(ones(m)),Xnul=ones(1,size(Y,2)),itol=1e-4,tol0=1e-4,tol=1e-4,ρ=0.001)
       
       p=Int(size(XX.X,1)/cross); n=size(Y,2); q=size(Z,2)
        if (!LOCO)
    #permutation without LOCO
          Tg,λg=K2Eig(Kg)
          λc,tNul, NulKc, X1=Scan0(cross,n,m,Tg,Λg,Y,XX;Z=Z,Xnul=Xnul,itol=itol,tol=tol,ρ=ρ)
          maxLODs, H1par= permutation(nperm,cross,p,q,X1,tNul,NulKc,λg,λc;tol0=tol0,tol=tol,ρ=ρ)
          maxLODs=convert(Array{Float64,1},maxLODs)
            
        else #loco test
            maxLODs = zeros(nperm);H1par=[]
            Tg,Λg=K2Eig(Kg,true)
            Λc, tNuls, NulKc, X1 = scan0loco(cross,n,m,Tg,Λg,Y,XX;Xnul=Xnul,Z=Z,itol=itol,tol=tol,ρ=ρ)
            for l=1:nperm
                maxLODs[l],H1par0 = locoPermutation(cross,p,q,n,m,X1,XX.chr,tNuls,NulKc,Λg,Λc;tol0=tol0,tol=tol,ρ=ρ)
                H1par=[H1par;H1par0]
                if (mod(l,100)==0)
                    println("Scan for $(l)th permutation is done.")
                end
            end
        end
            
    cutoff= quantile(maxLODs,1.0.-pval)

    return maxLODs, H1par, cutoff

end







