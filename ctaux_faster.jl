module CTAUX_faster
    export ctaux_main
    using Random
    using FFTW
    using Dierckx
    using LinearAlgebra
    using InteractiveUtils

    struct Vertex
        tau::Float64        
        spin::Float64              
    end

    mutable struct Vertices
        num::Int64
        values::Array{Vertex,1}
        indices::Array{Int64,1}
        currentk::Int64
    end

    struct Greenfunction
        spl::Dierckx.Spline1D
    end



    mutable struct QMC_variables
        count::Int64
        vertices::Vertices
        averaged_sign::Float64   
        P::Array{Int64,1} #Distribution of the orders
        N::Array{Float64,3}
        Dr::Array{Float64,2}
        Dc::Array{Float64,2}
        S::Array{Float64,2}
        Scount::Array{Float64,2}
        R_temp::Vector{Float64}
        L_temp::Vector{Float64}
        Mklm_temp::Array{Float64,3}
        ev_temp::Vector{Float64}
        g0l_temp::Matrix{Float64}
        ntemp_temp::Matrix{Float64}
        Mkl_temp::Vector{Float64}


        function QMC_variables(count,vertices,averaged_sign,P,N,Dr,Dc,S,Scount,mkink,norbs)
            R_temp = zeros(mkink)
            L_temp = zeros(mkink)
            Mklm_temp = zeros(mkink,mkink,norbs)
            ev_temp = zeros(Float64,norbs)
            g0l_temp = zeros(Float64,mkink,norbs)
            ntemp_temp = zeros(Float64,mkink,mkink)
            Mkl_temp = zeros(norbs)
            return new(count,vertices,averaged_sign,P,N,Dr,Dc,S,Scount,R_temp,L_temp,Mklm_temp,ev_temp,g0l_temp,ntemp_temp, Mkl_temp)
        end
    end

    struct QMC_parameters
        ntime::Int64
        mfreq::Int64
        mu::Float64
        K::Float64
        U::Float64
        beta::Float64
        mesh_tau::Array{Float64,1}
        mesh_omega::Array{Float64,1}
        L::Int64
        norb::Int64
        expgamma0::Tuple{Float64,Float64}
        delta::Array{ComplexF64,3}
        gtau0::Array{Greenfunction,1}    
        
    end

    struct QMC
        v::QMC_variables
        p::QMC_parameters
    end

    function test()
        beta=10.0
        U = 2.0
        mu=U/2
        K = 1.0
        V = 1.0 #Strength of the hybridization
        ctaux_main(U,K,beta,V,mu)
    end

    function ctaux_main(U,K,beta,V,mu;ntime=1024,mfreq=1024,mkink=1024,L=1,nthermal=1000,norb=2,mqs = 1000000)
        qmc = ctaux_initialize(U,K,beta,V,mu;ntime=ntime,mfreq=mfreq,mkink=mkink,L=L,norb= norb)
        @time ctaux_qmc!(qmc,nthermal,measurements=false)
        @time ctaux_qmc!(qmc,mqs,measurements=true)
    end

    function print_steps(istep,numsteps,averagek,averagesign,number_accept)
        if istep % (numsteps ÷ 10) == 0 
            println("Number of QMC steps: ",istep," of ",numsteps)
            println("Average order k: ",averagek/istep)
            println("Average sign :", averagesign/istep)
            println("Total","\t","accepted","\t","rate")
            println(istep,"\t",number_accept,"\t",number_accept/istep)
        end
    end

    function ctaux_qmc!(q::QMC,numsteps;measurements = true)

        number_accept = 0

        averagesign = 0.0
        averagek = 0.0

        for istep=1:numsteps
            pass,rsign = update!(q)
            averagesign += rsign
            averagek += get_currentk(q)
            number_accept += ifelse(pass,1,0)
            #if pass
            #    number_accept += 1
            #end
            
            

            print_steps(istep,numsteps,averagek,averagesign,number_accept)

            

            if measurements
                calc_S!(q)
            end
            
        end

        if measurements 
            ntime = q.p.ntime
            S = get_S(q)
            S ./= numsteps
            Scount = get_Scount(q)
            Scount .*= ntime / sum(Scount)
    
            for i = 1:ntime
                S[i, :] ./= Scount[i]
            end
        end
        println("Average sign = ",averagesign/numsteps)
        return
    end

    function print_vertices(q::QMC)
        for i=1:get_currentk(q)
           vertex = get_timeordered_vertex(q::QMC,i)
           println(vertex)
        end
        println("indices: ",get_indices(q))
    end

    function update!(q::QMC)
        r = rand()
        #println(r)
        
        if r > 0.5
            pass,rsign = qmc_insert!(q)
        else
#            print_vertices(q)
            pass,rsign = qmc_remove!(q)
            
        end
        #println("pass $pass")
        return pass,rsign
    end

    function qmc_insert!(q::QMC)        
        vertex,position = generate_insert_vertex(q)
        det_ratio = calc_insert_ratio!(q,vertex)
        #println(det_ratio)
        ratio = det_ratio[1]*det_ratio[2]*get_K(q)/(get_currentk(q)+1)

        #println("det_ratioup and down","$det_ratio")
        rsign = sign(ratio)
        pass = MP_test(ratio)
        #println("ratio: $ratio")

        if pass
            insert_vertex!(q,vertex,position)       
            currentk_add1!(q)
            update_insert_N!(q,position,det_ratio)

        end

        return pass,rsign
    end

    function update_insert_N!(q::QMC,position,det_ratio)
        update_insert_N!(q,position,1,det_ratio[1])
        update_insert_N!(q,position,2,det_ratio[2])
        return
    end

    function update_insert_N!(q::QMC,position,rspin,det_ratio)
        k = get_currentk(q)        
        if k == 1
            updateN_direct!(q,1/det_ratio,1,1,rspin)            
            return
        end
        lambda=det_ratio

        #R = zeros(Float64,k-1)
        R = get_R_temp(q,k-1)
        N = view(get_N(q),1:k-1,1:k-1,rspin)
        Dr = view(get_Dr(q),1:k-1,rspin)
        mul!(R,N',Dr)
        R .*= -1/lambda 
        #R = -N'*get_Dr(q)[1:k-1,rspin]/lambda
        #R = - nmat[1:k-1,1:k-1,rspin]'*Dr[1:k-1,rspin]/λ

        #L = zeros(Float64,k-1)
        L = get_L_temp(q,k-1)
        Dc = view(get_Dc(q),1:k-1,rspin)
        mul!(L,N,Dc)
        L .*= -1/lambda
        #L = - N*get_Dc(q)[1:k-1,rspin]/lambda
#        nmat[1:k-1,1:k-1,rspin]*Dc[1:k-1,rspin]/λ

        ntemp = get_ntemp_temp(q,k) #zeros(Float64,k,k)
        #nmatt = zeros(Float64,k-1,k-1)
        #nmatt=N

        for j in 1:k-1
            for i in 1:k-1
                N[i,j] += lambda*L[i]*R[j]
            end
        end
        as = get_indices(q)[position] #indexcon[vertex[3]]        
        ntemp[as,as] = 1/lambda
        #println("as ",as,"\t",k)
        ntemp[1:k-1,1:k-1] .= view(N,1:k-1,1:k-1)
        ntemp[k,1:k-1] .= view(R,1:k-1)
        ntemp[1:k-1,k] .= view(L,1:k-1)        

        updateN_direct!(q,ntemp,rspin)  

#        nmat[1:k,1:k,rspin] = ntemp[1:k,1:k]

        return
    end

    function qmc_remove!(q::QMC)
        if get_currentk(q) == 0
            pass = false
            rsign = 1.0
            return pass,rsign
        end
        vertex,position = generate_remove_vertex(q)
        det_ratio = calc_remove_ratio!(q,vertex,position)
        ratio = (get_currentk(q)/get_K(q))*det_ratio[1]*det_ratio[2]
        rsign = sign(ratio)
        #println("ratio_remove: $ratio")

        pass = MP_test(ratio)
        if pass
            index = get_indices(q)[position]
            remove_vertex!(q,vertex,position)
            currentk_remove1!(q)
            update_remove_N!(q,index,det_ratio)
        end

        return pass,rsign
    end

    function update_remove_N!(q::QMC,position,det_ratio)
        update_remove_N!(q,position,1,det_ratio[1])
        update_remove_N!(q,position,2,det_ratio[2])
        return
    end

    function update_remove_N!(q::QMC,index,rspin,det_ratio)
#        println("position: $position")
        as = index
        k = get_currentk(q)
        if k == 0
            return
        end


        lambda = det_ratio
        ntemp = get_ntemp_temp(q,k) #zeros(Float64,k,k)
        N = get_N(q)
        @inbounds for j in 1:k
            jj = j + ifelse(j - as>=0,1,0)
            for i in 1:k
                ii = i + ifelse(i - as >= 0,1,0)
                ntemp[i,j] = N[ii,jj,rspin] -N[ii,as,rspin]*N[as,jj,rspin]/lambda
            end
        end

        updateN_direct!(q,ntemp,rspin)
        return
    end    

    function calc_remove_ratio!(q,vertex,position)
        index = get_indices(q)[position]
        det_ratio_up = get_N(q)[index,index,1]
        det_ratio_down = get_N(q)[index,index,2]
        return det_ratio_up,det_ratio_down     
    end

    function MP_test(ratio)
        r = rand()
        pass = ifelse(min(1.0,abs(ratio)) > r,true,false)
        return pass
    end

    function calc_insert_ratio!(q::QMC,vertex)
        ispin = 1
        det_ratio_up = calc_insert_ratio_spin!(q,vertex,ispin)
        ispin = 2
        det_ratio_down = calc_insert_ratio_spin!(q,vertex,ispin)
        return det_ratio_up,det_ratio_down
    end

    function calc_insert_ratio_spin!(q::QMC,vertex,ispin)
        k = get_currentk(q)
        ssign = (-1)^(ispin-1)
        #println(get_expgamma0(q,1))
        ev =  ifelse(ssign*vertex.spin==1,get_expgamma0(q,1),get_expgamma0(q,2))

        Gtau0 = calc_Gf0_tau(q,1e-8,ispin)       
        #println("Gtau0 $Gtau0") 
        Dd = ev - Gtau0*(ev-1)  
        #println("Dd ", Dd)

        if k ==0    
            det_ratio = Dd
            return det_ratio
        end

        tau_j = vertex.tau
        sj = vertex.spin
        ev = ifelse(ssign*sj==1,get_expgamma0(q,1),get_expgamma0(q,2))     
        for i in 1:k
            tau_i =get_timeordered_vertex(q,i).tau 
            Gtau0 = calc_Gf0_tau(q,tau_i,tau_j,ispin)
            value = -Gtau0*(ev-1)
            updateDc!(q,value,i,ispin)
        end

        tau_i = vertex.tau
        for j=1:k
            vertex_j = get_timeordered_vertex(q,j)
            tau_j = vertex_j.tau
            s_j = vertex_j.spin
            Gtau0 = calc_Gf0_tau(q,tau_i,tau_j,ispin)
            ev = ifelse(ssign*s_j==1,get_expgamma0(q,1),get_expgamma0(q,2))
            value = -Gtau0*(ev-1)
            updateDr!(q,value,j,ispin)   
        end
        lambda = Dd - dot_Dr_N_Dc(q,k,ispin)

        det_ratio = lambda[1]
    
        return det_ratio
    end

    @inline function updateDc!(q,value,i,ispin)
        index = get_indices(q)[i]
        @inbounds q.v.Dc[index,ispin] = value
        return
    end

    @inline function updateDr!(q,value,i,ispin)
        index = get_indices(q)[i]
        q.v.Dr[index,ispin] = value
        return
    end

    @inline function updateN_direct!(q,values,ispin)
        n = size(values)
        @inbounds for j=1:n[2]
            for i=1:n[1]
                q.v.N[i,j,ispin] = values[i,j]
            end
        end
        #q.v.N[1:n[1],1:n[2],ispin] = values[1:n[1],1:n[2]]
        return
    end

    @inline function updateN_direct!(q,value::Float64,i,j,ispin)
        @inbounds q.v.N[i,j,ispin] = value
        return
    end

    function dot_Dr_N_Dc(q::QMC,k,ispin)
        Dr = get_Dr(q,k,ispin)
        Dc = get_Dc(q,k,ispin)
        N = view(get_N(q),1:k,1:k,ispin)
        temp = get_R_temp(q,k)
        mul!(temp,N,Dc)
        return dot(Dr,temp)

        #return dot(view(q.v.Dr,1:k,ispin),view(q.v.N,1:k,1:k,ispin)*view(q.v.Dc,1:k,ispin))
    end

    function remove!(q::QMC,measurements)
    end

    function ctaux_initialize(U,K,beta,V,mu;ntime=1024,mfreq=1024,mkink=1024,L=1,seednum=1234,norb=2)
        count = 0
        Random.seed!(seednum)
        gamma0 =acosh(1.0+(beta*U)/(2.0*K))
        expgamma0 = (exp(gamma0),exp(-gamma0))
        mesh_tau = calc_linear_mesh(0,beta,ntime)
        mesh_omega = calc_linear_mesh(π/beta,(π/beta)*(2mfreq-1.0),mfreq)

        delta = make_hybridization(mesh_omega,V,norb)
        Ef = -mu
        Gf0_omega = make_Gf0_omega(mesh_omega,norb,Ef,delta,U)

        Gf0_tau = make_Gf0_tau(Gf0_omega,mesh_omega,ntime,beta)
        Gf0_tau = - Gf0_tau

        #println(typeof(Gf0_tau))
        gtau0 = set_greenfunctions(Gf0_tau,mesh_tau) 

        vertices = Vertices(0,Vertex[],Int64[],0)

        #=
    struct QMC_parameters
        ntime::Int64
        mfreq::Int64
        mu::Float64
        K::Float64
        U::Float64
        beta::Float64
        mesh_tau::Array{Float64,1}
        mesh_omega::Array{Float64,1}
        L::Int64
        norb::Int64
        expgamma0::Tuple{Float64,Float64}
        delta::Array{ComplexF64,2}
        gtau0::Array{Greenfunction,1}        
    end
        =#
        p = QMC_parameters(ntime,mfreq,mu,K,U,beta,mesh_tau,mesh_omega,L,norb,expgamma0,delta,gtau0)
        #=
    mutable struct QMC_variables
        count::Int64
        vertices::Vertices   
    end

        =#
        averaged_sign = 0
        P = zeros(Int64,mkink)
        N = zeros(Float64,mkink,mkink,norb)
        Dr = zeros(Float64,mkink,norb)
        Dc = zeros(Float64,mkink,norb)
        #global S,Scount
        S = zeros(Float64,ntime,norb)
        Scount = zeros(Float64,ntime,norb)

        v = QMC_variables(count,vertices,averaged_sign,P,N,Dr,Dc,S,Scount,mkink,norb)
        qmc = QMC(v,p)

        return qmc
    end



    function make_hybridization(mesh_omega,V,norb)
        mfreq = length(mesh_omega)
        delta = zeros(ComplexF64,mfreq,norb,norb)
        for i=1:mfreq
            z = mesh_omega[i]*im
            for j=1:norb
                delta[i,j,j] = 1
            end
            delta[i,:,:] = delta[i,:,:]*(z-im*sqrt(1-z^2))/2
        end
        delta = delta*V^2        
        return delta
    end

    function make_Gf0_omega(mesh_omega,norb,Ef,delta,U)
        mfreq = length(mesh_omega)
        Gf0 = zeros(ComplexF64,mfreq,norb,norb) #Non-perturbative Green's function in the Matsubara space
        for i in 1:mfreq
            for j in 1:norb
                Gf0[i,j,j] = 1/(im*mesh_omega[i] - Ef - delta[i,j,j] - U/2)
            end
        end
        
        return Gf0
    end

    function make_Gf0_tau(Gf0_omega,mesh_omega,ntime,beta)
        return fft_ω2τ(Gf0_omega,mesh_omega,ntime,beta)
    end

    function calc_linear_mesh(xmin,xmax,n)
        x = zeros(Float64,n)
        for i in 1:n
            x[i] = (xmax-xmin)*(i-1)/(n-1)+xmin
        end        
        return x
    end 

    function generate_insert_vertex(q::QMC)
        tau = rand()*get_beta(q)
        spin = (-1)^(rand(1:2)-1)
#        x = rand(1:get_L(q))
        vertex = Vertex(tau,spin)

        index = 1
        currentk = get_currentk(q)
        indices = get_indices(q)
        values = get_values(q)

        if currentk > 0
            if tau < values[indices[1]].tau
                index = 1
            elseif  tau > values[indices[currentk]].tau #τ > τcon[indexcon[currentk]]
                index = currentk + 1
            else
                i = 1
                for j in 1:currentk
                    i += ifelse(values[indices[j]].tau < tau,1,0)
                end
                index = i
            end
        end     
        
        return vertex,index
    end

    function generate_remove_vertex(q::QMC)
        index = rand(1:get_currentk(q))
        vertex = get_values(q)[get_indices(q)[index]]
        return vertex,index
    end


    function insert_vertex!(q::QMC,vertex,position)
        insert_vertex!(q.v.vertices,vertex,position)
        return
    end

    function insert_vertex!(v::Vertices,tau,spin,position)
        insert_vertex!(v,Vertex(tau,spin),position)

        return
    end

    function insert_vertex!(v::Vertices,vertex,position)
        if v.currentk == 0          
            index = 1  
            insert!(v.values,1,vertex)
            insert!(v.indices,1,index)
            return
        end
        index = v.currentk+1
        insert!(v.indices,position,index)       #insert index at position 
        insert!(v.values,index,vertex)          #insert vertex at index
        return
    end

    function remove_vertex!(q::QMC,vertex,position)
        remove_vertex!(q.v.vertices,vertex,position)
        return
    end

    function remove_vertex!(v::Vertices,vertex,position)
        is = v.indices[position]
        deleteat!(v.values,is)
        deleteat!(v.indices,position)
        for i in 1:v.currentk-1
            v.indices[i] += ifelse(v.indices[i] >= is,-1,0)
        end

        return
    end

    function remove_vertex!(v::Vertices,tau,spin,x,position)
        remove_vertex!(v,Vertex(tau,spin,x),position)
        return
    end



    function fft_ω2τ(f_omega,mesh_omega,ntime,beta)
        mfreq = size(f_omega)[1]
        nl = size(f_omega)[2]
#        println("nl = $nl")
        f_tau = zeros(Float64,ntime,nl,nl)
        for i=1:nl
            for j=1:nl
                f_tau[:,i,j] = fft_backward(ntime,beta,f_omega[:,i,j],mesh_omega)
            end
        end

        for i=1:nl
            for itau=1:ntime
                if f_tau[itau,i,i] > 0
                    f_tau[itau,i,i] = 1e-6 #to avoid positive values
                end
            end
        end
        return f_tau
    end

    function fft_backward(ntime,beta,v_omega,mesh_omega)
        mfreq = length(mesh_omega)
        tail = calc_tails(v_omega,mesh_omega)
        gk = zeros(ComplexF64,2mfreq)
    
        for j=1:mfreq
            gk[2j] = v_omega[j] - tail/(im*mesh_omega[j])
        end
               
        fft!(gk)
  
        
        g_tau = zeros(Float64,ntime)
        g_tau[1:ntime]= real(gk[1:ntime])*(2/beta).-tail/2
     
        a = real(v_omega[mfreq])*mesh_omega[mfreq]/π
        g_tau[1] += a
        g_tau[ntime] += -a

        return g_tau
        
    end  

    function calc_tails(vω,mesh_omega) #to calculate a tail
        mfreq = length(mesh_omega)
        ntail = 128
    
        Sn = 0.0
        Sx = 0.0
        Sy = 0.0
    
        Sxx = 0.0
        Sxy = 0.0
    
        for j in mfreq-ntail:mfreq
            ωn = mesh_omega[j]
            Sn += 1
            Sx += 1/ωn^2
            Sy += imag(vω[j])*ωn
            Sxx += 1/ωn^4
            Sxy += imag(vω[j])*ωn/ωn^2
        end
            rtail = (Sx*Sxy-Sxx*Sy)/(Sn*Sxx - Sx*Sx)
    
        return rtail
    end


    @inline function calc_Gf0_tau(gtau0::Array{Greenfunction,1},tau,l)
        return gtau0[l].spl(tau)
    end

    @inline function calc_Gf0_tau(q::QMC,tau,l)
        return calc_Gf0_tau(q.p.gtau0,tau,l)
    end

    function calc_Gf0_tau(q::QMC,tau_i,tau_j,l)
        dtau = tau_i - tau_j
        if dtau < 0
            dtau += get_beta(q)
            Gtau0 = - calc_Gf0_tau(q,dtau,l) 
        elseif dtau == 0.0
            Gtau0 = calc_Gf0_tau(q,1e-8,l)
        else
            Gtau0 = calc_Gf0_tau(q,dtau,l)        
        end
        return Gtau0
    end

    @inline function get_tau(q::QMC,itau::Int)
        return @inbounds q.p.mesh_tau[itau]
    end

    @inline function get_timeordered_vertex(q::QMC,i)
        return @inbounds get_values(q)[get_indices(q)[i]]
    end

    function get_tau_from_vertices(q::QMC,i)
        return  get_timeordered_vertex(q,i).tau
    end

    function get_spin_from_vertices(q::QMC,i)
        return  get_timeordered_vertex(q,i).tau
    end

    @inline function get_omega(q::QMC,iomega::Int)
        return @inbounds q.p.mesh_omega[iomega]
    end

    @inline function get_expgamma0(q::QMC,upordown)
        i = ifelse(upordown == 1,1,2)
        return @inbounds q.p.expgamma0[i]
    end

    @inline function get_U(q::QMC)
        return q.p.U
    end

    @inline function get_K(q::QMC)
        return q.p.K
    end    

    @inline function get_L(q::QMC)
        return q.p.L
    end

    @inline function get_norb(q::QMC)
        return q.p.norb
    end

    @inline function get_beta(q::QMC)
        return q.p.beta
    end

    @inline function get_ntime(q::QMC)
        return q.p.ntime
    end

    @inline function get_mfreq(q::QMC)
        return q.p.mfreq
    end

    @inline function get_count(q::QMC)
        return q.v.count
    end

    @inline function nextcount!(q::QMC)
        q.v.count += 1
        return
    end

    @inline function get_N(q::QMC)
        return q.v.N
    end

    @inline function get_S(q::QMC)
        return q.v.S
    end

    @inline function get_Scount(q::QMC)
        return q.v.Scount
    end


    @inline function get_Dr(q::QMC)
        return q.v.Dr
    end

    @inline function get_Dr(q::QMC,k,iorb)
        return view(q.v.Dr,1:k,iorb)
    end

    @inline function get_Dc(q::QMC)
        return q.v.Dc
    end

    @inline function get_Dc(q::QMC,k,iorb)
        return view(q.v.Dc,1:k,iorb)
    end



    @inline function get_vertices(q::QMC)
        return q.v.vertices
    end

    @inline function get_currentk(q::QMC)
        return get_vertices(q).currentk
    end


    @inline function currentk_add1!(q::QMC)
        get_vertices(q).currentk += 1
        return
    end    

    @inline function currentk_remove1!(q::QMC)
        get_vertices(q).currentk += -1
        return
    end    

    @inline function get_indices(q::QMC)
        return get_vertices(q).indices
    end

    @inline function get_values(q::QMC)
        return get_vertices(q).values
    end

    @inline function get_expgamma0(q::QMC,i)
        return @inbounds q.p.expgamma0[i]
    end

    @inline function get_L_temp(q::QMC,k)
        return view(q.v.L_temp,1:k)
    end

    @inline function get_R_temp(q::QMC,k)
        return view(q.v.R_temp,1:k)
    end

    @inline function get_Mklm_temp(q::QMC,k)
        return view(q.v.Mklm_temp,1:k,1:k,:)
    end

    @inline function get_ev_temp(q::QMC)
        return q.v.ev_temp
    end

    @inline function get_g0l_temp(q::QMC,k)
        return view(q.v.g0l_temp,1:k,:)
    end

    @inline function get_ntemp_temp(q::QMC,k)
        return view(q.v.ntemp_temp,1:k,1:k)
    end

    @inline function get_Mkl_temp(q::QMC)
        return q.v.Mkl_temp
    end

    

    function set_greenfunctions(Gf0_tau,mesh_tau)        
        nl = size(Gf0_tau)[2]
        println(nl)
        gtau0 = Greenfunction[]        
#        gtau0 = Array{Greenfunction}{undef,nl}
        for l=1:nl
            spl = Dierckx.Spline1D(mesh_tau, Gf0_tau[:,l,l])
            push!(gtau0,Greenfunction(spl))
        end
        return gtau0

    end


    function calc_S!(q::QMC)
        currentk = get_currentk(q)
        norbs = get_norb(q)
        Mklm = get_Mklm_temp(q,currentk)
        #Mklm = zeros(Float64,currentk,currentk,norbs)
        ev = get_ev_temp(q) #zeros(Float64,norbs)
        g0l = get_g0l_temp(q,currentk)#zeros(Float64,currentk,norbs)
        β = get_beta(q)
        ntime = get_ntime(q)
        S = get_S(q)
        Scount = get_Scount(q)
        Mkl = get_Mkl_temp(q)
        #global γ

        for k in 1:currentk
            kk = get_indices(q)[k]
            vertex = get_timeordered_vertex(q,k)
#            kk = indexcon[k]
            τk = vertex.tau #τcon[kk]
            sk = vertex.spin #spincon[kk]
            for sigma in 1:norbs
                ssign = (-1)^(sigma-1)
                ev[sigma] =  ifelse(ssign*sk==1,get_expgamma0(q,1),get_expgamma0(q,2))
                #ev[sigma] = exp(γ*ssign*sk)
            end
            for l in 1:currentk                
                ll = get_indices(q)[l] #indexcon[l]
                τl = get_vertices(q).values[ll].tau #τcon[ll]
                for sigma in 1:norbs
                    Mklm[kk,ll,sigma] = (ev[sigma]-1)*get_N(q)[kk,ll,sigma]
                end
            end
        end

        for l in 1:currentk
            ll = get_indices(q)[l]
#            ll = indexcon[l]
            τl = get_vertices(q).values[ll].tau
#            τl = τcon[ll]
            for sigma in 1:norbs
                g0l[ll,sigma] = calc_Gf0_tau(q,τl,sigma) #Gτ0spl[sigma](τl)
            end
        end

#        global ntime,β
        #Mkl = zeros(Float64,norbs)
        Mkl .= 0.0
        ξ = β/ntime
        id = 3
        for k in 1:currentk
            kk = get_indices(q)[k]
            #kk = indexcon[k]
            τk = get_vertices(q).values[kk].tau
            #τk = τcon[kk]
            iτ = ceil(Int64,ntime*τk/β)
            Mkl[:] .= 0.0#zeros(Float64,norbs)
            for l in 1:currentk
                ll = get_indices(q)[l]
#                ll = indexcon[l]
                for sigma in 1:norbs                
                    Mkl[sigma] += Mklm[kk,ll,sigma]*g0l[ll,sigma]
                end
            end
            
            for j in -id:id
                jj = iτ+j
                if 1 <= jj <= ntime
                    τd = get_tau(q,jj)
                    #τd = τmesh[jj]
                    f = gauss(ξ,τd,τk)
                    for sigma in 1:norbs
                        S[jj,sigma] += f*Mkl[sigma]
                    end
                    Scount[jj] += f
                end
            end

        end
        

    end

    


    function gauss(ξ, x, x0)
        f = (1 / (sqrt(2π) * ξ)) * exp(-(x - x0)^2 / (2 * ξ^2))
        return f
    end

    
end

