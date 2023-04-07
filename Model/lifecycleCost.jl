function lifecycleCost(Q, xL,xU,xts,xa,m, track, gen, stor, ngFT, qsys)
    #===========================================================================
    # This subroutine calculates the the total lifecycle savings by implementing
    # a solar hybridization strategy: fLCS = NG cost avoided - investment cost
    # Written by: M. D. Stuber, May 17, 2018, 
    # Modified by Justin Rastinejad July 2022,
    # Earlier version of this code is used in the paper: Stuber (2018), DOI: 10.3390/pr6070076
    # This code is used in the paper: Rastinejad, Putnam, Stuber (2023)
    # INPUTS
    #   Q: vector of heat outputs from solar array for each time point
    #   X: vector of optimization variable bounds
    #   xts: thermal storage decision variable (h)
    #   xa: aperture area decision variable (m2)
    #   m: capital model to use (=1 for convex, =2 for nonconvex)
    #   track: tracking system =1 for 1-axis, =0 for fixed
    #   gen: solar technology PV or PTC
    #   stor: storage type EES or TES
    #   ngFT: natural gas price in $/ft3
    #   qsys: mean IPH system demand in MW
    # OUTPUTS
    #   flsc: total lifecycle savings if you implement solar
    ===========================================================================#
    
    if gen == "PTC"
        Cts0 = 45.14 #TES cost prefactor
        Ca0  = 425.0 #PTC cost prefactor
        expGen = 0.92 #exponential value for PTC costs
        expES = 0.91 #exponential value for energy costs
    end

    if gen == "PV"
        if track==0
            Ca0  = 200.18 # 200.18 PV cost prefactor
            expGen = 0.9617 # exponential value for the 0-track PV costs
        elseif track==1
            Ca0  = 223.49 #223.49 #PV1 cost prefactor
            expGen = 0.9586 #exponential value for the 1-axis PV costs
        end
        if stor =="EES"
            Cts0 = 736.38 #EES cost prefactor. Originally 736.38 on 12/5/22
            expES = 0.935 #exponential value for EES
        elseif stor == "TES"
            Cts0 =  45.14 #TES cost prefactor
            expES = 0.91 #exponential value for TES
        end
    end
    # economic parameters of the model
    r_th = 0.01 # conventional fuel price inflation
    r_d = 0.10 # discount rate #inflation sorta deal
    r_c = 0.065 # capital APR (loan interest basically)
    term = 10::Int # loan term [years]
    L = 30::Int # total project lifetime [years]
    c_tax = 0 # $ per metric ton of CO2
    emis = 53.07 # kg per mmbtu for natural gas
    ng = (ngFT + c_tax/1000*emis)/293.1/1.037 #9.87 Commercial NG cost $/mmbtu -> $/kWh 
    #taken July 2021 (https://www.eia.gov/dnav/ng/hist/n3020us3m.htm) (if data in $/1000ft^3, get rid of the /1.037 in line above)
    #ng = (4.73 + c_tax/1000*emis)/293.1/1.037 # Industrial NG cost $/mmbtu -> $/kWh taken for July2021 from (https://www.eia.gov/dnav/ng/hist/n3035us3m.htm)

    ##################################################
    t = range(1,stop=8760,length=8760)
    qmean = qsys*1000.0 # mean thermal power demand of IPH [kW]
    qdev = 0.00 # deviation fraction from the mean for ramp up/down (∈[0,1], 0:constant, >0:transient load)
    qp_peak = qmean*(1.0+qdev)# peak IPH demand
    qproc = qmean*(1.0.+qdev*sin.((t.-7.0)/12.0*pi)) # dynamic IPH demand [kW] vector
    qsol = Q*xa/1000.0 # total solar power [kW] vector
    stLoad = qsol-qproc # total storage load [kW] vector
    # initialize arrays
    qstored = zeros(typeof(xts),length(Q)+1) #total thermal energy stored
    qex   = zeros(typeof(xts),length(Q))# excess heat, h(q_s-q_l) in paper
    qloss = copy(qex)# heat loss, q_l in paper
    qanc  = ones(typeof(xts),length(Q))*-qmean# ancillary heat, q_a in paper
    qng = copy(qanc)*-1.0# heat from NG system, q_ng in paper
    dt = 1 # 1h discrete time points, h in paper
    # set up the energy balance equations and calculate SFs
    SFnum=0.0# initialize numerator of solar fraction
    η = 0.85                                      #Li-Ion Battery Effiiency
    ql = copy(qex);                                 #Losses due to EES efficiency
    DoD = 0.8;                                      #Depth of discharge, original value is 0.8
    if stor == "TES"
        for i in 1:length(Q)
            qex[i]=dt*(stLoad[i])
            qstored[i+1]=minD1(qp_peak*xts,maxD1(qstored[i]+qex[i],0.0))
            qanc[i]=qstored[i]+qex[i]-qstored[i+1]
            qloss[i]=maxD1(0.0,qanc[i])
            qng[i] = -minD1(0.0,qanc[i])
            SFnum += qsol[i]-maxD1(0.0,qstored[i]+qex[i]-qp_peak*xts)
        end
        SF = SFnum/sum(qproc)# smooth solar fraction
        Cp = [ng*(1+r_th)^(i-1) for i=1:L]*sum(qproc) # opex w/o solar (gas cost)
        hL = xL[1] # lower bound on storage (h)
        aL = xL[2] # lower bound on area (m^2)
        hU = xU[1] # upper bound on storage (h)
        aU = xU[2] # upper bound on area (m^2)
        Cts = Cts0*((hL*qp_peak)^(expES)+((hU*qp_peak)^(expES)-(hL*qp_peak)^(expES))/(hU-hL)*(xts-hL))
        Ca = Ca0*((aL)^(expGen)+(aU^(expGen)-aL^(expGen))/(aU-aL)*(xa-aL))
        #
        if m == :convex # convex capital model
            Ccap0 = Cts+Ca # (16) in paper
        else # nonconvex capital model
            Ccap0 = Ca0*xa^(expGen) + Cts0*(xts*qp_peak)^(expES) # (14) in paper
        end
        # vectorize annual capital costs
        Ccap = [r_c*Ccap0*(1.0+r_c/12.0)^(12*term)/((1+r_c/12)^(12*term)-1) for i=1:L]
        if term<L Ccap[term+1:L,1] = zeros(L-term,1)  end # set capital cost (debt service) =0 after financing is paid
        disc = [(1+r_d)^(-i) for i = 1:L]# vectorize discount factor for time value
        flcs = (SF*Cp - Ccap)'*disc # (13) in paper
        return flcs
    elseif stor == "EES"
        for i in 1:length(Q)
            qex[i]=dt*(stLoad[i])
            ql[i] = (1.0/η-1.0)*maxD1(minD1(η*qstored[i],-qex[i]),0.0)# no if-statements version
            qstored[i+1] = minD1(qp_peak*xts,maxD1((qstored[i]+qex[i]-ql[i]),0.0))
            qanc[i]=qstored[i]+qex[i]-ql[i]-qstored[i+1]
            qloss[i]= maxD1(0,qanc[i])
            qng[i] = -minD1(0,qanc[i])
        end
        SF = 1.0-sum(qng)/sum(qproc)
        Cp = [ng*(1+r_th)^(i-1) for i=1:L]*sum(qproc) # opex w/o solar (gas cost)
        hL = xL[1] # lower bound on storage (h)
        aL = xL[2] # lower bound on area (m^2)
        hU = xU[1] # upper bound on storage (h)
        aU = xU[2] # upper bound on area (m^2)
        Cts = Cts0*((hL*qp_peak)^(expES)+((hU*qp_peak)^(expES)-(hL*qp_peak)^(expES))/(hU-hL)*(xts/DoD-hL))
        Ca = Ca0*((aL)^(expGen)+(aU^(expGen)-aL^(expGen))/(aU-aL)*(xa-aL))
        #
        if m == :convex # convex capital model
            Ccap0 = Cts+Ca # (16) in paper
        else # nonconvex capital model
            Ccap0 = Ca0*xa^(expGen) + Cts0*(xts/DoD*qp_peak)^(expES) # (14) in paper
        end
        # vectorize annual capital costs
        Ccap = [r_c*Ccap0*(1.0+r_c/12.0)^(12*term)/((1+r_c/12)^(12*term)-1) for i=1:L]
        if term<L Ccap[term+1:L,1] = zeros(L-term,1)  end # set capital cost (debt service) =0 after financing is paid
        disc = [(1+r_d)^(-i) for i = 1:L]# vectorize discount factor for time value
        flcs = (SF*Cp - Ccap)'*disc # (13) in paper
        return flcs
    end
end