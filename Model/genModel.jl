function genModel(DNI, DHI, Temp, tz, lat,long, track, gen)

    # Made by Justin Rastinejad 7/17/21
    # PV model based on https://www.scirp.org/journal/paperinformation.aspx?paperid=83262
    # PV modeling based on SunPower's E19-320 PV module
    
    #initializing vectors
    if  gen == "PV"
        display("photovoltaic Power/GHI")
    

        t = range(1,stop=length(DNI),length=length(DNI)) # time vector
        GHI = zeros(length(DNI)) #Global Horizatal Irradiation
        T_c = zeros(length(DNI)) # operating Temperature
        inc = ones(length(DNI)) # initialize incidence angle vector
        C = zeros(length(GHI)) #constant
        E_g = zeros(length(GHI)) #bandgap
        R_s = zeros(length(GHI)) #series resistance
        I_0 = zeros(length(GHI)) 
        X1 = zeros(length(GHI)) #breaking up calc pt1
        X2 = zeros(length(GHI)) #breaking up calc pt2
        X3 = zeros(length(GHI)) #breaking up calc pt3
        I_max = zeros(length(GHI)) #current at max power
        V_max = zeros(length(GHI)) #voltage at max power
        Ps = zeros(length(GHI)) #power output in Watts
    
        # calculate the incidence angles for the desired location
        for i = 1:length(t) # loop over all times
            inc[i] = solarAngles(t[i],tz,lat,long, track) # get incidence angles at this time
            #corrected error on 1/12/22 where it wasn't properly using SolarAngles
            if DNI[i]<0 inc[i]=0 end # reset to no-sun condition if dni data is negative
            GHI[i] = DHI[i] + DNI[i] *cos(inc[i]/180*pi)
            T_c[i] = Temp[i] + 273.15 #Operating Temperature at time (K)
            if GHI[i]<0 GHI[i]=0 end #reset to zero if angles make it funky
        end

        #determining I_ph
        alpha= 0.035 #temperature coefficient for closed circuit current (given as I_sc*mu_I_sc)
        I_sc = 6.24 #short circuit current (A)
        T_r = 298.15 #reference temperature (K)
        G_r = 1000 #referece irradiation (W/m^2)
        I_phs = zeros(length(GHI))
        for i = 1:length(I_phs)
            I_phs[i] = (I_sc + alpha*(T_c[i] - T_r))*(GHI[i]/G_r) #light current (A)
        end

        #determining band gap
        alph = 4.730 * 10^(-4) #band gap parameter (e*V/K^2) Table 2 MatLab paper
        beta = 636 #band gap parameter (K) Table 2  MatLab paper
        E_g0 = 1.17 #band gap at T=0K (eV) Table 2  MatLab paper
    
        #determining I_0
        V_oc = 64.8 #open circuit voltage (V)
        q = 1.60218*10^(-19) #Electron charge (Coulumb)
        k = 1.38066*10^(-23) #Boltzmann's constant (J/K)
        N = 96 #number of cells in series 72
        V_mp = 54.7 #maximum power voltage (V)
        I_mp = 5.86 #maximum power current (A)
        n = (q*(2*V_mp - V_oc))/ (N*k*T_r*(I_mp/(I_sc - I_mp)+log((I_sc-I_mp)/I_sc))) #ideality factor (matches paper woooo)
    
        #Calculating power
        for i = 1:length(I_phs)
            C[i] = n*N*k*T_c[i]/q
            E_g[i] = E_g0 - alph*T_c[i]^2/(T_c[i]+beta) #MatLab paper, but same answer as this paper
            R_s[i] = (N*n*k*T_c[i]/q * log(1-(I_mp/I_sc)) + V_oc - V_mp) / I_mp #series resistance (ohm) (matches paper at 25 C)
            I_0[i] = (I_sc)/(exp(q*V_oc/(n*k*N*T_c[i]))-1) * (T_c[i]/T_r)^3 * exp((q*E_g[i]*(1/T_r - 1/T_c[i]))/(n*k)) #saturation current of diode (A) (matches paper at 25 C)
            X1[i] = (C[i] + 2*R_s[i]*(I_0[i] + I_phs[i]))
            X2[i] = -1* (I_0[i] + I_phs[i]) * (C[i]*log(1 + I_phs[i]/I_0[i])+ 2*C[i] + 2*R_s[i]*(I_0[i] + I_phs[i]))
            X3[i] = C[i]* log(1 + I_phs[i]/I_0[i]) * (I_0[i] + I_phs[i])^2
            I_max[i] = (-X2[i] - sqrt(X2[i]^2 - 4*X1[i]*X3[i])) / (2*X1[i])
            V_max[i] = n*N*k*T_c[i]/q * log(1 + (I_phs[i] - I_max[i])/I_0[i]) - R_s[i] *I_max[i] 
            Ps[i] = V_max[i]*I_max[i] / 1.63 #power output per every 1.63 square meters
        end
    
        #Inefficiencies
        eR = 0.985 #reflection
        eS = 0.95 #soiling (sand/dirt getting on panels)
        eO = 0.97 #inverter 
        eW = 0.99 #Wiring
        eH = 0.99 #electricity to resistive heating
        Boost = 1.0 #used for determining necessary increase in efficiency to be competitive with PTC =1 by default

        for i = 1:length(GHI)
            Ps[i] = Ps[i] * eR * eS *eO * eW * eH * Boost
        end
        display(mean(Ps))
        display(mean(GHI))
        return Ps, GHI #power at SRC should be about 320 max
    end
    if gen == "PTC"
        display("PTC Power/GHI")
        r = 0.935 # mirror reflectance
        gct = 0.963 # glass cover transmittance
        abs = 0.96 # absorptance of receiver
        intF = 1 # intercept factor
        ra = 90 # rim angle [deg]
        apertureL = 6 # aperture length [m]
        mirrorL = 109.333 # mirror length [m]
        es = 0.98 # heat collector shadowing error factor
        et = 0.994 # tracking error factor
        eg = 0.98 # geometry error factor
        el = 0.96 # general losses error factor
        er = 0.88 # mirror reflectivity
        ed = er/r # dirt on mirror error factor
        edh = (1+ed)/2 # dirt on heat collector error factor

        absOD = 0.05 # outer-diameter of absorber tube [m]
        focalD = apertureL/(4*tan(ra/2*pi/180)) # focal distance [m]
        cr = apertureL/(pi*absOD) # concentration ratio
        apertureA = apertureL*mirrorL #aperture area [m^2]


        t = range(1,stop=length(DNI),length=length(DNI)) # time vector
        inc=ones(length(DNI)) # initialize incidence angle vector
        sun = zeros(length(DNI)) # is the sun up?
        # calculate the incidence angles for the desired location
        for i = 1:length(t) # loop over all times
            inc[i] = solarAngles(t[i],tz,lat,long, 1) # get incidence angles at this time. hard coded tracking for PTCs (last entry is a 1)
            if DNI[i]<0 inc[i]=0 end # reset to no-sun condition if dni data is negative
        end
        # area losses
        #=
        lossA = tan.(inc*pi/180)*(2/3*apertureL*focalD+apertureL*
                            focalD*(apertureL^2/(48*focalD^2)))
        =#
        # optical efficiency
        optEff = (cos.(inc*pi/180)+8.84E-4*inc*pi/180-5.369E-5*(inc*pi/180).^2)*
            es*et*eg*ed*edh*el*r*gct
        QabsA = optEff.*DNI*abs # W/m2 absorbed power density
        display(mean(QabsA))
        display(mean(DNI))
        return QabsA, DNI
        end

    if gen != "PV" || "PTC"
        display("Invalid Collector input, please use either PV or PTC")
    end
end

