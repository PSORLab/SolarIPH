function pvAngles(h,tzone,lat,long)
  #=============================================================================
  This was bad and I no longer use it -Justin 8/2/21 
  Based off of https://www.homerenergy.com/products/pro/docs/latest/how_homer_calculates_the_radiation_incident_on_the_pv_array.html
    
    h: hour of year (1:8760)
    tzone: time-zone
    long: longitude
    lat: latitude
   OUTPUT
    ia: solar incidence angle
  =============================================================================#
  input_file = "Zenith Angle 2020 Data.csv"
  Data = CSV.File(input_file) |> DataFrame
  ZAs = convert(Array{Float64,1}, Data[:,10]) #Zenith Angle Data
  theta = zeros(length(h))
  for i = 1:length(h)
    n = floor(h[i]/24+1) #current day
    LT = h[i]%24 
    decl = 23.45*pi/180*sin(2*pi*(284+n)/365) #Solar declination angle
    B = 2*pi*(n-81)/365
    EoT = 9.87*sin(2*B)- 7.53*cos(B) - 1.5*sin(B) #3.82*(0.000075 + 0.001868*cos(B) - 0.32077*sin(B) - 0.014615*cos(2*B) - 0.04089*sin(2*B)) #Equation of time
    TC = 4*(long-15*tzone) #Time Correction Factor
    LST = LT + TC/60 #Local Solar Time
    ha = (LST - 12)*15/180*pi #hour angle (zero at noon)
    theta[i] = 180/pi*acos(cos(ZAs[i]*pi/180)*cos(ha) + cos(decl)*(sin(ha)^2)) #taken from Meinel and Meinel, 1976 for a horizontal N-S axis with E-W tracking
    display(theta[i])
  end
  
  display(theta)
end;
