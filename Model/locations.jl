function locations(place)
    #Code used to simplify the location choosing process. Enter the name of the city as a string and it will output which csv file to use and the timezone, Lattitude, and Longitude
    #Made by Justin Rastinejad on 05/28/2022
    
    if place == "Firebaugh"
        #avg GHI: 
        input_file = "TMYFirebaugh_2020.csv"
        tz = -8
        Lat = 36.85
        Long = -120.46
        #ng = (13.01)/293.1/1.037 #california industrial rate 
        #https://www.eia.gov/dnav/ng/ng_pri_sum_dcu_SCA_m.htm
    end

    if place == "Aurora"
        input_file = "TMYAurora_2020.csv"
        tz = -7
        Lat = 39.73
        Long = -104.83
    end

    if place == "Renton"
        input_file = "TMYRenton_2020.csv"
        tz = -8     
        Lat = 47.49
        Long = -122.22
    end

    if place =="Weston"
        input_file = "TMYWeston_2020.csv"
        tz = -5
        Lat = 42.37
        Long = -71.30
    end
    
    return input_file, tz, Lat, Long
end
