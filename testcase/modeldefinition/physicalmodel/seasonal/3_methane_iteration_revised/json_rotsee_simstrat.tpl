ptf #
{
    "Input": {
        "Initial conditions": "Rotsee_Input/rotsee_initialconditions.dat",
        "Grid": "Rotsee_Input/rotsee_grid.dat",
        "Morphology": "Rotsee_Input/rotsee_morphology.dat",
        "Forcing": "Rotsee_Input/rotsee_forcing.dat",
        "Absorption": "Rotsee_Input/rotsee_absorption.dat",
        "Inflow": "Rotsee_Input/rotsee_Qinp.dat",
        "Outflow": "Rotsee_Input/rotsee_Qout.dat",
        "Inflow temperature": "Rotsee_Input/rotsee_Tinp.dat",
        "Inflow salinity": "Rotsee_Input/rotsee_Sinp.dat"
    },
    "Output": {
        "Path": "Rotsee_Output/",
        "Depths": "Rotsee_Input/rotsee_outputdepths.dat",
        "OutputDepthReference": "surface",
        "Times": "Rotsee_Input/rotsee_outputtime.dat"
    },
    "ModelConfig": {
        "MaxLengthInputData": 1000,
        "CoupleAED2": false,
        "TurbulenceModel": 1.0,
        "StabilityFunction": 1.0,
        "FluxCondition": 1.0,
        "Forcing": 3.0,
        "UseFilteredWind": false,
        "SeicheNormalization": 2.0,
        "WindDragModel": 3.0,
        "InflowPlacement": 0.0,
        "PressureGradients": 0.0,
        "IceModel": 0,
        "SnowModel": 0
    },
    "Simulation": {
        "Timestep s": 300.0,
        "Start d": 16147.0,
        "End d": 17198.0,
        "DisplaySimulation": 0
    },
    "ModelParameters": {
        "lat": #lat       #,
        "p_air": #p_air     #,
        "a_seiche": #a_seiche  #,
        "q_nn": #q_NN      #,
        "f_wind": #f_wind    #,
        "c10": #C10       #,
        "cd": #CD        #,
        "hgeo": #fgeo      #,
        "k_min": #k_min     #,
        "p_radin": #p1        #,
        "p_windf": #p2        #,
        "beta_sol": #beta      #,
        "beta_snowice": 0.4,
        "albsw": #albsw     #,
        "ice_albedo": 0.3,
        "snow_albedo": 0.77,
        "freez_temp": 0.05
    }
}