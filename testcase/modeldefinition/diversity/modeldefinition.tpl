ptf #
{
    "Simulation": {
        "Name": "2016",
        "Start": "01.01.2016",
        "End": "31.12.2016",
        "dt": 0.01,
        "repetitions": 3,
        "thinning": 10
    },
    "Morphology": {
        "Maxdepth": 16.0,
        "Bathymetry": {
            "zareas": [
                0.0,
                2.0,
                4.0,
                6.0,
                8.0,
                10.0,
                12.0,
                14.0,
                16.0
            ],
            "areas": [
                494199.0,
                440930.0,
                393912.0,
                348221.0,
                301994.0,
                240935.0,
                153063.0,
                69403.0,
                18691.0
            ]
        },
        "Gridcells": 500
    },
    "Parameter": {
        "sediment_exchange_k": #sediment_exchange_k#,
        "sediment_diffusivity": #sediment_diffusivity#,
        "sediment_vmax": #sediment_vmax#,
        "OM_sedimentation_rate": #OM_sedimentation_rate#,
        "sediment_kI_o2": #sediment_kI_o2#,
        "o2_consumption_rate": #o2_consumption_rate#
    },
    "Forcing": {
        "Weather": "E:/polybox/SURF/MobTrait/results/diversity/rotsee_forcing.dat",
        "fwind": 0.5,
        "Inflow": "E:/polybox/SURF/MobTrait/results/diversity/rotsee_inflow.dat",
        "Simstrat calibration": "E:/polybox/SURF/MobTrait/results/physicalmodel/seasonal/3_methane_iteration_revised/Rotsee_Output/"
    }
}