from . import cams, wrfout

GLOBAL_MODELS = {
    "cams_eac4": cams.CAMS_EAC4,
    "cams_global_forecasts": cams.CAMS_Global_Forecasts,
    "cams_eac4_pl": cams.CAMS_EAC4_Pressure,
    "cams_global_forecasts_pl": cams.CAMS_Global_Forecasts_Pressure,
    "wrfout": wrfout.WRFOutput,
}
