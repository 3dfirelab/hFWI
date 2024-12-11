import numpy as np
from types import SimpleNamespace
import pdb 

def DryingFactor_(Latitude, Month):
    return DryingFactor(np.array([Latitude]), Month, SimpleNamespace(**{'adjust_DryingFactor':'NSH'}))[0]
    #return DryingFactor(np.array([Latitude]), Month, SimpleNamespace(**{'adjust_DryingFactor':'original'}))[0]

def DayLength_(Latitude, numb_day, Month):
    return DayLength(np.array([Latitude]), numb_day, Month, SimpleNamespace(**{'adjust_DayLength':'bins'}))[0]
    #return DayLength(np.array([Latitude]), numb_day, Month, SimpleNamespace(**{'adjust_DayLength':'original'}))[0]


def DryingFactor(Latitude, Month, cfg):
    """
    Calculates the Drying Factor.
    NB: this parameter is called "Day-length adjustment in DC". This name is not used here because of the adjustments on the parameter "Effective day-length" used in DMC.

    Option for the drying factor from cfg.adjust_DryingFactor:
        - 'original': values for the Northern hemisphere applied everywhere (Wagner et al, 1987: https://cfs.nrcan.gc.ca/pubwarehouse/pdfs/19927.pdf)
        - 'NSH': the values for the Southern hemisphere are those for the northern shifted by 6 months (https://github.com/buckinha/pyfwi/blob/master/pyFWI/FWIFunctions.py)
        - 'NSHeq': the same idea is applied, but near the equator, one same value is applied for all months (https://rdrr.io/rforge/cffdrs/src/R/dcCalc.R)

    Parameters
    ----------
    Latitude: array
        Latitude in decimal degrees of the location (degree North)
    Month: int
        numeral month, from 1 to 12
    cfg: class configuration_FWI_CMIP6
        Class used to carry information on which calculations have to be performed, specifically on options for adjustments.

    Returns
    ----------
    DryingFactor: array
        Drying Factor (1)
    """

    if len(np.array(Month).shape) > 0:
        raise Exception(
            "This code has been rewritten so that only the month must be a scalar... dont know what to do if there are multiple values."
        )

    # LfN are the original values from the Canadian Fire Weather Index
    LfN = [-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6]
    # LfS are the same values shifted by 6 months
    LfS = [6.4, 5.0, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8]
    # Lfeq is the average
    Lfeq = 1.4

    if cfg.adjust_DryingFactor == "original":
        retVal = LfN[Month - 1] * np.ones(Latitude.shape)

    elif cfg.adjust_DryingFactor == "NSH":
        # setting drying factor at different latitudes
        retVal = np.nan * np.ones(Latitude.shape)
        ind_LatG0 = np.where(Latitude > 0)
        retVal[ind_LatG0] = LfN[Month - 1]
        ind_LatL0 = np.where(Latitude <= 0)
        retVal[ind_LatL0] = LfS[Month - 1]

    elif cfg.adjust_DryingFactor == "NSHeq":
        # setting drying factor at different latitudes
        retVal = np.nan * np.ones(Latitude.shape)
        ind_LatG0 = np.where(Latitude > 20)
        retVal[ind_LatG0] = LfN[Month - 1]

        ind_LatG0 = np.where((Latitude > -20) & (Latitude <= 20))
        retVal[ind_LatG0] = Lfeq

        ind_LatL0 = np.where(Latitude <= -20)
        retVal[ind_LatL0] = LfS[Month - 1]

    else:
        raise Exception("Unknown name for the type of drying factor")

    return retVal


def DayLength(Latitude, numb_day, MONTH, cfg):
    """
    Calculates the effective Day Length
    Option for the effective day length cfg.adjust_DayLength:
        - 'original': values adapted for Canadian latitudes, depends on the month (Wagner et al, 1987: https://cfs.nrcan.gc.ca/pubwarehouse/pdfs/19927.pdf)
        - 'bins': depends on 4 bins of latitudes and the month (https://github.com/buckinha/pyfwi/blob/master/pyFWI/FWIFunctions.py)
        - 'continuous': depends continuously on latitude and the day of the year (https://github.com/NCAR/fire-indices & https://www.ncl.ucar.edu/Document/Functions/Crop/daylight_fao56.shtml)


    Parameters
    ----------
    Latitude: array
        Latitude in decimal degrees of the location (degree North)
    numb_day: scalar
        Number of the day in the year
    Month: int
        numeral month, from 1 to 12
    cfg: class configuration_FWI_CMIP6
        Class used to carry information on which calculations have to be performed, specifically on options for adjustments.

    Returns
    ----------
    DayLength: array
        Day Length (h)
    """
    if cfg.adjust_DayLength == "original":
        DayLength46N = [6.5, 7.5, 9.0, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8.0, 7.0, 6.0]
        retVal = DayLength46N[MONTH - 1] * np.ones(Latitude.shape)

    elif cfg.adjust_DayLength == "bins":
        # preparing values of day length
        DayLength46N = [6.5, 7.5, 9.0, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8.0, 7.0, 6.0]
        DayLength20N = [7.9, 8.4, 8.9, 9.5, 9.9, 10.2, 10.1, 9.7, 9.1, 8.6, 8.1, 7.8]
        DayLength20S = [10.1, 9.6, 9.1, 8.5, 8.1, 7.8, 7.9, 8.3, 8.9, 9.4, 9.9, 10.2]
        DayLength40S = [11.5, 10.5, 9.2, 7.9, 6.8, 6.2, 6.5, 7.4, 8.7, 10.0, 11.2, 11.8]

        # setting day length at different latitudes
        retVal = np.nan * np.ones(Latitude.shape)
        ind_latG33L90 = np.where((33 < Latitude) & (Latitude <= 90))
        retVal[ind_latG33L90] = DayLength46N[MONTH - 1]
        ind_latG0L33 = np.where((0 < Latitude) & (Latitude <= 33))
        retVal[ind_latG0L33] = DayLength20N[MONTH - 1]
        ind_latGm30L0 = np.where((-30 < Latitude) & (Latitude <= 0))
        retVal[ind_latGm30L0] = DayLength20S[MONTH - 1]
        ind_latGm90Lm30 = np.where((-90 <= Latitude) & (Latitude <= -30))
        retVal[ind_latGm90Lm30] = DayLength40S[MONTH - 1]

    elif cfg.adjust_DayLength == "continuous":
        lat = Latitude * np.pi / 180  # degree -> radian
        sun_dec = 0.409 * np.sin(2 * np.pi / 365 * numb_day - 1.39)  # equation 24
        # preparing equation 25, with special cases handled
        val_for_arccos = -np.tan(lat) * np.tan(sun_dec)
        val_for_arccos[np.where(val_for_arccos < -1)] = -1
        val_for_arccos[np.where(val_for_arccos > 1)] = 1
        sunset_hour_angle = np.arccos(val_for_arccos)  # equation 25
        retVal = 24 / np.pi * sunset_hour_angle  # equation 34

    else:
        raise Exception("Unknown name for the type of day length")

    if np.any(np.isnan(retVal)):
        raise InvalidLatitude(Latitude)

    return retVal

'''
def daylength(dayOfYear, lat):
    """Computes the length of the day (the time between sunrise and
    sunset) given the day of the year and latitude of the location.
    Function uses the Brock model for the computations.
    For more information see, for example,
    Forsythe et al., "A model comparison for daylength as a
    function of latitude and day of year", Ecological Modelling,
    1995.
    Parameters
    ----------
    dayOfYear : int
        The day of the year. 1 corresponds to 1st of January
        and 365 to 31st December (on a non-leap year).
    lat : float
        Latitude of the location in degrees. Positive values
        for north and negative for south.
    Returns
    -------
    d : float
        Daylength in hours.
    """
    latInRad = np.deg2rad(lat)
    declinationOfEarth = 23.45*np.sin(np.deg2rad(360.0*(283.0+dayOfYear)/365.0))
    if -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) <= -1.0:
        return 24.0
    elif -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) >= 1.0:
        return 0.0
    else:
        hourAngle = np.rad2deg(np.arccos(-np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth))))
        return 2.0*hourAngle/15.0

'''
