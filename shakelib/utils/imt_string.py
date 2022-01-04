# local imports
import re


def oq_to_file(oqimt):
    """Convert openquake IMT nomenclature to filename friendly form.

    Examples:
    SA(1.0) (Spectral Acceleration at 1 second) -> psa1p0
    SA(0.3) (Spectral Acceleration at 0.3 second) -> psa0p3
    SA(15.0) (Spectral Acceleration at 15 seconds) -> psa15p0
    SA(3) (Spectral Acceleration at 3 seconds) -> psa3p0
    SA(.5) (Spectral Acceleration at 0.5 seconds) -> psa0p5

    Args:
        oqimt (str): Openquake IMT nomenclature string.
    Returns:
        str: Filename friendly IMT string.
    Raises:
        ValueError: when there is no corresponding filename-friendly
            IMT representation.
    """
    if oqimt in ["PGA", "PGV", "MMI"]:
        return oqimt.lower()
    float_pattern = r"[-+]?\d*\.\d+|\d+"
    periods = re.findall(float_pattern, oqimt)
    if not len(periods):
        fmt = 'IMT string "%s" has no file-name friendly representation.'
        raise ValueError(fmt % oqimt)
    period = periods[0]
    if period.find(".") < 0:
        integer = period
        fraction = "0"
    else:
        integer, fraction = period.split(".")
        if not len(integer):
            integer = "0"
    fileimt = f"psa{integer}p{fraction}"
    return fileimt


def file_to_oq(fileimt):
    """Convert filename friendly IMT form to openquake form.

    Examples:
    psa1p0 (Spectral Acceleration at 1 second) -> SA(1.0)
    psa0p3 (Spectral Acceleration at 0.3 second) -> SA(0.3)
    psa15p0 (Spectral Acceleration at 15 seconds) -> SA(15.0)

    Args:
        fileimt (str): Filename friendly IMT string.
    Returns:
        str: Openquake IMT nomenclature string.
    """
    if fileimt in ["pga", "pgv", "mmi"]:
        return fileimt.upper()
    if "p" not in fileimt:
        fmt = "%s is not a valid filename-friendly IMT string."
        raise ValueError(fmt % fileimt)
    integer, fraction = fileimt.replace("psa", "").split("p")
    if not len(fraction):
        fmt = "%s is not a valid filename-friendly IMT string."
        raise ValueError(fmt % fileimt)
    oqimt = f"SA({integer}.{fraction})"
    return oqimt
