# CalLMIP Phase 1b calibration-site table.
#
# Values extracted from the PLUMBER2 met files and verified in
# DATA_MANIFEST.md. UTC offsets (UTC = local - offset) were confirmed
# empirically from the SWdown sunrise/sunset midpoint at every site; they are
# NOT derivable from longitude (FR-Pue, NL-Loo, RU-Fyo, US-MMS differ from
# round(lon/15)) and are NOT stored in the met files.

const CALIBRATION_SITES = [
    (
        id = "CA-Qfo",
        lat = 49.6925,
        long = -74.3421,
        utc_offset = -5,
        zref = 24.0,
        met_years = (2004, 2010),
    ),
    (
        id = "CH-Dav",
        lat = 46.8153,
        long = 9.8559,
        utc_offset = 1,
        zref = 35.0,
        met_years = (1997, 2014),
    ),
    (
        id = "DE-Gri",
        lat = 50.9495,
        long = 13.5125,
        utc_offset = 1,
        zref = 3.0,
        met_years = (2004, 2014),
    ),
    (
        id = "DE-Hai",
        lat = 51.0792,
        long = 10.4530,
        utc_offset = 1,
        zref = 43.5,
        met_years = (2000, 2012),
    ),
    (
        id = "DE-Tha",
        lat = 50.9636,
        long = 13.5669,
        utc_offset = 1,
        zref = 42.0,
        met_years = (1998, 2014),
    ),
    (
        id = "DK-Sor",
        lat = 55.4859,
        long = 11.6446,
        utc_offset = 1,
        zref = 57.0,
        met_years = (1997, 2014),
    ),
    (
        id = "FI-Hyy",
        lat = 61.8475,
        long = 24.2950,
        utc_offset = 2,
        zref = 23.0,
        met_years = (1996, 2014),
    ),
    (
        id = "FR-Pue",
        lat = 43.7414,
        long = 3.5958,
        utc_offset = 1,
        zref = 11.0,
        met_years = (2000, 2014),
    ),
    (
        id = "IT-Lav",
        lat = 45.9562,
        long = 11.2813,
        utc_offset = 1,
        zref = 33.0,
        met_years = (2005, 2014),
    ),
    (
        id = "IT-MBo",
        lat = 46.0147,
        long = 11.0458,
        utc_offset = 1,
        zref = 2.5,
        met_years = (2003, 2012),
    ),
    (
        id = "IT-Noe",
        lat = 40.6062,
        long = 8.1512,
        utc_offset = 1,
        zref = 3.0,
        met_years = (2004, 2014),
    ),
    (
        id = "NL-Loo",
        lat = 52.1666,
        long = 5.7436,
        utc_offset = 1,
        zref = 27.0,
        met_years = (1997, 2013),
    ),
    (
        id = "RU-Fyo",
        lat = 56.4615,
        long = 32.9221,
        utc_offset = 3,
        zref = 48.0,
        met_years = (2003, 2014),
    ),
    (
        id = "US-MMS",
        lat = 39.3232,
        long = -86.4131,
        utc_offset = -5,
        zref = 48.0,
        met_years = (1999, 2014),
    ),
    (
        id = "US-NR1",
        lat = 40.0329,
        long = -105.5464,
        utc_offset = -7,
        zref = 26.0,
        met_years = (1999, 2014),
    ),
    (
        id = "US-SRG",
        lat = 31.7894,
        long = -110.8277,
        utc_offset = -7,
        zref = 3.25,
        met_years = (2009, 2014),
    ),
    (
        id = "US-SRM",
        lat = 31.8214,
        long = -110.8660,
        utc_offset = -7,
        zref = 6.4,
        met_years = (2004, 2014),
    ),
    (
        id = "US-Ton",
        lat = 38.4316,
        long = -120.9660,
        utc_offset = -8,
        zref = 23.0,
        met_years = (2001, 2014),
    ),
    (
        id = "US-Var",
        lat = 38.4133,
        long = -120.9507,
        utc_offset = -8,
        zref = 3.0,
        met_years = (2001, 2014),
    ),
    (
        id = "US-Whs",
        lat = 31.7438,
        long = -110.0522,
        utc_offset = -7,
        zref = 4.0,
        met_years = (2008, 2014),
    ),
    (
        id = "US-Wkg",
        lat = 31.7365,
        long = -109.9419,
        utc_offset = -7,
        zref = 6.4,
        met_years = (2005, 2014),
    ),
]

const CALIBRATION_SITE_IDS = [s.id for s in CALIBRATION_SITES]

# Validation sites (protocol Table 2), met from ME-org "Phase 1b Validation".
# UTC offsets verified 2026-08-18 from SWdown sunrise/sunset midpoints (the AU
# sites sit on the +9.5 half-hour zone); AU-Tum, US-Ha1, US-UMB report HOURLY
# met (handled by resample_to_30min). Provenance: OzFlux (AU), LaThuile
# (US-FPe), FLUXNET2015 (rest).
const VALIDATION_SITES = [
    (
        id = "AU-ASM",
        lat = -22.2830,
        long = 133.2490,
        utc_offset = 9.5,
        zref = 11.6,
        met_years = (2011, 2017),
    ),
    (
        id = "AU-How",
        lat = -12.4952,
        long = 131.1500,
        utc_offset = 9.5,
        zref = 23.0,
        met_years = (2003, 2017),
    ),
    (
        id = "AU-Stp",
        lat = -17.1507,
        long = 133.3502,
        utc_offset = 9.5,
        zref = 4.8,
        met_years = (2010, 2017),
    ),
    (
        id = "AU-Tum",
        lat = -35.6566,
        long = 148.1517,
        utc_offset = 10,
        zref = 70.0,
        met_years = (2002, 2017),
    ),
    (
        id = "CH-Cha",
        lat = 47.2102,
        long = 8.4104,
        utc_offset = 1,
        zref = 2.0,
        met_years = (2006, 2014),
    ),
    (
        id = "US-FPe",
        lat = 48.3077,
        long = -105.1019,
        utc_offset = -7,
        zref = 3.5,
        met_years = (2000, 2006),
    ),
    (
        id = "US-GLE",
        lat = 41.3665,
        long = -106.2399,
        utc_offset = -7,
        zref = 23.0,
        met_years = (2009, 2014),
    ),
    (
        id = "US-Ha1",
        lat = 42.5378,
        long = -72.1715,
        utc_offset = -5,
        zref = 30.0,
        met_years = (1992, 2012),
    ),
    (
        id = "US-Me2",
        lat = 44.4523,
        long = -121.5574,
        utc_offset = -8,
        zref = 29.0,
        met_years = (2002, 2014),
    ),
    (
        id = "US-UMB",
        lat = 45.5598,
        long = -84.7138,
        utc_offset = -5,
        zref = 50.0,
        met_years = (2000, 2014),
    ),
]

const VALIDATION_SITE_IDS = [s.id for s in VALIDATION_SITES]

site_info(id) =
    only(filter(s -> s.id == id, [CALIBRATION_SITES; VALIDATION_SITES]))
