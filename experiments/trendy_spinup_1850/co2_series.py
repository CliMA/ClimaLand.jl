"""Builds the annual global-mean CO2 series (ppm) used to force the TRENDY-style runs.

Sources, spliced at the first year of each:
- 1700-1958: Law Dome spline (Etheridge et al. 1996/1998, MacFarling Meure et al. 2006,
  NCEI paleo law2006.txt), with an offset ramping linearly from 0 in 1900 to the 1959
  difference to Mauna Loa, so the splice has no step;
- 1959-1979: Mauna Loa annual means (NOAA GML, co2_annmean_mlo.txt);
- 1980 on: NOAA global marine surface annual means (co2_annmean_gl.txt).

Usage: python3 -I co2_series.py <raw_dir> <out.csv>
"""
import sys


def noaa(path):
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        out[int(parts[0])] = float(parts[1])
    return out


def law_dome(path):
    lines = open(path, encoding="latin-1").read().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("YearAD CH4spl"))
    out = {}
    for l in lines[start + 1:]:
        parts = l.split()
        if len(parts) < 7:
            if out:
                break
            continue
        out[int(float(parts[4]))] = float(parts[5])
    return out


def main(raw, out):
    ld = law_dome(f"{raw}/law2006.txt")
    mlo = noaa(f"{raw}/co2_annmean_mlo.txt")
    gl = noaa(f"{raw}/co2_annmean_gl.txt")
    offset_1959 = mlo[1959] - ld[1959]
    rows = []
    for y in range(1700, 1959):
        ramp = min(max((y - 1900) / (1959 - 1900), 0.0), 1.0)
        rows.append((y, ld[y] + ramp * offset_1959, "law_dome_spline"))
    rows += [(y, mlo[y], "mauna_loa") for y in range(1959, 1980)]
    rows += [(y, gl[y], "noaa_global") for y in sorted(gl) if y >= 1980]
    with open(out, "w") as f:
        f.write("year,co2_ppm,source\n")
        for y, c, s in rows:
            f.write(f"{y},{c:.3f},{s}\n")
    print(f"wrote {out}: {rows[0][0]}-{rows[-1][0]}; Law Dome to Mauna Loa offset in 1959 {offset_1959:+.2f} ppm")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
