# Laskar et. al. (2004) data

This is repackaged data from http://vo.imcce.fr/insola/earth/online/earth/earth.html, specifically:
- [`INSOLN.LA2004.BTL.ASC`](http://vo.imcce.fr/insola/earth/online/earth/La2004/INSOLN.LA2004.BTL.ASC)
- [`INSOLP.LA2004.BTL.ASC`](http://vo.imcce.fr/insola/earth/online/earth/La2004/INSOLP.LA2004.BTL.ASC)

From the following paper:
> Laskar, J., Robutel, P., Joutel, F., Gastineau, M., Correia, A.C.M., Levrard, B. 
> A long term numerical solution for the insolation quantities of the Earth.
> A&A 428, 261-285 (2004), DOI: 10.1051/0004-6361:20041335

The two files are combined into one, with columns
1. Time from J2000 in 1000 years
2. eccentricity
3. obliquity (radians)
4. longitude of perihelion from moving equinox (radians)
    - modified from the original files to be continuous over time
