# Zenith Angle Equations
These equations are from Tapio Schneider and Lenka Novak's textbook draft, chapter 3.

## Mean Anomaly
The mean anomaly $M$ at current time $t$ is
```math
M = \frac{2\pi (t - t_0)}{Y_a} + M_0,
```
where $t_0$ is the time at the epoch (J2000), defined as January 1, 2000 at 12hr UTC, 
$M_0$ is the mean anomaly at the epoch, and $Y_a$ is the length of the anomalistic year.

## True Anomaly
The true anomaly $A$ is given by
```math
A = M + \left( 2e - \frac{1}{4}e^{3} \right) \sin(M) + \frac{5}{4} e^2 \sin(2M) + \frac{13}{12} e^3 \sin(3M).
```

## True Longitude
The true longitude is the sum of the true anomaly and the longitude of perihelion ($\varpi$),
```math
L = A + \varpi.
```

## Declination
The sine of the declination angle is
```math
\sin \delta = \sin \gamma \sin L,
```
where $\gamma$ is the orbital obliquity.

## Equation of Time 
The equation of time corrects for the difference between apparent solar time (local solar noon)
and mean solar time due to the obliquity of the planet and eccentricity of the orbit,
```math
\Delta t = -2 e \sin(M) + \tan^2(\gamma/2) \sin(2M+2\varpi).
```

## Hour Angle
The hour angle $\eta$ measures the longitude difference between the subsolar point position at time $t$ and the current location.
The hour angle is calculated with respect to the prime meridian and then the longitude is added to calculate the hour angle at the given location.
```math
\eta = \eta_UTC + \lambda = \frac{2\pi (t+\Delta t)}{T_d} + \lambda,
```
where $t$ is the current time referenced to the epoch $t_0$, $T_d$ is the length of a day, and $\lambda$ is the longitude.

## Zenith Angle
The zenith angle at any time is given by,
```math
\cos \theta = \cos \phi \cos \delta \cos \eta + \sin \phi \sin \delta,
```
where $\phi$ is the latitude.

## Sunrise/Sunset Angle
The sunrise/sunset angle is the hour angle $\eta_d$ at which the sun rises or sets,
```math
\cos \eta_d = - \tan \phi \tan \delta.
```

## Diurnally averaged Insolation
The diurnally averaged insolation is obtained from the averaged cosine of the zenith angle, computed in terms of the sunrise/sunset hour angle as
```math
\overline{\cos \theta} = \frac{1}{\pi} \left( \eta_d \sin \phi \sin \delta + \cos \phi \cos \delta \cos \eta_d \right).
```

## Azimuth Angle
The azimuth angle is
```math
\zeta = \frac{3\pi}{2} - \arctan \left( \frac{\sin \eta}{\cos \eta \sin \phi - \tan \delta \cos \phi} \right).
```
The azimuth is defined as 0 to the East and increasing counter-clockwise, such that at local solar noon when $\eta=0$, then $\zeta = \frac{3\pi}{2}$.

## Planet-Sun Distance
The distance between the planet and the sun depends on the eccentricity of the orbit $e$, true anomaly $A$, and mean planet-sun distance $d_0$ through
```math
d = \frac{1-e^2}{1+e\cos A} d_0
```
For the Earth, $d_0 = 1$ AU.