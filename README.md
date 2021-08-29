# Exoplanet Mass limits for Radial Velocity data 
What does the program do?
1) It accepts input in the form of time series radial velocity data for a particular star.
2) It computes the Lomb-Scargle periodogram of this time-series.
3) It adds an artificial signal to the original data, as if generated by an orbiting planet/companion object. You control the planet-mass and orbital frequency.
4) MassSlider.py creates an interactive widget which allows you to smoothly change the parameters of the artificial signal (mass, frequency and phase) and see the changed periodogram in real time.
![](img/slider.png)
5) GradMap.py computes the detection probability (fraction of phases for which the power of modified data exceeds the false-alarm probability threshold) for a range of masses and frequencies.
![](img/map.png)
6) It neatly saves all the data in folders and creates a gradient 2D map of detection probability.
