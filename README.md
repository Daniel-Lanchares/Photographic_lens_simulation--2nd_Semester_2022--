# Photographic_lens_simulation--2nd_Semester_2022--
Project from the "Optics" course. 2nd semester of the 2021/2022 school-year.

By ray-tracing using the paralax aproximation $(\sin(\theta) \approx \theta$, the angle between the ray and the optical axis $\left. \right)$ we can smodel rays as slope-height vectors, $(m, h)$, and optical elements as matrices. Contrary to spatial representations, translations are also modelled as matrices acting on these vectors. The refraction indices appear color coded for better understanding, with darker blues corresponding to larger values:

$$
\vec{r} = \left( 
\begin{align}
&m\\
&h
\end{align}
\right)
\text{ }
M = \left( 
\begin{matrix}
 1 & 0\\
\frac{1-n_{01}}{R_{01}} & n_{01}
\end{matrix}
\right)
\text{            }
T = \left( 
\begin{matrix}
 1 & d\\
0 & 1
\end{matrix}
\right)
$$

$$
\begin{align*}
m &: \text{slope of the ray} &n_{01} &: \text{refraction index quotient }n_0/n_1 &d &: \text{translated distance projection}\\
h &: \text{height of the ray} &R_{01} &: \text{signed radius of optical element} & & \text{onto horizontal (optical axis)}
\end{align*}
$$

With careful tuning this simple formalism allows us to mimic real life lenses with surprising accuracy, obtaining close aproximations of their focal lenght and maximum apperture, as seen in this example of a 90mm f/2.8 double Gauss lens:

![Double Gauss](<https://raw.githubusercontent.com/Daniel-Lanchares/Photographic_lens_simulation--2nd_Semester_2022--/main/Pictures/Apperture simulations/Diafragma f2.8.png>)

This particular example was constructed to explain the effect of wide appertures on objects increasingly far from the focusing plane (here at infinity), causing them to mix in an out of focus blur commonly refered to as bokeh.

We can also use this lens setup to test focus itself. In this next image the focus has been set to its physical minimum (if we try to focus on objects closer than 1720mm the apperture would run into a lens element), which is made apparent by the red plane representing the focusing plane of objects and infinity, which no longer coincides with the image sensor / analog film. With a keen eye one can observe a tiny bit of lens breathing on the zoomed in section, as the lilac rays are no longer at the edge of the sensor. This problem becomes much more aparent in zoom lenses, as the number of moving lens groups increases.

![Double Gauss near focused](<https://raw.githubusercontent.com/Daniel-Lanchares/Photographic_lens_simulation--2nd_Semester_2022--/main/Pictures/Focusing simulations/Enfoque 1720mm.png>)

While fixed length lenses are still popular and preferable in some situations most lenses today have some sort of optical zoom technology (even those on our phones, if often complemented by digital zooming). For this demonstration we have a replica of a 35-70mm lens, thouch as it is probably clear no attempt was made to model a physical apperture, instead focusing on the focal length only.

![Zoom 35mm](<https://raw.githubusercontent.com/Daniel-Lanchares/Photographic_lens_simulation--2nd_Semester_2022--/main/Pictures/35-70 Zoom lens/1-Zoom 35.png>) 
![Zoom 70mm](<https://raw.githubusercontent.com/Daniel-Lanchares/Photographic_lens_simulation--2nd_Semester_2022--/main/Pictures/35-70 Zoom lens/3-Zoom 70.png>)

The change in cone of vision is represented by the dotted red lines, with the green rays demonstrating that an object that covers about half the image sensor height at 35mm will cover its entirety at 70mm. The various lens groups of this particular design appear denoted, though most zoom lenses use them as well. The colector and master groups usually stay still while the moving group, well..., moves and compensator also moves to ensure the focus remains unchanged.

Finally we took advantage of our flexible coding aproach to introduce chromatic aberrations and a demonstration of how it is corrcted. By introducing a dependency on wavelength to our refraction indices (and promoting our matrices to matrix functions under the hood) we can model dispersion, and how certain lens combinations can mitigate its effects. Now the colour of each ray represents its actual colour approximately.

![Acromat lens combo](<https://raw.githubusercontent.com/Daniel-Lanchares/Photographic_lens_simulation--2nd_Semester_2022--/main/Pictures/Chromatic aberration correction.png>)
