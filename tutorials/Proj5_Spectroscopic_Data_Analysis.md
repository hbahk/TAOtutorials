# Project 5: Spectroscopic Data Analysis

Techniques in Astronomical Observation Project \#5<br>
Term: 1.5 week<br>
Prepared by Hyeonguk Bahk<br>

After the Big Bang, hydrogen underwent two significant phase changes. Recombination occurred at redshift $z \approx 1100$, allowing protons and electrons to combine into neutral hydrogen as the universe cooled down. Subsequently, reionization took place at $6 < z < 20$, driven by the emergence of the first stars and quasars. The study of the reionization era is crucial for understanding the early phase of star/galaxy formation. This era is characterized by the condensation of matter and the emission of high-energy photons, which ionize the surrounding neutral hydrogen. To investigate reionization at $z \approx 6$, we examine the fraction of neutral hydrogen in the diffuse intergalactic medium (IGM) using the absorption strength near the Lyman-alpha (Ly $\alpha$) line: i.e. **the Gunn-Peterson trough**. Gamma-ray burst (GRB) afterglows, with their bright luminosities detectable at $z \approx 6$, provide unique advantages for this analysis. GRBs have non-thermal synchrotron spectra with a simple power law, are located in less biased regions, and have short durations. Therefore, they are excellent tools for probing the ionization fraction during this critical period of the universe.

1. Download the provided raw spectral data of GRB 130606A or a similar astronomical target. The data should include the spectra of the target and a standard star for calibration. Then obtain a flux-calibrated 1D spectrum as follows.
    - a. Apply standard data reduction steps to the raw spectrum, including bias correction, flat fielding, and cosmic ray removal.
    - b. Extract the 1D spectrum from the 2D spectral image. This involves tracing the spectrum on the CCD, summing or averaging the flux in the spatial direction, and subtracting the background.
    - c. Perform wavelength calibration using known spectral lines from the calibration lamp images.
    - d. Flux-calibrate the extracted 1D spectrum using observations of a standard star. Convert the observed pixel counts to physical flux units.
    - e. Are there any additional steps required before applying scientific analysis to this target?

2. Measure the redshift of this source $z_{\rm s}=5.9131$, by identifying metal lines Si II $\lambda$1260.42, O I $\lambda$1302.17, and C II $\lambda$1334.53. Does your result match the provided value?

3. The continuum of GRB afterglow can be represented by the power-law model, as
   \begin{equation}
   \log f_{\nu} = \alpha + \nu \log \beta,    
   \end{equation}
   where $f_\nu$ is flux density at specific frequency $\nu$. The absorption profile due to the damped Ly $\alpha$ (DLA) region in host galaxy or in intervening IGM, can be obtained from the optical depth
   \begin{equation}
   \tau_{DLA}(\nu) = N_{HI} \sigma(\nu(1 + z_{DLA})),
   \end{equation}
   where $N_{\rm HI}$ is the column density of the neutral hydrogen, $\sigma(\nu)$ is the Ly $\alpha$ cross-section for the rest frame and $z_{\rm DLA}$ is the redshift of the intervening DLA region. The exact formula for $\sigma(\nu)$ is
   \begin{equation}
   \sigma(\nu) = \frac{3\lambda_{\alpha}^2{\Lambda_{\alpha}^2}}{8\pi} \frac{(\nu/\nu_{\alpha})^4}{4\pi^2(\nu - \nu_{\alpha})^2 + (\Lambda_{\alpha}^2/4)(\nu/\nu_{\alpha})^6},
   \end{equation}
   where $\lambda_{\alpha} = 1215.67$Å is the wavelength of the rest frame Ly $\alpha$ emission line, and $\Lambda_\alpha = 6.25 \times 10^8 {\rm \ s}^{-1}$ is the decay constant for Ly $\alpha$. Finally, the GP effect due to the IGM at the damping wing with optical depth at observed wavelength $\lambda = \lambda_\alpha (1 + z_{\rm s}) + \Delta \lambda$ is
   \begin{equation}
   \tau(\Delta\lambda) = \tau_{GP} R_{\alpha} \left(\frac{1 + \delta}{\pi}\right)^{3/2} (I(x_2) - I(x_1)),
   \end{equation}
   where $\delta = \lambda/[\lambda_\alpha (1+z_{\rm s})$, $x_1 = (1+z_{\rm IGM,l})/[(1+z_{\rm s})(1+\delta)]$, $x_2 = (1+z_{\rm IGM,u})/[(1+ z_{\rm s})(1 + \delta)]$, $R_\alpha = \Lambda_\alpha/(4\pi\nu_\alpha)$, $\tau_{\rm GP}=3.97 \times10^5f_{\rm HI}[(1 + z_{\rm s})/7]^{3/2}$,  and
   \begin{equation}
   I(x) = \frac{x^{9/2}}{1 - x} + \frac{9}{7}x^{7/2} + \frac{9}{5}x^{5/2} + 3x^{3/2} + 9x^{1/2} - \frac{9}{2} \ln\left(\frac{1 + x^{1/2}}{1 - x^{1/2}}\right).
   \end{equation}
   Here IGM is assumed to extend in a redshift range from $z_{\rm IGM,l}$ to $z_{\rm IGM,u}$ with the neutral hydrogen fraction $f_{\rm HI}$.
    - a. Measure the neutral fraction $f_{\rm HI}$ of IGM, considering two components of DLA absorption from the host galaxy and IGM extended in the redshift range of $5.67<z<z_{\rm s}$. Use the maximum likelihood estimation technique, or equivalently, the least square minimization technique. You are encouraged to try other methods, of course.
    - b. Discuss the constraining power on the value of $f_{\rm HI}$ and possible systematics of this analysis.
    - c. Discuss the implications of the obtained result.

## References
Madau, P., \& Rees, M. J. 2000 ApJ, 542, L69 <br>
Miralda-Escude, J. 1998, ApJ, 501, 15 <br>
Gunn, J. E., \& Peterson, B. A. 1965, ApJ, 142, 1633 <br>
Peebles, P. J. E. 1993, Principles of Physical Cosmology (Princeton: Princeton University Press) <br>
Schutz, B. 2009, A First Course in General Relativity (Cambridge University Press) <br>
Spergel, D. N., et al. 2003, ApJS, 148, 175 <br>
Totani, T., et al. 2006, PASJ, 58, 485 <br>
## Usefil Links
[An introduction to analysis of single dispersion spectra with IRAF](https://astro.snu.ac.kr/~hhwang/as_mon_1.pdf)