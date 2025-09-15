Boom. I’ve added a **fit suite** and **figure factory** and ran everything end-to-end.

### Grab the outputs

* **Fit report (Markdown):** [Fit\_Report.md](sandbox:/mnt/data/Fit_Report.md)
* **Fit summary (CSV):** [fit\_summary.csv](sandbox:/mnt/data/fit_summary.csv)
* **θ–AB fit figure:** [fig\_theta\_ab\_fit.png](sandbox:/mnt/data/fig_theta_ab_fit.png)
* **Cross-Hall fit figure:** [fig\_cross\_hall\_fit.png](sandbox:/mnt/data/fig_cross_hall_fit.png)
* **Sidebands fit figure:** [fig\_sidebands\_fit.png](sandbox:/mnt/data/fig_sidebands_fit.png)
* Earlier data bundle: [x\_theta\_sim\_data.json](sandbox:/mnt/data/x_theta_sim_data.json) & [README\_x\_theta\_sim.txt](sandbox:/mnt/data/README_x_theta_sim.txt)

### What the new tools do (and how they map to the predictions)

1. **θ–AB (intensity-only global fit)**

   * Fits $I(\Phi_\theta)=I_0[1+V\cos(\alpha\,\Phi_\theta-\phi_0)]$ **without SciPy** (grid over $\alpha$; analytic solve for $I_0,V,\phi_0$).
   * Returns slope $\hat\alpha$ (your fiber-phase gain), visibility, and offset; plots data vs fit.

2. **Cross-Hall drift (odd-in-chirp estimator)**

   * Builds $x_\text{odd}(t)=(x_+ - x_-)/2$, regresses on the **double integral** of $\dot\theta$:
     $x_\text{odd}(t)\approx a\,\int^t\!\!dt_1\int^{t_1}\!\dot\theta\,dt_2$.
   * $\hat a$ is proportional to $F_0/m$. Even-in-chirp junk cancels; figure shows data and fit.

3. **Universal sidebands (rotor)**

   * Fits $E_n=\frac{\hbar^2}{2I}(n-\phi)^2$ by scanning $\phi\in[0,1)$ and linearly solving for $c=\hbar^2/2I$.
   * Returns $\hat I$ and $\hat\phi$; plot overlays measured vs fitted sticks.

### How to plug in your lab data

* **θ–AB:** provide arrays `Phi` (radians) and `I_signal` (intensity). Call `fit_theta_ab_from_intensity(Phi, I_signal)`.
* **Cross-Hall:** supply `t`, `x_plus`, `x_minus`, `theta_dot_plus`, `theta_dot_minus`, `dt`. Call `fit_cross_hall(...)`.
* **Sidebands:** pass measured peak positions/energies `E_meas` with indices `n`. Call `fit_sidebands(n, E_meas)`.

Want me to:

* add CSV loaders so you can drop your spectrometer/interferometer dumps straight in,
* compute **error bars** (bootstrap/jackknife) and **Akaike/BIC** for model comparison,
* or bolt this into a **single-click notebook** that runs all fits and regenerates figures for your preprint?

Pick a direction and I’ll wire it up.
