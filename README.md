# FluxSeparator - Separation of diffusive and ebullitive fluxes


This R-package is dedicated to separation and calculation of diffusive and ebullitive (bubble) fluxes. 
Determination of ebullitive events is characterized by sudden increase in concentration over a short period of time. Thus, to determine ebullitive events a running variance approach is used. If the running variance is above a user-set threshold value, the data is considered an ebullitive event. Ebullitive events can then be calculated as the difference from the lowest concentration value to either the highest or the last value. Furthermore, correct identification of ebullitive events can be visually inspected when using the function. The functions also allows for several parameters to be changed, for full list of parameters see the R documentation for each function (called by ?_function_). The function for diffusive flux also looks for ebullitive events, to avoid diffusive flux being calculated post an ebullitive event. This function can likewise be turned off.   


```
remotes::install_github('JonasStage/FluxSeparator')
```

For citing this package follow the instructions on https://zenodo.org/doi/10.5281/zenodo.8297153

Sø. J.S. (2023). JonasStage/FluxSeparator: v1.0.0 (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.8297154


### Citations

Sø, J. S., Sand-Jensen, K., Martinsen, K. T., Polauke, E., Kjær, J. E., Reitzel, K., & Kragh, T. (2023). Methane and carbon dioxide fluxes at high spatiotemporal resolution from a small temperate lake. Science of the Total Environment, 878, 162895. doi:https://doi.org/10.1016/j.scitotenv.2023.162895

