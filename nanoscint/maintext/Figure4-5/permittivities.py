def eps_YAG(wl):
    # wl wavelength in microns
    # From Sellmeier coefficients
    eps = 1 + 2.28200*wl**2./(wl**2.-0.01185) + 3.27644*wl**2./(wl**2.-282.734)    
    return eps 
