

void conductionQ(int nz, double z[], double dt, double Qn, double Qnp1,
		 double T[], double ti[], double rhoc[], double emiss,
		 double Fgeotherm, double *Fsurf);

void conductionT(int nz, double z[], double dt, double T[], double Tsurf,
		 double Tsurfp1, double ti[], double rhoc[], double Fgeotherm,
		 double *Fsurf);

double flux_noatm(double R, double decl, double latitude, double HA,
		  double surfaceSlope, double azFac);

void heatflux_from_temperature(int nz, double z[], double T[],
			       double k[], double H[]);

void setgrid(int nz, double z[], double zmax, double zfac);

void tridag(double a[], double b[], double c[], double r[], double T[],
	    unsigned long nz);
