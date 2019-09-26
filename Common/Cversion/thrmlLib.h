

void conductionQ(int, double [], double, double, double,
		double [], double [], double [], double, double, double *);

double flux_noatm(double R, double decl, double latitude, double HA,
		  double surfaceSlope, double azFac);

void setgrid(int nz, double z[], double zmax, double zfac);

void tridag(double a[], double b[], double c[], double r[], double T[], unsigned long nz);
