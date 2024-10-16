class planetool:

    def magnetic_field(self, R, m):
        pi = 3.14159
        mu_0 = 4 * pi * 10**(-7)
        # R: planet's radius
        # B: magnetic field intensity
        B = ((self.mu_0) / (4 * self.pi)) * (self.m / R**3)
        return B

    def coriolis(self, w, phi, mh, vh):
        radian_phi = np.radians(self.phi)
        # Calculate sin(phi)
        sin_phi = np.sin(radian_phi)
        # Fc: coriolis force
        Fc = 2 * mh * vh * w * sin_phi
        return Fc

    def tidal_force(self, Mc, m, rc, R):
        G = 6.67 * 10**(-11)
        # Ftital: tidal force
        F_total = 2 * ((self.G * self.Mc * self.m) / (self.rc**3)) * R
        return F_total

    def wave_energy(self, R, ro, A, m):
        G = 6.67 * 10**(-11)
        g = self.G * self.m / (R**2)
        # Ed: wave energy
        Ed = (1 / 2) * ro * g * (A**2)
        return Ed

    def pressure(self, R, m, P0, ro, h):
        G = 6.67 * 10**(-11)
        g = self.G * self.m / (R**2)
        # P: pressure
        P = P0 + self.h * ro * self.g
        return P

    def heat_flow(self, k, x, y, z):
        x, y, z = sp.symbols('x y z')
        f = x**2 + y**3 + z**2
        df_dx = sp.diff(f, x)  # derivative in x direction
        df_dy = sp.diff(f, y)  # derivative in y direction
        df_dz = sp.diff(f, z)  # derivative in z direction
        # q: heat flow
        q = (-self.k) * (df_dx + df_dy + df_dz)
        return q

    def gibbs_energy(self, DeltaH, T, DeltaS):
        # DeltaG: change in Gibbs free energy
        DeltaG = self.DeltaH - (self.T * self.DeltaS)
        return DeltaG

    def enthalpy_change(self, H_products, H_reactants):
        # DeltaH: change in enthalpy
        DeltaH = sum(H_products) - sum(H_reactants)
        return DeltaH

    def surface_area(self, R):
        pi = 3.14159
        # A: surface area
        A = 4 * self.pi * (R**2)
        return A

    def ocean_volume(self, R, ocean_fraction):
        pi = 3.14159
        # V_ocean: ocean volume
        V_ocean = (4 / 3) * self.pi * (R**3) * self.ocean_fraction
        return V_ocean

    def escape_velocity(self, R, M):
        G = 6.67 * 10**(-11)
        # V_escape: escape velocity
        V_escape = math.sqrt((2 * self.G * M) / R)
        return V_escape
