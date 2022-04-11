# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:36:58 2022

@author: micha
"""
from math import sin, cos, sqrt, degrees, atan, pi, tan, radians, atan2, asin, acos
import numpy as np


   
class Transformacje:
    def __init__(self, model: str = "wgs84"):
        #https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.ecc2 = (2 * self.flattening - self.flattening ** 2)
        print(model,self.b)
        
        
        
    def deg2rad(self, kat):
        '''
        Funkcja przelicza kąt w stopniach na radiany.
        
        Argumenty:
        ----------
        kat - kąt podany w stopniach | typ: float lub int
        
        Wyniki:
        -------
        katr - kąt przeliczony z argumentu 'kat' na radiany | typ: float
        
        Dodatkowy opis:
        ---------------
        '''
        katr = kat * pi / 180
        return(katr)
    
    def rad2deg(self, kat):
        '''
        Funkcja przelicza kąt w radianach na stopnie.
        
        Argumenty:
        ----------
        kat - kąt podany w radianach | typ: float lub int
        
        Wyniki:
        -------
        katd - kąt przeliczony z argumentu 'kat' na stopnie dziesiętne | typ: float
        
        Dodatkowy opis:
        ---------------
        '''
        from math import pi
        katd = kat * 180 / pi
        return(katd)
    
    def deg2dms(self, kat):
        """
        Funkcja przelicza stopnie dziesiętne na stopnie, minuty i sekundy.

        Parameters
        ----------
        kat : TYPE
            stopnie dziesiętne | typ: float lub int

        Returns
        -------
        DMS:
        
        Lista złożona w kolejności ze stopni, minut i sekund | typ: lista
        W liście DMS wartości:
        
        indeks = 0: stopnie | typ: int
        indeks = 1: minuty  | typ: int
        indeks = 2: sekundy | typ: float
        
        """
        
        st = abs(kat)
        mn, sek =divmod(st * 3600, 60)
        st, mn = divmod(mn, 60)
        sek = round(sek, 5)
        mn = int(mn)
        st = int(st)
        DMS = [st, mn, sek]
        
        return(DMS)
        
    def xyz2plh(self, X, Y, Z):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.

        Parameters:
        ----------
        X : TYPE: FLOAT
            Współrzędna X w układzie ortokartezjańskim
        Y : TYPE: FLOAT
            Współrzędna Y w układzie ortokartezjańskim
        Z : TYPE: FLOAT
            Współrzędna Z w układzie ortokartezjańskim

        Returns
        -------
        phi : TYPE: FLOAT
            [stopnie dziesiętne] - szerokość geodezyjna
        lam : TYPE: FLOAT
            [stopnie dziesiętne] - długośc geodezyjna.
        hel : TYPE FLOAT
            [metry] - wysokość elipsoidalna
            
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec

        """
        
    
        
        r   = sqrt(X**2 + Y**2)           # promień
        phi_prev =atan(Z/(r*(1-self.ecc2)))
        phi_next = 0
        epsilon = 0.000001/206265
        while abs(phi_prev - phi_next)>epsilon:
            phi_prev = phi_next
            N    = self.a/ sqrt(1 - self.ecc2 * (sin(phi_prev))**2)
            hel  = (r/cos(phi_prev))- N
            phi_next = atan((Z/r) *(((1 - self.ecc2 * N/(N + hel))**(-1))))
        phi = phi_next
        lam = atan(Y/X)
        N = self.a/ sqrt(1 - self.ecc2 * (sin(phi))**2);

        return degrees(phi), degrees(lam), hel
        return phi, lam, hel
    
    def plh2xyz(self, phi, lam, hel):
        
        '''
        Funkcja przelicza ze współrzędnych krzywoliniowych na współrzędne prostokątne.
        
        Parameters:
        ----------
        
        phi - szerokość geograficzna punktu | typ: lista
        lam - długość geograficzna punktu   | typ: lista
        hel - wysokość punktu               | typ: float lub int

        Returns:
        -------
        X - współrzędna prostokątna X punktu | typ: float
        Y - współrzędna prostokątna Y punktu | typ: float
        Z - współrzędna prostokątna Z punktu | typ: float
        
        '''

        phi = self.deg2rad(phi)
        lam = self.deg2rad(lam)
        
        
        N =  self.a / (sqrt(1 - self.ecc2 * (sin(phi)) ** 2))
        
        X = (N + hel) * cos(phi) * cos(lam)
        Y = (N + hel) * cos(phi) * sin(lam)
        Z = (N * (1 - self.ecc2) + hel) * sin(phi)
        
        return(X, Y, Z)
    
    """ Przeliczenie na NEU"""
    def średnia(self, wartosci):
        """
        Funkcja liczy średnią wartość z elementów w liscie
        
        Parameters:
        ----------
        wartosci : [float] : lista wartosci
        
        Returns:
        --------
        srednia : [float] : średnia arytmetyczna elementów z listy 
        
        """
        suma = 0
        ilość = 0
        for wartosc in wartosci:
            suma += wartosc
            ilość += 1
        srednia = float(suma / ilość)
        return(srednia)
    
    def Rneu(self, phi, lam):
        """
        Funkcja, która, przyjmujac współrzedne krzywoliniowe utworzy macierz obrotu 
        potrzebną do przeliczenia współrzędnych do układu współrzędnych neu
    
        INPUT:
        ----------
        phi : [float] : wspołrzędna fi punktu początkowego układu lokalnego
        lam : [float] :wspołrzędna l punktu początkowego układu lokalnego
    
        OUTPUT:
        -------
        R : [array of float64] : macierz obrotu
    
        """
        N=[(-sin(phi) * cos(lam)), (-sin(phi) * sin(lam)), (cos(phi))]
        E=[(-sin(lam)), (cos(lam)),  (0)]
        U=[( cos(phi) * cos(lam)), ( cos(phi) * sin(lam)), (sin(phi))]
        R=np.transpose(np.array([N,E,U]))
        return (R, N, E, U)
    
    def NEU(self, R,v):
        """
        Funckja obliczająca wektor w układzie neu
    
        Parameters:
        -----------
        R : R : [array of float64] : macierz obrotu
        v : [array of float64] : wektor w układzie XYZ
        
        Returns:
        -------
        delta_neu : [array of float64] : współrzedne topocentryczne (North , East (E), Up (U))
    
        """
        NEU=np.zeros(v.shape)
        for a in range(v.shape[0]):
            for b in range(3):
                for c in range(3):
                    NEU[a,c]+=v[a,b]*R[c,b]
        return (NEU)
    
    """ Przeliczenie na xy2000"""
    def fl2xy2000(self, phi, lam):
        """
        Funkcja przelicza współrzędne geodezyjne na współrzędne układu 2000.

        Parameters
        ----------
        phi : TYPE : [float] : Szerokość geodezyjna [stopnie]
        lam : TYPE : [float] : Długość geodezyjna [stopnie]

        Returns
        -------
        x2000 : TYPE : [float] : współrzędna X w układzie 2000 [m]
        y2000 : TYPE : [float] : współrzędna Y w układzie 2000 [m]

        """
        if (lam > 13.5 and lam) < 16.5:
            strefa = 5
            lam0 = 15
        elif (lam > 16.5 and lam < 19.5):
            strefa = 6
            lam0 = 18
        elif (lam > 19.5 and lam < 22.5):
            strefa = 7
            lam0 = 21
        elif (lam > 22.5 and lam < 25.5):
            strefa = 8
            lam0 = 24
        phi = self.deg2rad(phi)
        lam = self.deg2rad(lam)
        lam0 = self.deg2rad(lam0)
        b2 = (self.a**2)*(1 - self.ecc2)
        ep2 = ((self.a**2 - b2))/(b2)
        t = atan(phi)
        n2 = ep2 * ((cos(phi))**2)
        N =  self.a / (sqrt(1 - self.ecc2 * (sin(phi)) ** 2))
        A0 = 1-(self.ecc2/4) - ((3*(self.ecc2**2))/64)-((5*(self.ecc2**3))/256)
        A2 = (3/8)*(self.ecc2+(self.ecc2**2/4)+((15*(self.ecc2**3))/128))
        A4 = (15/256)*((self.ecc2**2)+((3*(self.ecc2**3))/4))
        A6 = (35*(self.ecc2**3))/3072
        sigma = self.a * (A0*(phi) - A2 * sin(2*phi) + A4 * sin(4*phi) - A6 * sin(6 * phi))
        dlam = lam - lam0
        xgk = sigma + ((dlam**2)/2) * N * sin(phi) * cos(phi) * (1 + ((dlam**2)/12) * ((cos(phi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((dlam**4)/360) * ((cos(phi))**4) * (61-58 * (t**2) + t**4 + 270 * n2 - 330 * n2 * (t**2)))
        ygk = dlam * N * cos(phi) * (1 + ((dlam**2)/6) * ((cos(phi))**2) * (1 - t**2 + n2) + ((dlam**4)/120) * ((cos(phi))**4) * (5 - 18 * (t**2) + t**4 + 14 * n2 - 58 * n2 * (t**2)))
        m = 0.999923 #skala PL-2000
        x2000 = xgk * m
        y2000 = ygk * m + (strefa * 1000000) + 500000
        return x2000, y2000
    
    """ Przeliczenie na xy1992"""
    def fl2xy1992(self, phi, lam):
        """
        Funkcja przelicza współrzędne geodezyjne na współrzędne układu 1992.

        Parameters
        ----------
        phi : TYPE : [float] : Szerokość geodezyjna [stopnie]
        lam : TYPE : [float] : Długość geodezyjna [stopnie]

        Returns
        -------
        x1992 : TYPE : [float] : współrzędna X w układzie 1992 [m]
        y1992 : TYPE : [float] : współrzędna Y w układzie 1992 [m]

        """
        
        phi = self.deg2rad(phi)
        lam = self.deg2rad(lam)
        lam0 = self.deg2rad(19)
        b2 = (self.a**2)*(1 - self.ecc2)
        ep2 = ((self.a**2 - b2))/(b2)
        t = tan(phi)
        n2 = ep2 * ((cos(phi))**2)
        N =  self.a / (sqrt(1 - self.ecc2 * (sin(phi)) ** 2))
        A0 = 1-(self.ecc2/4) - ((3*(self.ecc2**2))/64)-((5*(self.ecc2**3))/256)
        A2 = (3/8)*(self.ecc2+(self.ecc2**2/4)+((15*(self.ecc2**3))/128))
        A4 = (15/256)*((self.ecc2**2)+((3*(self.ecc2**3))/4))
        A6 = (35*(self.ecc2**3))/3072
        dlam = lam - lam0
        sigma = self.a * (A0*(phi) - A2 * sin(2*phi) + A4 * sin(4*phi) - A6 * sin(6 * phi))
        xgk = sigma + ((dlam**2)/2) * N * sin(phi) * cos(phi) * (1 + ((dlam**2)/12) * ((cos(phi))**2) * (5 - t**2 + 9 * n2 + 4 * (n2**2)) + ((dlam**4)/360) * ((cos(phi))**4) * (61-58 * (t**2) + t**4 + 270 * n2 - 330 * n2 * (t**2)))
        ygk = dlam * N * cos(phi) * (1 + ((dlam**2)/6) * ((cos(phi))**2) * (1 - t**2 + n2) + ((dlam**4)/120) * ((cos(phi))**4) * (5 - 18 * (t**2) + t**4 + 14 * n2 - 58 * n2 * (t**2)))
        m = 0.9993
        x1992 = xgk * m - 5300000
        y1992 = ygk * m + 500000
        return x1992, y1992
    """ Wyznaczenie kąta Azymutu """
    def Azymut(self, N, E):
        """   
        Funkcja wyznacza kąt azymutu na podstawie współrzędnych topocentrycznych
        
        Parameters
        -------
        N  : [float] : wpółrzedna topocentryczna N (north) [m]
        E  : [float] : wpółrzedna topocentryczna E (east) [m]
        
        Returns
        -------
        Az : [float] : azymut [rad]
        
        """  
        Az = atan2(E, N)
        
        return Az
    
    """ Wyznaczenie kąta elewacji """
    def Kąt_Elewacji(self, N, E, U):
        """   
        Funkcja wyznacza kąt kąt elewacji na podstawie współrzędnych topocentrycznych
        
        Parameters
        -------
        N  : [float] : wpółrzedna topocentryczna N (north) [m]
        E  : [float] : wpółrzedna topocentryczna E (east) [m]
        U  : [float] : wpółrzedna topocentryczna U (up) [m] 
       
        Returns
        -------
        el : [float] : kąt elewacji [rad]
        
        """  
        el = acos(U/sqrt(N**2 + E**2 + U**2))
        
        return el
    
    """ Wyznaczenie odległości 2D """
    def odl2D(self,X_sr, Y_sr, X, Y):
        """   
        Funkcja wyznacza odległość na płaszczyźnie na podstawie współrzędnych płaskich prostokątnych początek przyjmujęmy jako średnia X oraz średnia Y
        
        Parameters
        -------
        X_sr  : [float] : średnia współrzędna X wszystkich punktów [m]
        Y_sr  : [float] : średnia współrzędna Y wszystkich punktów [m]
        X  : [float] : współrzędna X punktu  [m]
        Y  : [float] : współrzędna Y punktu  [m]
       
        Returns
        -------
        odległość2D : [float] : odległość na płaszczyźnie (2D) [m]
        
        """     
        odległość2D = []
        for X, Y in zip(X, Y):
            d = sqrt((X_sr - X)**2 + (Y_sr - Y)**2)
            odległość2D.append(d)
        
        return odległość2D
    
    """ Wyznaczenie odległości 3D """
    def odl3D(self,N_sr, E_sr, U_sr, N, E, U):
            
        '''
        Funkcja liczy długość wektora przyjmując elementy z listy jako współrzędne końcowe wektora, a współrzędne średnie jako współrzędne początkowe.

        Argumenty:
        ----------
        N_sr - współrzędna N średnia                     | typ: float lub int
        E_sr - współrzędna E średnia                     | typ: float lub int
        U_sr - współrzędna U średnia                     | typ: float lub int
        N - lista współrzędnych N dla kolejnych punktów | typ: lista wartości float lub int
        E - lista współrzędnych E dla kolejnych punktów | typ: lista wartości float lub int
        U - lista współrzędnych U dla kolejnych punktów | typ: lista wartości float lub int

        Wyniki:
        -------
        Dlugosc - lista odległośći punktów od X, Y, Z średnich | typ: lista wartości float
    
        Dodatkowy opis:
        ---------------
    
        ---
        '''
        odległość = []
    
        for N, E, U in zip(N, E, U):
            odl = sqrt((N_sr - N) ** 2 + (E_sr - E) ** 2 + (U_sr - U) ** 2)
            odległość.append(odl)
        return(odległość)
    
    
    
if __name__ == "__main__":
    #utworzenie obiektu
    geo = Transformacje(model = "grs80")
    plik = "wsp_inp.txt"
    tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)
    dane = np.ones((12,3))
    for i, n in enumerate(tablica):
        phi, lam, hel = geo.xyz2plh(tablica[i,0], tablica[i,1], tablica[i,2])
        dane[i,:] = [phi, lam, hel]
        
    print(dane) #phi, lam, hel
    
    dane2 = np.ones((12,3))
    for i, n in enumerate(dane):
        X, Y, Z = geo.plh2xyz(dane[i,0], dane[i,1], dane[i,2])
        dane2[i,:] = [X, Y, Z]
        
    print(dane2) # X, Y, Z
    
    phi_sr = geo.średnia(dane[:,0])
    lam_sr = geo.średnia(dane[:,1]) 
    X_sr = geo.średnia(dane2[:,0])
    Y_sr = geo.średnia(dane2[:,1])
    Z_sr = geo.średnia(dane2[:,2])
    [phi_sr, lam_sr, hel_sr] = geo.xyz2plh(X_sr, Y_sr, Z_sr)
    R, N, E, U = geo.Rneu(phi_sr,lam_sr)
    
    ii=0
    v=np.array(np.zeros((dane2.shape[0],3)))
    for ii in range(0, dane2.shape[0]):
        v[ii,0]=X_sr - dane2[ii,0]
        v[ii,1]=Y_sr - dane2[ii,1]
        v[ii,2]=Z_sr - dane2[ii,2]
        ii += 1

    neu=geo.NEU(R, v)
    print(neu) # N, E, U
        
    dane3 = np.ones((12,2))
    for i, n in enumerate(dane):
        x2000, y2000 = geo.fl2xy2000(dane[i,0], dane[i,1])
        dane3[i,:] = [x2000, y2000]
        
    print(dane3) # x2000, y2000
    
    dane4 = np.ones((12,2))
    for i, n in enumerate(dane):
        x1992, y1992 = geo.fl2xy1992(dane[i,0], dane[i,1])
        dane4[i,:] = [x1992, y1992]
        
    print(dane4) # x1992, y1992
    
    
    dane5 = np.ones((12,1))
    for i, n in enumerate(neu):
        Az = geo.Azymut(neu[i,0], neu[i,1])
        dane5[i,:] = [Az]
        
    print(dane5) # Azymut
    
    dane6 = np.ones((12,1))
    for i, n in enumerate(neu):
        el = geo.Kąt_Elewacji(neu[i,0], neu[i,1], neu[i,2])
        dane6[i,:] = [el]
        
    print(dane6) # Kąt elewacji
    
    
    odległość2D = geo.odl2D(X_sr, Y_sr, dane2[:,0], dane2[:,1])
    print(odległość2D)
    c = np.array([odległość2D])
    d = c.T
    print(d) #odległość 2D
    
    N_sr = geo.średnia(neu[:,0])
    E_sr = geo.średnia(neu[:,1])
    U_sr = geo.średnia(neu[:,2])
    odległość = geo.odl3D(N_sr, E_sr, U_sr, neu[:,0], neu[:,1], neu[:,2])
    a = np.array([odległość])
    b = a.T
    print(b) #odległość 3D

    DANE = np.hstack((dane, dane2, neu, dane3, dane4, dane5, dane6, d, b))
    print(DANE)
    
    
    np.savetxt("wsp_out.txt", DANE, delimiter='  ', fmt = ['%10.8f', '%10.8f', '%10.5f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.3f', '%10.3f', '%10.3f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], header = 'Konwersja współrzednych geodezyjnych \\ Michał Karaban', comments=' phi [st]     lambda[st]     hel[m]          X[m]              Y[m]              Z[m]          N[m]         E[m]         U[m]         x2000[m]        y2000[m]    x1992[m]    y1992[m]      Az[rad]     El[rad]     Odl2D[m]   Odl3D[m] \n ' )
    
