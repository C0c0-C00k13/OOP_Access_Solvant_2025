"""Script de la classe point_atome.
"""

import Atom

class point_atome:   
    """Génère un objet de la classe point_atome.
    
    Attributs
    ------
    atom_center : Calc_Atom
        Atome auquel le point est relié
    x_pt : float
        Coordonées de l'atome sur l'axeX
    y_pt : float
        Coordonées de l'atome sur l'axY
    z_pt : float
        Coordonées de l'atome sur l'axeZ
    
    Methodes
    ------
    __init___ : intialise un point.
    calcul_distance : renvoie la distance entre un point entre un atome.
    __str__ : affiche les attributs de classe du point.  
    """
    def __init__(self, atom_center, x_pt, y_pt, z_pt):
        """Initalise un objet de la classe point_atome.
        
        Parameters
        ------
        atom_center: Calc_Atom
            Atome auquel le point est relié
        x_pt : float
            Coordonées de l'atome sur l'axeX
        y_pt : float
            Coordonées de l'atome sur l'axY
        z_pt : float
            Coordonées de l'atome sur l'axeZ
            
        Returns
        ------
        point_atome
            Objet de type point_atome
        """
        self.atom_center = atom_center 
        self.x_pt = x_pt
        self.y_pt = y_pt
        self.z_pt = z_pt

    def __str__(self):
        """
        """
        return f"Atome {print(self.atom_center)},\t Point coords[{self.x_pt},{self.y_pt},{self.z_pt}]"

    def calcul_distance(self, atome):
        """
        """
        help = "Methodes pemettant de calculer la distance entre 2 atomes"
        if isinstance(atome, Calc_Atom):
            return pow((pow((self.x_pt - atome.x), 2) + pow((self.y_pt - atome.y),2) + pow((self.z_pt - atome.z),2)), 0.5)

if __name__ == "__main__":
    protein1 = ("2C8Q","./Data/insuline.pdb")
    list_atome = PDBRetrieve_Atoms(protein1[0],protein1[1])
    
    new_liste_atom = []
    for grp_atom in list_atome:
        for atome in grp_atom:
            # print(f"Coordonnées : {atome.coord}; Full name : {atome.fullname}; Element : {atome.element}")
            new_atom = Calc_Atom(atome.fullname, atome.element, atome.coord)
            new_liste_atom.append(new_atom)
    new_liste_atom[0]
    atome1, atome2 = new_liste_atom[0], new_liste_atom[1]
    point1 = point_atome(atome1,1.5,1.6,0.2)
    print(point1)
    print(point1.calcul_distance(atome2))