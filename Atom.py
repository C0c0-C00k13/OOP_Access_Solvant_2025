import lecturePDB

# Section 1 : Creation d'objet Atome
class Calc_Atom:
    
    """Ceci est la classe atome.
    ------
    Attributes
    atom_fullname : str
    atom_type : str
    x : float
    y : float
    z : float
    ------
    Methodes
    __init___ : intialise un atome.
    calcul_distance : renvoie la distance entre un point entre 2 atomes.
    __str__ : affiche la valeur des classes.  
    """
    
    def __init__(self, atom_fullname, atom_type, coord):
        self.atom_fullname = atom_fullname
        self.atom_type = atom_type
        self.x = coord[0]
        self.y = coord[1]
        self.z = coord[2]

    def calcul_distance(self, another_atome):
        help = "Methodes pemettant de calculer la distance entre 2 atomes"
        if isinstance(another_atome, Calc_Atom):
            return pow((pow((self.x - another_atome.x), 2) + pow((self.y - another_atome.y),2) + pow((self.z - another_atome.z),2)), 0.5)

    def __str__(self):
        """Redéfinition du comportement avec print()."""
        return f"Atome {self.atom_fullname}; type : {self.atom_type}; coords[{self.x},{self.y},{self.z}]"


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

    # Calcule de la distance entre 2 atomes
    print(f'Distance de atome 1 à atome 2 : {atome1.calcul_distance(atome2)}')
    
    # Utilisation de 'print' sur un abjet de classe Atome utilisation de la methode __str__
    print(atome1)
    