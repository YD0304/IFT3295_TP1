"""IFT3295 - TP1 - Algorithme de chevauchement de sequences (Q1.5, Ex. 2.1).

Ce fichier est A COMPLETER. Regles:
* Completez uniquement le corps des fonctions ci-dessous.
* Ne modifiez ni les noms, ni les signatures, ni les types, ni les constantes.
* Typez tous vos parametres et retours (types simples: int, str, list, ...).
* Le programme principal est ``main.py`` (fourni): il appelle les fonctions
  de ce module, il ne faut donc rien executer ici directement.
"""

MATCH: int = 4
MISMATCH: int = -4
INDEL: int = -8


def chevauchement_maximal(x: str, y: str) -> tuple[int, str, str, int]:
    """Calcule le meilleur chevauchement ordonne entre deux sequences.

    Args:
        x (str): Premiere sequence a aligner.
        y (str): Deuxieme sequence a aligner.

    Returns:
        tuple[int, str, str, int]: Score maximal, deux lignes alignees et
        longueur du chevauchement.

    Examples:
        >>> chevauchement_maximal("ACCA", "CACGC")
        (8, 'CA', 'CA', 2)
        >>> chevauchement_maximal("CACGC", "ACCA")
        (4, 'ACGC', 'AC-C', 4)
    """
    n: int = len(x)
    m: int = len(y)

    # V[i][0] = 0 : le debut de x peut etre ignore gratuitement.
    # V[0][j] = j * INDEL : le debut de y doit etre aligne (pas de saut gratuit).
    v: list[list[int]] = [[0] * (m + 1) for _ in range(n + 1)]
    # Pointeurs: 0 = diagonale, 1 = haut (x contre '-'), 2 = gauche ('-' contre y)
    ptr: list[list[int]] = [[0] * (m + 1) for _ in range(n + 1)]

    for j in range(1, m + 1):
        v[0][j] = j * INDEL
        ptr[0][j] = 2

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s: int = MATCH if x[i - 1] == y[j - 1] else MISMATCH
            diag: int = v[i - 1][j - 1] + s
            haut: int = v[i - 1][j] + INDEL
            gauche: int = v[i][j - 1] + INDEL
            best: int = diag
            p: int = 0
            if haut > best:
                best = haut
                p = 1
            if gauche > best:
                best = gauche
                p = 2
            v[i][j] = best
            ptr[i][j] = p

    # La fin de x doit etre consommee; la fin de y est libre:
    # on prend le maximum de la derniere ligne.
    j_best: int = 0
    for j in range(m + 1):
        if v[n][j] > v[n][j_best]:
            j_best = j
    score: int = v[n][j_best]

    # Remontee jusqu'a atteindre la colonne 0 (debut de y aligne).
    i: int = n
    j = j_best
    a: list[str] = []
    b: list[str] = []
    while j > 0:
        p = ptr[i][j] if i > 0 else 2
        if p == 0:
            a.append(x[i - 1])
            b.append(y[j - 1])
            i -= 1
            j -= 1
        elif p == 1:
            a.append(x[i - 1])
            b.append("-")
            i -= 1
        else:
            a.append("-")
            b.append(y[j - 1])
            j -= 1

    ligne_x: str = "".join(reversed(a))
    ligne_y: str = "".join(reversed(b))
    
    return score, ligne_x, ligne_y, len(ligne_x)


def matrice_chevauchements(reads: list[str]) -> list[list[int]]:
    """Construit la matrice des scores de chevauchement de toutes les paires.

    Args:
        reads (list[str]): Sequences a comparer, par exemple les 20 reads de
            ``reads.fq``.

    Returns:
        list[list[int]]: Matrice carree dont l'element ``[i][j]`` est le score
        maximal de la paire ordonnee ``(reads[i], reads[j])``. La diagonale
        contient des zeros.
    """
    raise NotImplementedError  # TODO
