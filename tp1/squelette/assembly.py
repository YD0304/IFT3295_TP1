"""IFT3295 - TP1 - Assemblage de fragments (Ex. 2.2 et 2.3).

Ce fichier est A COMPLETER. Regles:
* Completez uniquement le corps des fonctions ci-dessous.
* Ne modifiez ni les noms, ni les signatures, ni les types.
* Vous pouvez utiliser NetworkX pour representer et manipuler le graphe.
* Typez tous vos parametres et retours (types simples: int, str, list, ...).
* Le programme principal est ``main.py`` (fourni): il appelle les fonctions
  de ce module, il ne faut donc rien executer ici directement.
"""

import networkx as nx

from overlap import chevauchement_maximal

SEUIL_DEFAUT: int = 80

# Les noeuds du graphe sont les entiers 0..n-1 (indices des reads); chaque
# arete porte l'attribut entier ``poids`` (= score du chevauchement).


def graphe_chevauchements(
    scores: list[list[int]], seuil: int = SEUIL_DEFAUT
) -> nx.DiGraph:
    """Construit le graphe oriente des chevauchements pertinents.

    Args:
        scores (list[list[int]]): Matrice carree des scores de chevauchement.
        seuil (int): Score minimal requis pour conserver une arete.

    Returns:
        nx.DiGraph: Graphe dont les noeuds sont les indices ``0..n-1``.
        Une arete ``(i, j)`` existe si ``scores[i][j] >= seuil``; son attribut
        ``poids`` vaut le score correspondant.
    """
    n = len(scores)
    graphe = nx.DiGraph() #Directed graph
    graphe.add_nodes_from(range(n))
    for i in range(n):
        for j in range(n):
            if i != j and scores[i][j] >= seuil: #seuil et pas diagonale
                graphe.add_edge(i, j, poids=scores[i][j])
    return graphe

def reduction_transitive(graphe: nx.DiGraph) -> nx.DiGraph:
    """Calcule la reduction transitive du graphe de chevauchement.

    Supprime une arete ``(u, v)`` lorsqu'un autre chemin simple de ``u`` a
    ``v`` subsiste. Pour un cycle de deux noeuds, conserve uniquement l'arete
    de plus grand poids.

    Args:
        graphe (nx.DiGraph): Graphe oriente dont les aretes portent un attribut
            entier ``poids``.

    Returns:
        nx.DiGraph: Copie reduite du graphe d'entree. Le graphe d'entree n'est
        pas modifie.
    """
    graphe_reduit = graphe.copy()

    # Supprimer les aretes de poids plus faible dans les cycles de deux noeuds
    # Traiter les paires de noeuds bidirectionnelles (cycles de taille 2)
    # On parcourt une copie des aretes pour ne pas modifier l'iterable
    for u, v in list(graphe_reduit.edges()):
        if graphe_reduit.has_edge(u, v) and graphe_reduit.has_edge(v, u):
            poids_uv = graphe_reduit[u][v]["poids"]
            poids_vu = graphe_reduit[v][u]["poids"]
            if poids_uv < poids_vu:
                graphe_reduit.remove_edge(u, v)
            elif poids_vu < poids_uv:
                graphe_reduit.remove_edge(v, u)
            else:
                # si egaux, garder une seule orientation: conserver (min,max)
                a, b = sorted((u, v))
                graphe_reduit.remove_edge(b, a)

    # Supprimer les aretes transitives: pour chaque arete (u,v), si un autre
    # chemin simple de u a v existe, on peut supprimer l'arete.
    for u, v in list(graphe_reduit.edges()):
        # l'arete pourrait avoir ete supprimee lors du traitement precedent
        if not graphe_reduit.has_edge(u, v):
            continue
        attributs = dict(graphe_reduit[u][v])
        graphe_reduit.remove_edge(u, v)
        if nx.has_path(graphe_reduit, u, v):
            # il existe un autre chemin, on laisse l'arete supprimee
            continue
        # aucun autre chemin: remettre l'arete
        graphe_reduit.add_edge(u, v, **attributs)

    return graphe_reduit


def ordre_assemblage(graphe: nx.DiGraph) -> list[int]:
    """Trouve l'ordre d'assemblage des reads dans le graphe reduit.

    Args:
        graphe (nx.DiGraph): Graphe reduit produit par
            ``reduction_transitive``.

    Returns:
        list[int]: Indices des reads dans un chemin du graphe qui visite
        chaque noeud exactement une fois. L'enonce garantit l'existence de
        ce chemin apres reduction.
    """
    #Edge case: graphe vide
    if len(graphe) == 0:
        return []

    #noed de depart: degre entrant = 0 (il y a exactement un tel noeud)
    depart: int | None = None
    for noeud in graphe:
        if graphe.in_degree(noeud) == 0:
            depart = noeud
            break

    ordre = [depart]
    courant = depart
    while len(ordre) < len(graphe):
        successeurs = list(graphe.successors(courant))
        courant = successeurs[0]
        ordre.append(courant)

    return ordre


def sequence_finale(reads: list[str], ordre: list[int]) -> tuple[str, list[int]]:
    """Assemble les reads en une sequence de fragment genomique.

    Args:
        reads (list[str]): Sequences des reads a assembler.
        ordre (list[int]): Indices des reads dans l'ordre d'assemblage.

    Returns:
        tuple[str, list[int]]: Tuple contenant la sequence assemblee et les
        longueurs des chevauchements entre les paires de reads consecutifs,
        dans l'ordre.
    """

   # Edge case: ordre vide
    if not ordre:
        return "", []

    # Assemble les reads selon l'ordre donne
    fragment = reads[ordre[0]]
    longueurs: list[int] = []

    # On parcourt les indices de l'ordre pour assembler les reads
    for k in range(len(ordre) - 1):
        precedent = ordre[k]
        suivant = ordre[k + 1]
        # On calcule le chevauchement maximal entre les deux reads
        _, _, _, chevauchement = chevauchement_maximal(
            reads[precedent], reads[suivant]
        )
        longueurs.append(chevauchement)
        fragment += reads[suivant][chevauchement:] #On ajoute la partie non chevauchante du read suivant

    return fragment, longueurs