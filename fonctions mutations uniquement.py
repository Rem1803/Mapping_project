def recherche_avec_indels_et_sub(self, motif, k=1):
        """Cherche le motif avec au plus k erreurs (substitutions, insertions, délétions).

        Args:
            motif : séquence à chercher (str).
            k : nombre maximal d'erreurs autorisées sur l'ensemble du motif.

        Returns:
            dict {position de début dans le texte : nombre minimal d'erreurs}
        """
        resultats = {}
        deja_vu = {}  # (i, debut, fin) -> plus petit nombre d'erreurs déjà exploré pour cet état
        alphabet = [c for c in self.debut if c != "$"]

        def explorer(i, debut, fin, erreurs):
            # i : indice de la prochaine lettre du motif à aligner (on recule vers 0)
            # [debut, fin[ : intervalle des suffixes qui correspondent à la partie déjà alignée

            # Mémoïsation : inutile de ré-explorer un état déjà atteint avec moins d'erreurs.
            # Cela évite les alignements redondants et garantit l'arrêt de la récursion.
            cle = (i, debut, fin)
            if deja_vu.get(cle, k + 1) <= erreurs:
                return
            deja_vu[cle] = erreurs

            # Tout le motif est aligné : chaque suffixe de l'intervalle est une occurrence
            if i < 0:
                for r in range(debut, fin):
                    pos = self.suffix_array[r]
                    resultats[pos] = min(erreurs, resultats.get(pos, k + 1))
                return

            # Branche « lettre en trop dans le motif » : on saute motif[i],
            # le génome ne bouge pas (l'intervalle reste le même).
            if erreurs + 1 <= k:
                explorer(i - 1, debut, fin, erreurs + 1)

            for lettre in alphabet:
                nouveau_debut = self.debut[lettre] + self.apparition[lettre][debut]
                nouvelle_fin = self.debut[lettre] + self.apparition[lettre][fin]
                if nouveau_debut >= nouvelle_fin:
                    continue  # cette lettre ne prolonge aucun suffixe de l'intervalle

                # Match ou substitution : le motif ET le génome avancent d'une lettre
                cout = erreurs + (lettre != motif[i])
                if cout <= k:
                    explorer(i - 1, nouveau_debut, nouvelle_fin, cout)

                # Branche « lettre manquante dans le motif » : le génome avance,
                # mais on n'a pas consommé de lettre du motif (i reste le même).
                if erreurs + 1 <= k:
                    explorer(i, nouveau_debut, nouvelle_fin, erreurs + 1)

        explorer(len(motif) - 1, 0, self.n, 0)
        return dict(sorted(resultats.items()))